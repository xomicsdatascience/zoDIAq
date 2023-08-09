import pandas as pd
import numpy as np
from numba import njit
from enum import Enum
import matplotlib.pyplot as plt

class Increment(Enum):
    NEITHER = 0
    LIBRARY = 1
    QUERY = 2


@njit
def calculate_parts_per_million_relative_difference(referenceMz, targetMz):
    return (referenceMz - targetMz) * (1e6) / referenceMz


@njit
def is_within_tolerance(ppm, tolerance):
    return abs(ppm) <= tolerance


@njit
def determine_smallest_peak_outside_ppm_tolerance(libMz, queryMz, ppmTolerance):
    ppm = calculate_parts_per_million_relative_difference(libMz, queryMz)
    if is_within_tolerance(ppm, ppmTolerance):
        return Increment.NEITHER
    elif libMz > queryMz:
        return Increment.QUERY
    else:
        return Increment.LIBRARY


@njit
def match_query_peak_to_all_succeeding_library_peaks_within_tolerance(
    baselineLibraryIdx, libraryPeaks, queryPeak, ppmTolerance
):
    tempLibraryIdx = baselineLibraryIdx + 0
    mzIdx, intensityIdx, tagIdx = 0, 1, 2
    data = []
    while tempLibraryIdx < len(libraryPeaks):
        ppm = calculate_parts_per_million_relative_difference(
            libraryPeaks[tempLibraryIdx][mzIdx], queryPeak[mzIdx]
        )
        if not is_within_tolerance(ppm, ppmTolerance):
            return data
        data.append(
            [
                libraryPeaks[tempLibraryIdx][tagIdx],
                libraryPeaks[tempLibraryIdx][intensityIdx],
                queryPeak[tagIdx],
                queryPeak[intensityIdx],
                queryPeak[mzIdx],
                ppm,
            ]
        )
        tempLibraryIdx += 1
    return data


@njit
def numba_enhanced_matching_of_library_to_query_pooled_spectra(
    libraryPeaks, queryPeaks, ppmTolerance
):
    baselineLibraryIdx, baselineQueryIdx = 0, 0
    mzIdx, intensityIdx, tagIdx = 0, 1, 2
    data = []
    while baselineLibraryIdx < len(libraryPeaks) and baselineQueryIdx < len(queryPeaks):
        incrementation = determine_smallest_peak_outside_ppm_tolerance(
            libraryPeaks[baselineLibraryIdx][mzIdx],
            queryPeaks[baselineQueryIdx][mzIdx],
            ppmTolerance,
        )
        if incrementation == Increment.LIBRARY:
            baselineLibraryIdx += 1
            continue
        if incrementation == Increment.QUERY:
            baselineQueryIdx += 1
            continue

        data.extend(
            match_query_peak_to_all_succeeding_library_peaks_within_tolerance(
                baselineLibraryIdx,
                libraryPeaks,
                queryPeaks[baselineQueryIdx],
                ppmTolerance,
            )
        )
        baselineQueryIdx += 1
    return data


def match_library_to_query_pooled_spectra(libraryPeaks, queryPeaks, ppmTolerance):
    """
    Peaks between library and query spectra are compared, and peaks with an m/z ppm tolerance
        are returned as matches.


    Parameters
    ----------
    libraryPeaks : list
        Each value in the list is a 3-length tuple with the following characteristics:
            m/z - the m/z value of the peak.
            intensity - the intensity value of the peak.
            tag - an int representing the library spectrum to which the peak belongs.
                This characteristics allows for pooling multiple library spectra
                together and being able to parse them out from each other after the
                analysis.

    queryPeaks : list
        Uses the same format as libraryPeaks, but the tag identifies query spectra
            instead.

    ppmTolerance : int
        The ppm tolerance allowed between library and query spectra peaks m/z values
            to be considered a match.

    Returns
    -------
    matchDf : pandas DataFrame
        A dataframe containing the library peak to query peak matches. Data includes
            the library/query tags of the peaks, their intensities, and the ppm
            difference between the two.
    """
    libraryArray = np.array(libraryPeaks)
    queryArray = np.array(queryPeaks)
    data = numba_enhanced_matching_of_library_to_query_pooled_spectra(
        libraryArray, queryArray, ppmTolerance
    )
    matchDf = pd.DataFrame(
        data,
        columns=[
            "libraryIdx",
            "libraryIntensity",
            "queryIdx",
            "queryIntensity",
            "queryMz",
            "ppmDifference",
        ],
    )
    matchDf[["libraryIdx", "queryIdx"]] = matchDf[["libraryIdx", "queryIdx"]].astype(
        int
    )
    return matchDf


def eliminate_low_count_matches(matches, minNumMatches=3):
    return (
        matches.groupby(["libraryIdx", "queryIdx"])
        .filter(lambda x: x.shape[0] >= minNumMatches)
        .reset_index(drop=True)
    )


def eliminate_matches_below_fdr_cutoff(matches, groupsAboveCutoff):
    return matches.groupby(["libraryIdx", "queryIdx"]).filter(
        lambda x: x.name in groupsAboveCutoff
    )

def create_ppm_histogram(ppms, offset, tolerance, histogramFile):
    binHeights, bins = create_standardized_histogram(ppms)
    barReductionForVisibility = 0.7
    binWidth = barReductionForVisibility * (bins[1] - bins[0])
    plt.clf()
    plt.bar(bins, binHeights, align="center", width=binWidth)
    plt.axvline(x=offset, color="black", linestyle="dashed", linewidth=4)
    plt.axvline(x=offset - tolerance, color="red", linestyle="dashed", linewidth=4)
    plt.axvline(x=offset + tolerance, color="red", linestyle="dashed", linewidth=4)
    plt.suptitle("offset: " + str(offset) + ", tolerance: " + str(tolerance))
    plt.savefig(histogramFile)


def filter_matches_by_ppm_offset_and_tolerance(matchDf, offset, tolerance):
    ppmLowerBound = offset - tolerance
    ppmUpperBound = offset + tolerance
    matchDf = matchDf[
        (matchDf["ppmDifference"] > ppmLowerBound)
        & (matchDf["ppmDifference"] < ppmUpperBound)
    ]
    return eliminate_low_count_matches(matchDf)

def calculate_ppm_offset_tolerance(ppms, numStandardDeviations):
    if numStandardDeviations:
        return calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(
            ppms, numStandardDeviations
        )
    else:
        return calculate_ppm_offset_tolerance_using_tallest_bin_peak(ppms)


def calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(
    ppms, numStandardDeviations
):
    return np.mean(ppms), np.std(ppms) * numStandardDeviations


def create_standardized_histogram(ppms, numBins=200):
    binHeights, bins = np.histogram(ppms, bins=numBins)
    bins = normalize_bin_position_from_left_to_center_of_each_bin(bins)
    return binHeights, bins


def calculate_ppm_offset_tolerance_using_tallest_bin_peak(ppms):
    binHeights, bins = create_standardized_histogram(ppms)
    tallestBinIdx = max(range(len(binHeights)), key=binHeights.__getitem__)
    nearestNoiseBinIdx = identify_index_of_max_distance_to_noise_from_tallest_bin(
        binHeights, tallestBinIdx
    )
    offset = bins[tallestBinIdx]
    tolerance = abs(offset - bins[nearestNoiseBinIdx])
    return offset, tolerance


def normalize_bin_position_from_left_to_center_of_each_bin(bins):
    return (bins[:-1] + bins[1:]) / 2


def identify_index_of_max_distance_to_noise_from_tallest_bin(binHeights, tallestBinIdx):
    averageNoise = np.mean(binHeights[:10] + binHeights[-10:])
    binsBeforeTallest = binHeights[:tallestBinIdx]
    beforeNearestNoiseBinIdx = (
        len(binsBeforeTallest) - np.argmin(binsBeforeTallest[::-1] >= averageNoise) - 1
    )
    binsAfterTallest = binHeights[tallestBinIdx:]
    afterNearestNoiseBinIdx = (
        np.argmin(binsAfterTallest >= averageNoise) + tallestBinIdx
    )
    nearestNoiseBinIdx = max(
        [beforeNearestNoiseBinIdx, afterNearestNoiseBinIdx],
        key=lambda x: abs(tallestBinIdx - x),
    )
    return nearestNoiseBinIdx