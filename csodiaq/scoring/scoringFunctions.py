import numpy as np
from scipy.spatial.distance import cosine
import matplotlib.pyplot as pyplot
from csodiaq.identification.matchingFunctions import eliminate_low_count_matches


def score_library_to_query_matches(matches):
    scoreDf = (
        matches.groupby(["libraryIdx", "queryIdx"])
        .apply(
            lambda x: (
                calculate_cosine_similarity_score(
                    x["libraryIntensity"], x["queryIntensity"]
                ),
                calculate_macc_score(x["libraryIntensity"], x["queryIntensity"]),
            )
        )
        .reset_index(name="scores")
    )
    scoreDf["cosineScore"], scoreDf["maccScore"] = zip(*scoreDf["scores"])
    del scoreDf["scores"]
    return scoreDf.sort_values("maccScore", ascending=False).reset_index(drop=True)


def calculate_cosine_similarity_score(vectorA, vectorB):
    return 1 - cosine(vectorA, vectorB)


def calculate_macc_score(vectorA, vectorB):
    return len(vectorA.index) ** (1 / 5) * calculate_cosine_similarity_score(
        vectorA, vectorB
    )

def determine_index_of_fdr_cutoff(isDecoyArray, fdrCutoff=1e-2):
    if isDecoyArray[0]:
        raise ValueError(
            "None of the library peptides were identified in the query spectra (highest score was a decoy)."
        )
    if 1 not in isDecoyArray:
        return len(isDecoyArray)
    fdrs = calculate_fdr_rates_of_decoy_array(isDecoyArray)
    lastDecoyIdx = np.argmax(fdrs > fdrCutoff)
    return lastDecoyIdx


def calculate_fdr_rates_of_decoy_array(isDecoyArray):
    isDecoyArrayIdx = np.arange(1, len(isDecoyArray) + 1)
    decoyCumSum = np.array(isDecoyArray).cumsum()
    return decoyCumSum / isDecoyArrayIdx


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


def create_ppm_histogram(ppms, offset, tolerance, histogramFile):
    binHeights, bins = create_standardized_histogram(ppms)
    barReductionForVisibility = 0.7
    binWidth = barReductionForVisibility * (bins[1] - bins[0])
    pyplot.clf()
    pyplot.bar(bins, binHeights, align="center", width=binWidth)
    pyplot.axvline(x=offset, color="black", linestyle="dashed", linewidth=4)
    pyplot.axvline(x=offset - tolerance, color="red", linestyle="dashed", linewidth=4)
    pyplot.axvline(x=offset + tolerance, color="red", linestyle="dashed", linewidth=4)
    pyplot.suptitle("offset: " + str(offset) + ", tolerance: " + str(tolerance))
    pyplot.savefig(histogramFile)


def filter_matches_by_ppm_offset_and_tolerance(matchDf, offset, tolerance):
    ppmLowerBound = offset - tolerance
    ppmUpperBound = offset + tolerance
    matchDf = matchDf[
        (matchDf["ppmDifference"] > ppmLowerBound)
        & (matchDf["ppmDifference"] < ppmUpperBound)
    ]
    return eliminate_low_count_matches(matchDf)
