import pandas as pd
import numpy as np

from csodiaq.identification.matchingFunctions import (
    match_library_to_query_pooled_spectra,
    calculate_parts_per_million_relative_difference,
    is_within_tolerance,
    eliminate_low_count_matches,
    eliminate_matches_below_fdr_cutoff,
    calculate_ppm_offset_tolerance_using_mean_and_standard_deviation,
    calculate_ppm_offset_tolerance_using_tallest_bin_peak,
    filter_matches_by_ppm_offset_and_tolerance,
    identify_index_of_max_distance_to_noise_from_tallest_bin,
)


def test__matchingFunctions__calculate_parts_per_million_relative_difference():
    hundredToNineNinePpm = 10000
    thousandToNineNineNinePpm = 1000
    assert (
        calculate_parts_per_million_relative_difference(100, 99) == hundredToNineNinePpm
    )
    assert (
        calculate_parts_per_million_relative_difference(1000, 990)
        == hundredToNineNinePpm
    )
    assert (
        calculate_parts_per_million_relative_difference(1000, 999)
        == thousandToNineNineNinePpm
    )


def test__matchingFunctions__is_within_tolerance():
    tolerance = 1.0
    passValue1 = 0.9
    passValue2 = 1.0
    failValue = 1.1
    assert is_within_tolerance(passValue1, tolerance)
    assert is_within_tolerance(passValue2, tolerance)
    assert not is_within_tolerance(failValue, tolerance)


def find_target_mz(ppm, referenceMz):
    return referenceMz - (ppm * referenceMz) / 1e6


def test__matchingFunctions__match_library_to_query_pooled_spectra():
    ppmTolerance = 10.0
    ppmInsideTolerance = 10.0 - 1
    ppmOutsideTolerance = 10.0 + 1
    numPeaks = 5
    libraryIndex = 0
    matchingQueryIdx = 0
    mismatchingQueryIdx = 1
    genericIntensity = 100.0
    libPeaks = [(i * 100.0, genericIntensity, 0) for i in range(1, numPeaks + 1)]
    matchingMz = [
        find_target_mz(ppmInsideTolerance, libPeaks[i][0]) for i in range(numPeaks)
    ]
    mzMatchingPeaks = [
        (find_target_mz(ppmInsideTolerance, libPeaks[i][0]), genericIntensity, 0)
        for i in range(numPeaks)
    ]
    mzMismatchingPeaks = [
        (find_target_mz(ppmOutsideTolerance, libPeaks[i][0]), genericIntensity, 1)
        for i in range(numPeaks)
    ]
    mzPeaks = sorted(mzMatchingPeaks + mzMismatchingPeaks)
    expectedMatches = pd.DataFrame(
        [
            [
                libraryIndex,
                genericIntensity,
                matchingQueryIdx,
                genericIntensity,
                ppmInsideTolerance,
            ]
        ]
        * numPeaks,
        columns=[
            "libraryIdx",
            "libraryIntensity",
            "queryIdx",
            "queryIntensity",
            "ppmDifference",
        ],
    )
    expectedMatches.insert(loc=4, column="queryMz", value=matchingMz)
    matches = match_library_to_query_pooled_spectra(libPeaks, mzPeaks, ppmTolerance)
    matches["ppmDifference"] = [float(int(x)) for x in matches["ppmDifference"]]
    assert expectedMatches.equals(matches)


def test__matchingFunctions__eliminate_low_count_matches():
    highCountLibIdx = 0
    lowCountLibIdx = 1
    queryIdx = 0
    minNumMatches = 3
    highCountMatches = [[highCountLibIdx, queryIdx]] * minNumMatches
    lowCountMatches = [[lowCountLibIdx, queryIdx]] * (minNumMatches - 1)
    columns = ["libraryIdx", "queryIdx"]
    input = pd.DataFrame(highCountMatches + lowCountMatches, columns=columns)
    expectedOutput = pd.DataFrame(highCountMatches, columns=columns)
    output = eliminate_low_count_matches(input)
    assert expectedOutput.equals(output)


def test__matchingFunctions__eliminate_matches_below_fdr_cutoff():
    aboveCutoffLibIdx = 0
    belowCutoffLibIdx = 1
    queryIdx = 0
    minNumMatches = 3
    aboveCutoffMatches = [[aboveCutoffLibIdx, queryIdx]] * minNumMatches
    belowCutoffMatches = [[belowCutoffLibIdx, queryIdx]] * minNumMatches
    columns = ["libraryIdx", "queryIdx"]
    input = pd.DataFrame(aboveCutoffMatches + belowCutoffMatches, columns=columns)
    groupsAboveCutoff = [(aboveCutoffLibIdx, queryIdx)]
    expectedOutput = pd.DataFrame(aboveCutoffMatches, columns=columns)
    output = eliminate_matches_below_fdr_cutoff(input, groupsAboveCutoff)
    assert expectedOutput.equals(output)


def test__score_functions__calculate_ppm_offset_tolerance_using_mean_and_standard_deviation():
    mean = 10
    stdDev1 = 2
    stdDev2 = 4
    numbers = [7, 8, 9, 10, 11, 12, 13]
    (
        offset,
        tolerance,
    ) = calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(numbers, 1)
    assert offset == mean
    assert tolerance == stdDev1

    (
        offset,
        tolerance,
    ) = calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(numbers, 2)
    assert offset == mean
    assert tolerance == stdDev2


def test__score_functions__calculate_ppm_offset_tolerance_using_tallest_bin_peak():
    numBins = 200
    tallestBin = 50
    tallestBinQuantity = 100
    mediumLeftNeighboringBin = tallestBin - 1
    mediumRightNeighboringBin = tallestBin + 1
    numbers = list(range(-numBins // 2, numBins // 2))
    numbers += [tallestBin] * (tallestBinQuantity - 1)
    numbers += [mediumLeftNeighboringBin] * (tallestBinQuantity // 2 - 1)
    numbers += [mediumRightNeighboringBin] * (tallestBinQuantity // 2 - 1)
    expectedOffset = tallestBin
    expectedTolerance = 2
    offset, tolerance = calculate_ppm_offset_tolerance_using_tallest_bin_peak(numbers)
    assert abs(offset - expectedOffset) < 0.5
    assert abs(tolerance - expectedTolerance) < 0.5


def test__score_functions__identify_index_of_max_distance_to_noise_from_tallest_bin():
    noisePeakHeight = 1
    mediumPeakHeight = 50
    tallestPeakHeight = 100
    numNoisePeaks = 11
    numMediumPeaksLeft = 3
    numMediumPeaksRight = 2
    peaks = (
        [noisePeakHeight] * numNoisePeaks
        + [mediumPeakHeight] * numMediumPeaksLeft
        + [tallestPeakHeight]
        + [mediumPeakHeight] * numMediumPeaksRight
        + [noisePeakHeight] * numNoisePeaks
    )
    tallestPeakIdx = numNoisePeaks + numMediumPeaksLeft
    expectedIdx = numNoisePeaks - 1
    idx = identify_index_of_max_distance_to_noise_from_tallest_bin(
        np.array(peaks), tallestPeakIdx
    )
    assert expectedIdx == idx


def test__score_functions__filter_matches_by_ppm_offset_and_tolerance():
    libIdx = 0
    queryIdx = 0
    ppmOffset = 6
    matches = [[libIdx, queryIdx, (i * 5) + ppmOffset] for i in range(-5, 6)]
    columns = [
        "libraryIdx",
        "queryIdx",
        "ppmDifference",
    ]
    input = pd.DataFrame(matches, columns=columns)
    ppmTolerance = 11
    expectedOutput = input.iloc[3:-3,].reset_index(drop=True)
    output = filter_matches_by_ppm_offset_and_tolerance(input, ppmOffset, ppmTolerance)
    assert expectedOutput.equals(output)
