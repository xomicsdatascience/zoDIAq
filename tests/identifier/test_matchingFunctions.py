import pandas as pd

from csodiaq.identifier.matchingFunctions import match_library_to_query_pooled_spectra, calculate_parts_per_million_relative_difference, is_within_tolerance, eliminate_low_count_matches, eliminate_matches_below_fdr_cutoff


def test__matchingFunctions__calculate_parts_per_million_relative_difference():
    hundredToNineNinePpm = 10000
    thousandToNineNineNinePpm = 1000
    assert calculate_parts_per_million_relative_difference(100,99) == hundredToNineNinePpm
    assert calculate_parts_per_million_relative_difference(1000,990) == hundredToNineNinePpm
    assert calculate_parts_per_million_relative_difference(1000,999) == thousandToNineNineNinePpm

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
    libPeaks = [(i*100.0, genericIntensity, 0) for i in range(1, numPeaks+1)]
    matchingMz = [find_target_mz(ppmInsideTolerance, libPeaks[i][0]) for i in range(numPeaks)]
    mzMatchingPeaks = [(find_target_mz(ppmInsideTolerance, libPeaks[i][0]), genericIntensity, 0) for i in range(numPeaks)]
    mzMismatchingPeaks = [(find_target_mz(ppmOutsideTolerance, libPeaks[i][0]), genericIntensity, 1) for i in range(numPeaks)]
    mzPeaks = sorted(mzMatchingPeaks + mzMismatchingPeaks)
    expectedMatches = pd.DataFrame([[libraryIndex, genericIntensity, matchingQueryIdx, genericIntensity, ppmInsideTolerance]] * numPeaks, columns=["libraryIdx","libraryIntensity","queryIdx","queryIntensity","ppmDifference"])
    expectedMatches.insert(loc=4,
              column='queryMz',
              value=matchingMz)
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
    columns = ["libraryIdx","queryIdx"]
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
    columns = ["libraryIdx","queryIdx"]
    input = pd.DataFrame(aboveCutoffMatches + belowCutoffMatches, columns=columns)
    groupsAboveCutoff = [(aboveCutoffLibIdx, queryIdx)]
    expectedOutput = pd.DataFrame(aboveCutoffMatches, columns=columns)
    output = eliminate_matches_below_fdr_cutoff(input, groupsAboveCutoff)
    assert expectedOutput.equals(output)
