import pandas as pd

from csodiaq.identifier.matchingFunctions import find_all_library_matches_in_query_spectra, calculate_parts_per_million_relative_difference, is_within_tolerance


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

def test__matchingFunctions__find_all_library_matches_in_query_spectra():
    ppmTolerance = 10.0
    ppmInsideTolerance = 10.0 - 1
    ppmOutsideTolerance = 10.0 + 1
    numPeaks = 5
    libraryIndex = 0
    matchingQueryIdx = 0
    mismatchingQueryIdx = 1
    libPeaks = [(i*100.0, 100.0, 0) for i in range(1, numPeaks+1)]
    mzMatchingPeaks = [(find_target_mz(ppmInsideTolerance, libPeaks[i][0]), 100.0, 0) for i in range(numPeaks)]
    mzMismatchingPeaks = [(find_target_mz(ppmOutsideTolerance, libPeaks[i][0]), 100.0, 1) for i in range(numPeaks)]
    mzPeaks = sorted(mzMatchingPeaks + mzMismatchingPeaks)
    expectedMatches = pd.DataFrame([[libraryIndex, 100.0, matchingQueryIdx, 100.0, ppmInsideTolerance]] * numPeaks, columns=["libraryIdx","libraryIntensity","queryIdx","queryIntensity","ppmDifference"])
    matches = find_all_library_matches_in_query_spectra(libPeaks, mzPeaks, ppmTolerance)
    matches["ppmDifference"] = [float(int(x)) for x in matches["ppmDifference"]]
    assert expectedMatches.equals(matches)

