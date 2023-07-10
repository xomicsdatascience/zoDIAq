import pytest
import pandas as pd
import re
import numpy as np

from csodiaq.identifier.fdrFunctions import score_library_to_query_matches, calculate_cosine_similarity_score, calculate_macc_score, identify_all_decoys, determine_index_of_fdr_cutoff, calculate_ppm_offset_tolerance, calculate_ppm_offset_tolerance_using_mean_and_standard_deviation, calculate_ppm_offset_tolerance_using_tallest_bin_peak

@pytest.fixture
def vectorA():
    return pd.Series([1, 5, 10])

@pytest.fixture
def vectorB(vectorA):
    return vectorA + 1

def test__fdr_functions__calculate_cosine_similarity_score(vectorA, vectorB):
    AB = sum([vectorA[i] * vectorB[i] for i in range(len(vectorA))])
    A = sum([vectorA[i]**2 for i in range(len(vectorA))])
    B = sum([vectorB[i]**2 for i in range(len(vectorA))])
    expectedScore = AB / (A**0.5 * B**0.5)
    score = calculate_cosine_similarity_score(vectorA, vectorB)
    assert score == expectedScore

def test__fdr_functions__calculate_macc_score(vectorA, vectorB):
    cosineSimilarity = calculate_cosine_similarity_score(vectorA, vectorB)
    expectedScore = (len(vectorA)**(1/5)) * cosineSimilarity
    score = calculate_macc_score(vectorA, vectorB)
    assert score == expectedScore

def test__fdr_functions__score_library_to_query_matches(vectorA, vectorB):
    libraryIdx, queryIdx, ppmDiff = 0, 0, 0
    matchesDf = pd.DataFrame(index=vectorA.index,columns=["libraryIdx","libraryIntensity","queryIdx","queryIntensity","ppmDifference"])
    matchesDf["libraryIdx"] = [libraryIdx for i in vectorA.index]
    matchesDf["libraryIntensity"] = vectorA
    matchesDf["queryIdx"] = [queryIdx for i in vectorA.index]
    matchesDf["queryIntensity"] = vectorB
    matchesDf["ppmDifference"] = [ppmDiff for i in vectorA.index]
    maccScore = calculate_macc_score(vectorA, vectorB)
    expectedOutputDf  = pd.DataFrame(data=[[libraryIdx, queryIdx, maccScore]], columns=["libraryIdx","queryIdx","score"])
    outputDf = score_library_to_query_matches(matchesDf)
    assert expectedOutputDf.equals(outputDf)

def test__fdr_functions__identify_all_decoys():
    isNotDecoy, isDecoy = 0, 1
    targetLibraryIdx, decoyLibraryIdx, queryIdx, score = 0, 1, 0, 0

    scoreData = [
        [targetLibraryIdx, queryIdx, score],
        [decoyLibraryIdx, queryIdx, score],
    ]
    scoreDf = pd.DataFrame(scoreData, columns=["libraryIdx","queryIdx","score"])
    expectedOutput = np.array([
        isNotDecoy,
        isDecoy
    ])
    decoySet = set([decoyLibraryIdx])
    output = identify_all_decoys(decoySet, scoreDf)
    assert np.array_equal(output, expectedOutput)

def test__fdr_functions__determine_index_of_fdr_cutoff():
    fdrCutoff = 0.01
    numberOfNonDecoys = int(1/fdrCutoff)
    decoys = [1, 1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx


def test__fdr_functions__determine_index_of_fdr_cutoff__first_decoy_appears_before_fdr_cutoff():
    numberOfNonDecoys = 1
    decoys = [1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx

def test__fdr_functions__determine_index_of_fdr_cutoff__throws_error_when_top_score_is_decoy():
    numberOfNonDecoys = 0
    decoys = [1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    errorOutput = 'None of the library peptides were identified in the query spectra (highest score was a decoy).'
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)

def test__fdr_functions__calculate_ppm_offset_tolerance(): pass

def test__fdr_functions__calculate_ppm_offset_tolerance_using_mean_and_standard_deviation():
    mean = 10
    stdDev1 = 2
    stdDev2 = 4
    numbers = [7, 8, 9, 10, 11, 12, 13]
    offset, tolerance = calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(numbers, 1)
    assert offset == mean
    assert tolerance == stdDev1

    offset, tolerance = calculate_ppm_offset_tolerance_using_mean_and_standard_deviation(numbers, 2)
    assert offset == mean
    assert tolerance == stdDev2

def test__fdr_functions__calculate_ppm_offset_tolerance_using_tallest_bin_peak():
    numBins = 200
    tallestBin = 50
    tallestBinQuantity = 100
    mediumLeftNeighboringBin = tallestBin - 1
    mediumRightNeighboringBin = tallestBin + 1
    numbers = list(range(-numBins//2,numBins//2))
    numbers += [tallestBin] * (tallestBinQuantity - 1)
    numbers += [mediumLeftNeighboringBin] * (tallestBinQuantity // 2 - 1)
    numbers += [mediumRightNeighboringBin] * (tallestBinQuantity // 2 - 1)
    expectedOffset = tallestBin
    expectedTolerance = 2
    offset, tolerance = calculate_ppm_offset_tolerance_using_tallest_bin_peak(numbers)
    assert abs(offset - expectedOffset) < 0.5
    assert abs(tolerance - expectedTolerance) < 0.5
