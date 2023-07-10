import pytest
import pandas as pd
import re

from csodiaq.identifier.fdrFunctions import score_library_to_query_matches, calculate_cosine_similarity_score, calculate_macc_score, add_decoy_label_to_score_df, determine_index_of_fdr_cutoff

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

def test__fdr_functions__add_decoy_label_to_score_df():
    isNotDecoy, isDecoy = 0, 1
    targetLibraryIdx, decoyLibraryIdx, queryIdx, score = 0, 1, 0, 0

    scoreData = [
        [targetLibraryIdx, queryIdx, score],
        [decoyLibraryIdx, queryIdx, score],
    ]
    scoreDf = pd.DataFrame(scoreData, columns=["libraryIdx","queryIdx","score"])
    expectedOutputData = [
        [targetLibraryIdx, queryIdx, score, isNotDecoy],
        [decoyLibraryIdx, queryIdx, score, isDecoy],
    ]
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=["libraryIdx","queryIdx","score", "isDecoy"])
    decoySet = set([decoyLibraryIdx])
    outputDf = add_decoy_label_to_score_df(decoySet, scoreDf)
    assert expectedOutputDf.equals(outputDf)

def test__fdr_functions__determine_index_of_fdr_cutoff():
    fdrCutoff = 0.01
    numberOfNonDecoys = int(1/fdrCutoff)
    decoys = [1, 1]
    isDecoySeries = pd.Series([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx - 1


def test__fdr_functions__determine_index_of_fdr_cutoff__first_decoy_appears_before_fdr_cutoff():
    numberOfNonDecoys = 1
    decoys = [1]
    isDecoySeries = pd.Series([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx - 1

def test__fdr_functions__determine_index_of_fdr_cutoff__throws_error_when_top_score_is_decoy():
    numberOfNonDecoys = 0
    decoys = [1]
    isDecoySeries = pd.Series([0] * numberOfNonDecoys + decoys)
    errorOutput = 'None of the library peptides were identified in the query spectra (highest score was a decoy).'
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)

