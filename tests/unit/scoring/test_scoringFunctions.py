import pytest
import pandas as pd
import re
import numpy as np
from math import isclose

from zodiaq.scoring.scoringFunctions import (
    score_library_to_query_matches,
    calculate_cosine_similarity_score,
    calculate_macc_score,
    determine_index_of_fdr_cutoff,
    calculate_fdr_rates_of_decoy_array,
)


@pytest.fixture
def vectorA():
    return pd.Series([1, 5, 10])


@pytest.fixture
def vectorB(vectorA):
    return vectorA + 1


def test__score_functions__calculate_cosine_similarity_score(vectorA, vectorB):
    AB = sum([vectorA[i] * vectorB[i] for i in range(len(vectorA))])
    A = sum([vectorA[i] ** 2 for i in range(len(vectorA))])
    B = sum([vectorB[i] ** 2 for i in range(len(vectorA))])
    expectedScore = AB / (A**0.5 * B**0.5)
    score = calculate_cosine_similarity_score(vectorA, vectorB)
    assert isclose(score, expectedScore)


def test__score_functions__calculate_macc_score(vectorA, vectorB):
    cosineSimilarity = calculate_cosine_similarity_score(vectorA, vectorB)
    expectedScore = (len(vectorA) ** (1 / 5)) * cosineSimilarity
    score = calculate_macc_score(len(vectorA), cosineSimilarity)
    assert score == expectedScore


def test__score_functions__score_library_to_query_matches(vectorA, vectorB):
    libraryIdx, queryIdx, ppmDiff = 1, 0, 0
    matchesDf = pd.DataFrame(
        index=vectorA.index,
        columns=[
            "libraryIdx",
            "libraryIntensity",
            "queryIdx",
            "queryIntensity",
        ],
    )
    matchesDf["libraryIdx"] = [libraryIdx for i in vectorA.index]
    matchesDf["libraryIntensity"] = vectorA
    matchesDf["queryIdx"] = [queryIdx for i in vectorA.index]
    matchesDf["queryIntensity"] = vectorB
    cosineScore = calculate_cosine_similarity_score(np.sqrt(vectorA), np.sqrt(vectorB))
    expectedOutputDf = pd.DataFrame(
        data=[[libraryIdx, queryIdx, cosineScore]],
        columns=["libraryIdx", "queryIdx", "cosineScore"],
    )
    outputDf = score_library_to_query_matches(matchesDf)
    assert expectedOutputDf.equals(outputDf)

    lowScoreMatchesDf = matchesDf.copy()
    lowScoreMatchesDf["libraryIdx"] = [libraryIdx - 1 for i in vectorA.index]
    reverseVectorA = pd.Series(list(vectorA)[::-1])
    lowScoreMatchesDf["libraryIntensity"] = reverseVectorA
    lowCosineScore = calculate_cosine_similarity_score(
        np.sqrt(reverseVectorA), np.sqrt(vectorB)
    )
    unsortedMatchesDf = pd.concat([lowScoreMatchesDf, matchesDf])
    expectedOutputDf = pd.DataFrame(
        data=[
            [libraryIdx, queryIdx, cosineScore],
            [libraryIdx - 1, queryIdx, lowCosineScore],
        ],
        columns=["libraryIdx", "queryIdx", "cosineScore"],
    )
    sortedOutputDf = score_library_to_query_matches(unsortedMatchesDf)
    assert expectedOutputDf.equals(sortedOutputDf)


def test__score_functions__calculate_fdr_rates_of_decoy_array():
    numberOfNonDecoys = 100
    decoys = [1, 1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    expectedFdrs = [0] * numberOfNonDecoys
    expectedFdrs.append(1 / (numberOfNonDecoys + 1))
    expectedFdrs.append(2 / (numberOfNonDecoys + 2))
    fdrs = calculate_fdr_rates_of_decoy_array(isDecoySeries)
    np.testing.assert_array_equal(expectedFdrs, fdrs)


def test__score_functions__determine_index_of_fdr_cutoff():
    fdrCutoff = 0.01
    numberOfNonDecoys = int(1 / fdrCutoff)
    decoys = [1, 1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx


def test__score_functions__determine_index_of_fdr_cutoff__first_decoy_appears_before_fdr_cutoff():
    numberOfNonDecoys = 1
    decoys = [1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    lastDecoyIdx = numberOfNonDecoys + len(decoys) - 1
    assert indexCutoff == lastDecoyIdx


def test__score_functions__determine_index_of_fdr_cutoff__throws_error_when_top_score_is_decoy():
    numberOfNonDecoys = 0
    decoys = [1]
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    errorOutput = "None of the library peptides were identified in the query spectra (highest score was a decoy)."
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)


def test__score_functions__determine_index_of_fdr_cutoff__returns_original_df_when_no_decoys_found():
    numberOfNonDecoys = 10
    decoys = []
    isDecoySeries = np.array([0] * numberOfNonDecoys + decoys)
    indexCutoff = determine_index_of_fdr_cutoff(isDecoySeries)
    assert indexCutoff == numberOfNonDecoys
