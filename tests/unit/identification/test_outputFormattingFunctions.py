from csodiaq.identification.outputFormattingFunctions import (
    format_output_line,
    extract_metadata_from_match_and_score_dataframes,
    format_output_as_pandas_dataframe,
    identify_all_decoys,
)
import pandas as pd
import pytest
import numpy as np


@pytest.fixture
def identifierOutputData():
    return [
        "3",
        4,
        "testPeptide",
        "testProtein",
        0,
        0,
        1,
        8,
        "testIdentifiers",
        5,
        2,
        9,
        10,
        6,
        7,
        11,
    ]


def test__output_formatting_functions__format_output_line(identifierOutputData):
    libDict = {
        "peptide": "testPeptide",
        "proteinName": "testProtein",
        "isDecoy": 0,
        "precursorMz": 0,
        "precursorCharge": 1,
        "identification": "testIdentifiers",
        "peaks": list(range(2)),
    }
    queryDict = {
        "scan": "3",
        "precursorMz": 4,
        "peaksCount": 5,
        "CV": 6,
        "windowWidth": 7,
    }
    matchDict = {
        "cosineSimilarityScore": 8,
        "shared": 9,
        "ionCount": 10,
        "exclude_num": 11,
    }
    output = format_output_line(libDict, queryDict, matchDict)
    assert output == identifierOutputData


def test__output_formatting_functions__extract_metadata_from_match_and_score_dataframes():
    lib1Idx = 0
    lib2Idx = 1
    queryIdx = 0
    genericIntensity = 100.0
    genericPpmDifference = 10.0
    lib1PeakCount = 10
    lib1PeakCountAbovePrecursorMz = 3
    lib2PeakCount = 4
    precursorMz = 100.0
    abovePrecursorMz = precursorMz + 10
    belowPrecursorMz = precursorMz - 10
    lib1CosineScore = 1.0
    lib2CosineScore = 0.9
    lib1Score = 0.8
    lib2Score = 0.7
    lib1MatchAbovePrecursorMz = [
        [
            lib1Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            abovePrecursorMz,
            genericPpmDifference,
        ]
    ] * lib1PeakCountAbovePrecursorMz
    lib1MatchBelowPrecursorMz = [
        [
            lib1Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            belowPrecursorMz,
            genericPpmDifference,
        ]
    ] * (lib1PeakCount - lib1PeakCountAbovePrecursorMz)
    lib2Match = [
        [
            lib2Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            abovePrecursorMz,
            genericPpmDifference,
        ]
    ] * lib2PeakCount
    columns = [
        "libraryIdx",
        "libraryIntensity",
        "queryIdx",
        "queryIntensity",
        "queryMz",
        "ppmDifference",
    ]
    matchDf = pd.DataFrame(
        lib1MatchAbovePrecursorMz + lib1MatchBelowPrecursorMz + lib2Match,
        columns=columns,
    )

    scoreData = [
        [lib1Idx, queryIdx, lib1CosineScore],
        [lib2Idx, queryIdx, lib2CosineScore],
    ]
    scoreDf = pd.DataFrame(
        scoreData, columns=["libraryIdx", "queryIdx", "cosineScore"]
    )
    expectedOutput = {
        (lib1Idx, queryIdx): {
            "cosineSimilarityScore": lib1CosineScore,
            "shared": lib1PeakCount,
            "ionCount": genericIntensity * lib1PeakCountAbovePrecursorMz,
            "exclude_num": lib1PeakCount - lib1PeakCountAbovePrecursorMz,
        },
        (lib2Idx, queryIdx): {
            "cosineSimilarityScore": lib2CosineScore,
            "shared": lib2PeakCount,
            "ionCount": lib2PeakCount * genericIntensity,
            "exclude_num": 0,
        },
    }
    queryDict = {
        str(queryIdx): {
            "precursorMz": precursorMz,
        },
    }
    output = extract_metadata_from_match_and_score_dataframes(
        matchDf, scoreDf, queryDict
    )
    for idKey, metadataDict in expectedOutput.items():
        assert idKey in output
        for metadataType, metadataValue in metadataDict.items():
            assert metadataType in output[idKey]
            assert output[idKey][metadataType] == metadataValue


def test__output_formatting_functions__format_output_as_pandas_dataframe(
    identifierOutputData,
):
    expectedColumns = [
        "fileName",
        "scan",
        "MzEXP",
        "peptide",
        "protein",
        "isDecoy",
        "MzLIB",
        "zLIB",
        "cosine",
        "name",
        "Peak(Query)",
        "Peaks(Library)",
        "shared",
        "ionCount",
        "CompensationVoltage",
        "totalWindowWidth",
        "exclude_num",
    ]
    inputFileName = "dummyFile"
    expectedOutputDf = pd.DataFrame(
        [[inputFileName] + identifierOutputData], columns=expectedColumns
    )
    outputDf = format_output_as_pandas_dataframe(inputFileName, [identifierOutputData])
    assert expectedOutputDf.equals(outputDf)


def test__score_functions__identify_all_decoys():
    isNotDecoy, isDecoy = 0, 1
    targetLibraryIdx, decoyLibraryIdx, queryIdx = 0, 1, 0

    scoreData = [
        [targetLibraryIdx, queryIdx],
        [decoyLibraryIdx, queryIdx],
    ]
    scoreDf = pd.DataFrame(scoreData, columns=["libraryIdx", "queryIdx"])
    expectedOutput = np.array([isNotDecoy, isDecoy])
    decoySet = set([decoyLibraryIdx])
    output = identify_all_decoys(decoySet, scoreDf)
    assert np.array_equal(output, expectedOutput)
