from csodiaq.identifier.identifierFunctions import format_output_line, extract_metadata_from_match_and_score_dataframes
import pandas as pd

def test__identifier_functions__format_output_line():
    libDict = {
        "peptide":"testPeptide",
        "protein":"testProtein",
        "precursorMz":0,
        "precursorCharge":1,
        "identifier":"testIdentifiers",
        "peaks":2,
    }
    queryDict = {
        "scan":"3",
        "precursorMz":4,
        "peaks":5,
        "compensationVoltage":6,
        "windowWidth":7,
    }
    matchDict = {
        "cosine":8,
        "shared":9,
        "ionCount":10,
        "maccScore":11,
        #"excludeNum":12,
    }
    expectedOutput = [
        "3",
        4,
        "testPeptide",
        "testProtein",
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
        #12,
    ]
    output = format_output_line(libDict, queryDict, matchDict)
    assert output == expectedOutput

def test__identifier_functions__extract_metadata_from_match_and_score_dataframes():
    lib1Idx = 0
    lib2Idx = 1
    queryIdx = 0
    genericIntensity = 100.0
    genericPpmDifference = 10.0
    lib1PeakCount = 3
    lib2PeakCount = 4
    lib1CosineScore = 1.0
    lib2CosineScore = 0.9
    lib1Score = 0.8
    lib2Score = 0.7
    lib1Match = [[lib1Idx, genericIntensity, queryIdx, genericIntensity,
                         genericPpmDifference]] * lib1PeakCount
    lib2Match = [[lib2Idx, genericIntensity, queryIdx, genericIntensity, genericPpmDifference]] * lib2PeakCount
    columns = ["libraryIdx", "libraryIntensity", "queryIdx", "queryIntensity", "ppmDifference"]
    matchDf = pd.DataFrame(lib1Match + lib2Match, columns=columns)

    scoreData = [
        [lib1Idx, queryIdx, lib1CosineScore, lib1Score],
        [lib2Idx, queryIdx, lib2CosineScore, lib2Score],
    ]
    scoreDf = pd.DataFrame(scoreData, columns=["libraryIdx", "queryIdx", "cosineScore", "maccScore"])
    expectedOutput = {
        (lib1Idx, queryIdx): {
            "cosineSimilarityScore": lib1CosineScore,
            "peaks": lib1PeakCount,
            "ionCount": lib1PeakCount * genericIntensity,
            "maccScore": lib1Score,
            # TODO: calculate excludeNum
            #"excludeNum": 0
        },
        (lib2Idx, queryIdx): {
            "cosineSimilarityScore": lib2CosineScore,
            "peaks": lib2PeakCount,
            "ionCount": lib2PeakCount * genericIntensity,
            "maccScore": lib2Score,
            #"excludeNum": 0
        },
    }
    output = extract_metadata_from_match_and_score_dataframes(matchDf, scoreDf)
    for idKey, metadataDict in expectedOutput.items():
        assert idKey in output
        for metadataType, metadataValue in metadataDict.items():
            assert metadataType in output[idKey]
            assert output[idKey][metadataType] == metadataValue
