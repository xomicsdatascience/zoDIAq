import pandas as pd

def format_output_line(libMetadata, queMetadata, matchMetadata):
    return [
        queMetadata["scan"],
        queMetadata["precursorMz"],
        libMetadata["peptide"],
        libMetadata["proteinName"],
        libMetadata["precursorMz"],
        libMetadata["precursorCharge"],
        matchMetadata["cosineSimilarityScore"],
        libMetadata["identifier"],
        queMetadata["peaksCount"],
        libMetadata["peaks"],
        matchMetadata["shared"],
        matchMetadata["ionCount"],
        queMetadata["CV"],
        queMetadata["windowWidth"],
        matchMetadata["maccScore"],
        #matchMetadata["excludeNum"],
    ]

def extract_metadata_from_match_dataframe_groupby(group):
    return {
        "shared": len(group.index),
        "ionCount": sum(group["queryIntensity"]),
        #TODO: see TODO in extract_metadata_from_match_and_score_dataframes test
        #"extractNum": 0
    }

def extract_metadata_from_score_dataframe(df):
    #TODO: there has to be a cleaner way to do this
    maccDict = df.set_index(["libraryIdx", "queryIdx"])["maccScore"].to_dict()
    cosineDict = df.set_index(["libraryIdx", "queryIdx"])["cosineScore"].to_dict()
    outputDict = {
        k : {
            "maccScore": maccDict[k],
            "cosineSimilarityScore": cosineDict[k]
        } for k in maccDict.keys()}
    return outputDict


def extract_metadata_from_match_and_score_dataframes(matchDf, scoreDf):
    matchDict = {k: extract_metadata_from_match_dataframe_groupby(v) for k, v in matchDf.groupby(["libraryIdx","queryIdx"])}
    scoreDict = extract_metadata_from_score_dataframe(scoreDf)
    metadataDict = {
        k: {**matchDict[k], **scoreDict[k]} for k in scoreDict.keys()
    }
    return metadataDict

def format_output_as_pandas_dataframe(inputFileName, outputData):
    columns = [
        'scan',
        'MzEXP',
        'peptide',
        'protein',
        'MzLIB',
        'zLIB',
        'cosine',
        'name',
        'Peak(Query)',
        'Peaks(Library)',
        'shared',
        'ionCount',
        'CompensationVoltage',
        'totalWindowWidth',
        'MaCC_Score',
    ]
    outputDf = pd.DataFrame(outputData, columns=columns)
    outputDf.insert(0, 'fileName', [inputFileName] * len(outputDf.index))
    return outputDf

