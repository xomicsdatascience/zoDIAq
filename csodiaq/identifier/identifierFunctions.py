

def format_output_line(libMetadata, queMetadata, matchMetadata):
    return [
        queMetadata["scan"],
        queMetadata["precursorMz"],
        libMetadata["peptide"],
        libMetadata["protein"],
        libMetadata["precursorMz"],
        libMetadata["precursorCharge"],
        matchMetadata["cosine"],
        libMetadata["identifier"],
        queMetadata["peaks"],
        libMetadata["peaks"],
        matchMetadata["shared"],
        matchMetadata["ionCount"],
        queMetadata["compensationVoltage"],
        queMetadata["windowWidth"],
        matchMetadata["maccScore"],
        matchMetadata["excludeNum"],
    ]

def extract_metadata_from_match_dataframe_groupby(group):
    return {
        "peaks": len(group.index),
        "ionCount": sum(group["queryIntensity"])
    }

def extract_metadata_from_score_dataframe(df):
    return df.set_index(["libraryIdx", "queryIdx"])["score"].to_dict()


def extract_metadata_from_match_and_score_dataframes(matchDf, scoreDf):
    matchDict = {k: extract_metadata_from_match_dataframe_groupby(v) for k, v in matchDf.groupby(["libraryIdx","queryIdx"])}
    scoreDict = extract_metadata_from_score_dataframe(scoreDf)
    metadataDict = {
        k: {**matchDict[k], "score":scoreDict[k]} for k in scoreDict.keys()
    }
    return metadataDict
