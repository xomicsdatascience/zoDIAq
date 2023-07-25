import os.path
import pandas as pd
import numpy as np
from csodiaq.identifier.scoringFunctions import calculate_fdr_rates_of_decoy_array
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)

def format_output_line(libMetadata, queMetadata, matchMetadata):
    return [
        queMetadata["scan"],
        queMetadata["precursorMz"],
        libMetadata["peptide"],
        libMetadata["proteinName"],
        libMetadata["isDecoy"],
        libMetadata["precursorMz"],
        libMetadata["precursorCharge"],
        matchMetadata["cosineSimilarityScore"],
        libMetadata["identifier"],
        queMetadata["peaksCount"],
        len(libMetadata["peaks"]),
        matchMetadata["shared"],
        matchMetadata["ionCount"],
        queMetadata["CV"],
        queMetadata["windowWidth"],
        matchMetadata["maccScore"],
        matchMetadata["exclude_num"],
    ]

def extract_metadata_from_match_and_score_dataframes(matchDf, scoreDf, queryDict):
    matchDict = {k: extract_metadata_from_match_dataframe_groupby(v, queryDict[str(k[1])]) for k, v in matchDf.groupby(["libraryIdx","queryIdx"])}
    scoreDict = extract_metadata_from_score_dataframe(scoreDf)
    metadataDict = {
        k: {**matchDict[k], **scoreDict[k]} for k in scoreDict.keys()
    }
    return metadataDict

def extract_metadata_from_match_dataframe_groupby(group, queryMetadata):
    precursorMz = queryMetadata["precursorMz"]
    groupRowsAbovePrecursorMz = group[group["queryMz"] > precursorMz]
    return {
        "shared": len(group.index),
        "ionCount": sum(groupRowsAbovePrecursorMz["queryIntensity"]),
        "exclude_num": len(group.index) - len(groupRowsAbovePrecursorMz)
    }

def extract_metadata_from_score_dataframe(df):
    maccDict = df.set_index(["libraryIdx", "queryIdx"])["maccScore"].to_dict()
    cosineDict = df.set_index(["libraryIdx", "queryIdx"])["cosineScore"].to_dict()
    outputDict = {
        k : {
            "maccScore": maccDict[k],
            "cosineSimilarityScore": cosineDict[k]
        } for k in maccDict.keys()}
    return outputDict

def format_output_as_pandas_dataframe(inputFileName, outputData):
    columns = [
        'scan',
        'MzEXP',
        'peptide',
        'protein',
        'isDecoy',
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
        'exclude_num',
    ]
    outputDf = pd.DataFrame(outputData, columns=columns)
    outputDf.insert(0, 'fileName', [inputFileName] * len(outputDf.index))
    return outputDf

def create_outfile_header(outputDir, queryFile, correction):
    outputCsodiaqTag = 'CsoDIAq-file_'
    queryFileName = '.'.join(queryFile.split('/')[-1].split('.')[:-1])
    outputFile = outputCsodiaqTag + queryFileName
    outFileHeader = os.path.join(outputDir, outputFile)
    if correction != -1: outFileHeader += '_corrected'
    return outFileHeader

def drop_duplicate_values_from_df_in_given_column(df, column):
   return df.drop_duplicates(subset=column, keep="first").reset_index(
        drop=True
    )

def identify_leading_protein_fdrs_for_leading_proteins_below_fdr_cutoff(df, fdrCutoff = 0.01):
    df = drop_duplicate_values_from_df_in_given_column(df, column="leadingProtein")
    df["FDR"] = calculate_fdr_rates_of_decoy_array(df["isDecoy"])
    df = df[df["FDR"] < fdrCutoff]
    return dict(zip(df["leadingProtein"], df["FDR"]))

def organize_peptide_df_by_leading_proteins(peptideDf, leadingProteins):
    proteinToProteinGroup = create_dictionary_that_matches_individual_proteins_to_group_the_protein_belongs_to(leadingProteins)
    proteinDf = create_dataframe_where_peptides_match_to_one_or_more_leading_proteins(peptideDf, proteinToProteinGroup)
    return proteinDf

def create_dictionary_that_matches_individual_proteins_to_group_the_protein_belongs_to(proteinGroups):
    proteinToProteinGroup = {}
    for proteinGroup in proteinGroups:
        proteinToProteinGroup.update({
        proteinGroup[i]: proteinGroup for i in range(len(proteinGroup))
    })
    return proteinToProteinGroup

def create_dataframe_where_peptides_match_to_one_or_more_leading_proteins(peptideDf, proteinToProteinGroup):
    originalProteinGroups = peptideDf["protein"].apply(format_protein_string_to_list)
    peptideToLeadingProteinMatchIdx = []
    leadingProteinColumn = []
    for i in range(len(originalProteinGroups)):
        originalProteinGroup = originalProteinGroups[i]
        for originalProtein in originalProteinGroup:
            if originalProtein in proteinToProteinGroup.keys():
                peptideToLeadingProteinMatchIdx.append(i)
                leadingProteinColumn.append(format_protein_list_to_string(proteinToProteinGroup[originalProtein]))
    proteinDf = peptideDf.copy().iloc[peptideToLeadingProteinMatchIdx]
    proteinDf["leadingProtein"] = leadingProteinColumn
    proteinDf = proteinDf.drop_duplicates(keep="first").reset_index(
        drop=True
    )
    return proteinDf

def format_protein_string_to_list(proteinString):
    return proteinString.split("/")[1:]

def format_protein_list_to_string(proteinList):
    return f"{len(proteinList)}/{'/'.join(proteinList)}"

def determine_if_peptides_are_unique_to_leading_protein(df):
    uniqueValuesDf = df.groupby("leadingProtein").filter(lambda group: len(group) == 1)
    uniquePeptides = np.array([0]*len(df.index))
    uniquePeptides[uniqueValuesDf.index] = 1
    return list(uniquePeptides)