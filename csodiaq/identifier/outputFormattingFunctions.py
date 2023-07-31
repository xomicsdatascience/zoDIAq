import os.path
import pandas as pd
import numpy as np
from csodiaq.identifier.scoringFunctions import calculate_fdr_rates_of_decoy_array, determine_index_of_fdr_cutoff
from csodiaq.identifier.idpickerFunctions import identify_high_confidence_proteins
from csodiaq.utils import format_protein_list_to_string, format_protein_string_to_list

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
    matchDict = {
        k: extract_metadata_from_match_dataframe_groupby(v, queryDict[str(k[1])])
        for k, v in matchDf.groupby(["libraryIdx", "queryIdx"])
    }
    scoreDict = extract_metadata_from_score_dataframe(scoreDf)
    metadataDict = {k: {**matchDict[k], **scoreDict[k]} for k in scoreDict.keys()}
    return metadataDict


def extract_metadata_from_match_dataframe_groupby(group, queryMetadata):
    precursorMz = queryMetadata["precursorMz"]
    groupRowsAbovePrecursorMz = group[group["queryMz"] > precursorMz]
    return {
        "shared": len(group.index),
        "ionCount": sum(groupRowsAbovePrecursorMz["queryIntensity"]),
        "exclude_num": len(group.index) - len(groupRowsAbovePrecursorMz),
    }


def extract_metadata_from_score_dataframe(df):
    maccDict = df.set_index(["libraryIdx", "queryIdx"])["maccScore"].to_dict()
    cosineDict = df.set_index(["libraryIdx", "queryIdx"])["cosineScore"].to_dict()
    outputDict = {
        k: {"maccScore": maccDict[k], "cosineSimilarityScore": cosineDict[k]}
        for k in maccDict.keys()
    }
    return outputDict


def format_output_as_pandas_dataframe(inputFileName, outputData):
    columns = [
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
        "MaCC_Score",
        "exclude_num",
    ]
    outputDf = pd.DataFrame(outputData, columns=columns)
    outputDf.insert(0, "fileName", [inputFileName] * len(outputDf.index))
    return outputDf

def drop_duplicate_values_from_df_in_given_column(df, column):
    return df.drop_duplicates(subset=column, keep="first").reset_index(drop=True)

def generate_spectral_fdr_output_from_full_output(fullDf, fdrCutoff = 0.01):
    fdrs = calculate_fdr_rates_of_decoy_array(fullDf["isDecoy"])
    scoreDfCutoffIdx = np.argmax(fdrs > fdrCutoff)
    fullDf["spectralFDR"] = fdrs
    return fullDf.iloc[:scoreDfCutoffIdx, :]

def generate_peptide_fdr_output_from_full_output(fullDf, fdrCutoff = 0.01):
    peptideDf = drop_duplicate_values_from_df_in_given_column(fullDf, "peptide")
    fdrs = calculate_fdr_rates_of_decoy_array(peptideDf["isDecoy"])
    scoreDfCutoffIdx = np.argmax(fdrs > fdrCutoff)
    peptideDf["peptideFDR"] = fdrs
    return peptideDf.iloc[:scoreDfCutoffIdx, :]

def generate_protein_fdr_output_from_peptide_fdr_output(peptideDf): # no unit test made yet - tested in system test
    highConfidenceProteins = identify_high_confidence_proteins(peptideDf)
    proteinDf = organize_peptide_df_by_leading_proteins(peptideDf, highConfidenceProteins)
    proteinFdrDict = identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff(proteinDf)
    proteinDf = proteinDf[proteinDf["leadingProtein"].isin(proteinFdrDict.keys())]
    proteinDf["proteinCosine"] = proteinDf.groupby("leadingProtein")["cosine"].transform("max")
    proteinDf["leadingProteinFDR"] = proteinDf["leadingProtein"].apply(lambda x: proteinFdrDict[x])
    proteinDf["uniquePeptide"] = determine_if_peptides_are_unique_to_leading_protein(proteinDf)
    return proteinDf

def identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff(
    df, fdrCutoff=0.01
):
    df = drop_duplicate_values_from_df_in_given_column(df, column="leadingProtein")
    df["FDR"] = calculate_fdr_rates_of_decoy_array(df["isDecoy"])
    df = df[df["FDR"] < fdrCutoff]
    return dict(zip(df["leadingProtein"], df["FDR"]))


def organize_peptide_df_by_leading_proteins(peptideDf, leadingProteins):
    """
    Creates a new "protein" dataframe more reliant on the leading protein column
        than the peptide column.

    Extended Summary
    ----------------
    In spectral output, each library-query spectra match is accounted for. In the peptide output,
        duplicate library peptides are removed, preferentially keeping the peptides with the highest
        confidence. In the protein dataframe, using a set of proteins we are confident are present
        ("leadingProteins"), we reorder the peptide dataframe around these identified proteins.
        A peptide could belong to more than one identified protein. When that happens, the row
        corresponding to that peptide is duplicated, once for each protein it belongs to.

    Parameters
    ----------
    peptideDf : pandas DataFrame
        The output of the library-query matching workflow, where duplicate identified peptides
            are removed, preferentially keeping peptides with higher identification confidence.

    leadingProteins : set
        A set of proteins determined to be present as identified by the idpicker algorithm.

    Returns
    -------

    """
    proteinToProteinGroup = create_dictionary_that_matches_individual_proteins_to_group_the_protein_belongs_to(
        leadingProteins
    )
    proteinDf = create_dataframe_where_peptides_match_to_one_or_more_leading_proteins(
        peptideDf, proteinToProteinGroup
    )
    return proteinDf


def create_dictionary_that_matches_individual_proteins_to_group_the_protein_belongs_to(
    proteinGroups,
):
    """
    Creates a dictionary that matches a protein to the protein group it belongs to.

    Extended Summary
    ----------------
    The ID Picker algorithm groups proteins with identical peptide matches together.
        The protein column of the identification output represents proteins that COULD be
        represented by the peptide, whereas the leadingProtein column represents proteins
        that are determined to ACTUALLY be present. Sometimes the proteins in the protein
        column map to more than one leadingProtein group. This dictionary is created to make
        mapping between these columns straightforward.



    Parameters
    ----------
    proteinGroups : set
        A set of protein groups, represented by strings.

    Returns
    -------
    proteinToProteinGroup : dict
        As an example, the proteinGroup '3/protein1/protein2/protein3' contains 3 proteins.
        The following entries would therefore be included in the return dictionary.
            {
                'protein1': '3/protein1/protein2/protein3',
                'protein2': '3/protein1/protein2/protein3',
                'protein3': '3/protein1/protein2/protein3',
            }
    """
    proteinToProteinGroup = {}
    for proteinGroup in proteinGroups:
        proteinToProteinGroup.update(
            {proteinGroup[i]: proteinGroup for i in range(len(proteinGroup))}
        )
    return proteinToProteinGroup


def create_dataframe_where_peptides_match_to_one_or_more_leading_proteins(
    peptideDf, proteinToProteinGroup
):
    originalProteinGroups = peptideDf["protein"].apply(format_protein_string_to_list)
    peptideToLeadingProteinMatchIdx = []
    leadingProteinColumn = []
    for i in range(len(originalProteinGroups)):
        originalProteinGroup = originalProteinGroups[i]
        for originalProtein in originalProteinGroup:
            if originalProtein in proteinToProteinGroup.keys():
                peptideToLeadingProteinMatchIdx.append(i)
                leadingProteinColumn.append(
                    format_protein_list_to_string(
                        proteinToProteinGroup[originalProtein]
                    )
                )
    proteinDf = peptideDf.copy().iloc[peptideToLeadingProteinMatchIdx]
    proteinDf["leadingProtein"] = leadingProteinColumn
    proteinDf = proteinDf.drop_duplicates(keep="first").reset_index(drop=True)
    return proteinDf

def determine_if_peptides_are_unique_to_leading_protein(proteinDf):
    """
    Determines if the peptide in the peptide column is unique to the protein
        in the leadingProtein column.

    Parameters
    ----------
    proteinDf : pandas DataFrame
        The library-query match output ordered around leading proteins.

    Returns
    -------
    uniquePeptides : list
        A list containing binary values. The length of the list corresponds to the
            dataframe provided as input. 0 indicates the peptide is NOT unique, 1 that
            it is unique.
    """
    uniqueValuesDf = proteinDf.groupby("leadingProtein").filter(
        lambda group: len(group) == 1
    )
    uniquePeptides = np.array([0] * len(proteinDf.index))
    uniquePeptides[uniqueValuesDf.index] = 1
    return list(uniquePeptides)
