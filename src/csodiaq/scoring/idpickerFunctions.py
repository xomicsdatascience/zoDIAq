import numpy as np
import pandas as pd
from csodiaq.utils import format_protein_string_to_list


def identify_high_confidence_proteins(peptideDf):
    """
    This function runs the full ID Picker algorithm on the peptides identified in the
        matching workflow (and the proteins associated with those peptides).

    Extended Summary
    ----------------
    The ID Picker takes a list of peptides as well as a corresponding list of proteins
        each peptide is known to be part of. The purpose of the ID Picker algorithm is
        to determine the minimum number of these proteins that are present given the
        presence of these peptides.
    NOTE:
        The function is intentionally 4 lines long, each representing one of the 4 steps
            of the ID Picker algorithm.

    Parameters
    ----------
    peptideDf : pandas DataFrame
        The output of the library-query matching workflow, where library peptides were
        identified.

    Returns
    -------
    leadingProteins : set
        A set of proteins determined with high confidence to be present.

    """
    peptideProteinEdgeDf = initialize__format_peptide_protein_connections(peptideDf)
    peptideProteinEdgeDf = collapse__group_identically_connected_peptides_and_proteins(
        peptideProteinEdgeDf
    )
    peptideProteinEdgeDf["cluster"] = separate__identify_and_label_independent_clusters(
        peptideProteinEdgeDf
    )
    return reduce__identify_minimum_number_of_most_connected_proteins(
        peptideProteinEdgeDf
    )


def initialize__format_peptide_protein_connections(peptideDf, proteinColumn='protein'):
    """
    Determines every individual peptide-protein edge connection.

    Extended Summary
    ----------------
    The "initialize" step of the ID Picker algorithm.

    Parameters
    ----------
    peptideDf : pandas DataFrame
        The output of the library-query matching workflow, where library peptides were
        identified.

    Returns
    -------
    peptideProteinConnectionsDf : pandas DataFrame
        A 2-column dataframe representing single peptide-protein relationships.

    """
    peptideProteinConnections = []
    for i in range(len(peptideDf)):
        peptide = peptideDf["peptide"].loc[i]
        proteinGroup = peptideDf[proteinColumn].loc[i]
        for protein in format_protein_string_to_list(proteinGroup):
            peptideProteinConnections.append((peptide, protein))
    return pd.DataFrame(peptideProteinConnections, columns=["peptide", "protein"])


def collapse__group_identically_connected_peptides_and_proteins(
    peptideProteinConnectionsDf,
):
    """
    Identifies redundant peptide-protein connections and groups peptides/proteins accordingly.

    Extended Summary
    ----------------
    The "collapse" step of the ID Picker algorithm.

    Parameters
    ----------
    peptideProteinConnectionsDf : pandas DataFrame
        A 2-column dataframe representing single peptide-protein relationships.

    Returns
    -------
    peptideProteinConnectionsDf : pandas DataFrame
        A 2-column dataframe representing single peptide-protein relationships, where
            redundant relationships have been collapsed into single groups.
    """
    peptideProteinConnectionsDf["peptideGroup"] = group_nodes_by_identical_edges(
        peptideProteinConnectionsDf, isPeptideNodes=True
    )
    peptideProteinConnectionsDf["proteinGroup"] = group_nodes_by_identical_edges(
        peptideProteinConnectionsDf, isPeptideNodes=False
    )
    peptideProteinConnectionsDf = peptideProteinConnectionsDf[
        ["peptideGroup", "proteinGroup"]
    ]
    peptideProteinConnectionsDf = peptideProteinConnectionsDf.drop_duplicates(
        ["peptideGroup", "proteinGroup"]
    ).reset_index(drop=True)
    peptideProteinConnectionsDf.columns = ["peptide", "protein"]
    return peptideProteinConnectionsDf


def group_nodes_by_identical_edges(df, isPeptideNodes):
    if isPeptideNodes:
        mainNode = "peptide"
        connectedNode = "protein"
    else:
        mainNode = "protein"
        connectedNode = "peptide"
    tempDf = df.copy()
    tempDf["allConnections"] = df.groupby(mainNode)[connectedNode].transform(
        lambda x: [tuple(x.tolist())] * len(x)
    )
    return tempDf.groupby("allConnections")[mainNode].transform(
        lambda x: [tuple(sorted(set(x.tolist())))] * len(x)
    )


def separate__identify_and_label_independent_clusters(peptideProteinConnectionsDf):
    """
    Determines independent connectivity networks in the peptide-protein graph.

    Extended Summary
    ----------------
    The "collapse" step of the ID Picker algorithm.

    Parameters
    ----------
    peptideProteinConnectionsDf : pandas DataFrame
        A 2-column dataframe representing single peptide-protein relationships, where
            redundant relationships have been collapsed into single groups.

    Returns
    -------
    clusterColumn : list
        A list of cluster IDs. This list can be added as a new column to the input
            dataframe, indicating which cluster each peptide-protein relationship
            belongs to.
    """
    clusterColumn = np.array([-1] * len(peptideProteinConnectionsDf.index))
    clusters = extract_clusters_from_dataframe(peptideProteinConnectionsDf)
    for clusterNum in range(len(clusters)):
        clusterIdx = clusters[clusterNum]
        clusterColumn[clusterIdx] = clusterNum
    return clusterColumn


def extract_clusters_from_dataframe(df):
    clusters = []
    while len(df.index) > 0:
        (
            initialPeptideSet,
            initialProteinSet,
        ) = initialize_peptide_protein_sets_of_next_cluster(df)
        peptideSet, _ = identify_next_cluster_in_dataframe_recursively(
            df, initialPeptideSet, initialProteinSet
        )
        clusterDf = df[df["peptide"].isin(peptideSet)]
        clusters.append(clusterDf.index)
        df = df[~df.index.isin(clusterDf.index)]
    return clusters


def initialize_peptide_protein_sets_of_next_cluster(df):
    peptideSet = set([df.iloc[0]["peptide"]])
    proteinSet = set([df.iloc[0]["protein"]])
    return peptideSet, proteinSet


def identify_next_cluster_in_dataframe_recursively(df, oldPeptideSet, oldProteinSet):
    (
        newPeptideSet,
        newProteinSet,
    ) = identify_all_matching_peptide_proteins_in_cluster_from_old_set(
        df, oldPeptideSet, oldProteinSet
    )
    if newPeptideSet == oldPeptideSet and newProteinSet == oldProteinSet:
        return oldPeptideSet, newPeptideSet
    return identify_next_cluster_in_dataframe_recursively(
        df, newPeptideSet, newProteinSet
    )


def identify_all_matching_peptide_proteins_in_cluster_from_old_set(
    df, peptideSet, proteinSet
):
    df = df[df["peptide"].isin(peptideSet) | df["protein"].isin(proteinSet)]
    return set(df["peptide"]), set(df["protein"])


def reduce__identify_minimum_number_of_most_connected_proteins(
    peptideProteinConnectionsDf,
):
    """
    Identifies the presence of proteins that are represented by the given peptides.

    Extended Summary
    ----------------
    The "reduce" step of the ID Picker algorithm.

    Parameters
    ----------
    peptideProteinConnectionsDf : pandas DataFrame
        A 3-column dataframe representing single peptide-protein relationships and
            the cluster those relationships belong to.

    Returns
    -------
    leadingProteins : set
        A set of proteins determined with high confidence to be present.
    """
    leadingProteins = set()
    for _, clusterDf in peptideProteinConnectionsDf.groupby("cluster"):
        clusterDf[
            "originalProteinCount"
        ] = label_proteins_by_original_protein_count_for_breaking_ties(clusterDf)
        (
            sortedClusterDf,
            initialAcceptedProteins,
        ) = initialize_protein_identification_recursion_parameters(clusterDf)
        leadingProteins.update(
            identify_acceptable_proteins_recursively(
                sortedClusterDf, initialAcceptedProteins
            )
        )
    return leadingProteins


def initialize_protein_identification_recursion_parameters(clusterDf):
    sortedClusterDf = sort_dataframe_by_descending_protein_count(clusterDf)
    initialAcceptedProteins = set([sortedClusterDf.iloc[0]["protein"]])
    return sortedClusterDf, initialAcceptedProteins


def sort_dataframe_by_descending_protein_count(df):
    return (
        df.assign(
            currentProteinCount=df.groupby("protein")["protein"].transform("count")
        )
        .sort_values(
            by=["currentProteinCount", "originalProteinCount", "protein"],
            ascending=[False, False, True],
        )
        .drop(["currentProteinCount"], axis=1)
    )


def label_proteins_by_original_protein_count_for_breaking_ties(df):
    return df.groupby("protein")["protein"].transform("count")


def identify_acceptable_proteins_recursively(sortedClusterDf, acceptedProteinSet):
    acceptedProteinDf = sortedClusterDf[
        sortedClusterDf["protein"].isin(acceptedProteinSet)
    ]
    unclaimedPeptideDf = sortedClusterDf[
        ~sortedClusterDf["peptide"].isin(acceptedProteinDf["peptide"])
    ]
    if len(unclaimedPeptideDf.index) == 0:
        return set(acceptedProteinDf["protein"])
    unclaimedPeptideDf = sort_dataframe_by_descending_protein_count(unclaimedPeptideDf)
    nextProtein = unclaimedPeptideDf.iloc[0]["protein"]
    acceptedProteinSet.add(nextProtein)
    return identify_acceptable_proteins_recursively(sortedClusterDf, acceptedProteinSet)
