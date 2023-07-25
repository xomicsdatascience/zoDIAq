import numpy as np
import pandas as pd

def initialize__format_peptide_protein_connections(df):
    peptideProteinConnections = []

    for i in range(len(df)):
        peptide = df["peptide"].loc[i]
        proteinGroup = df["protein"].loc[i]

        for protein in proteinGroup.split("/")[1:]:
            peptideProteinConnections.append((peptide, protein))
    return pd.DataFrame(peptideProteinConnections, columns=["peptide","protein"])

def collapse__group_identically_connected_peptides_and_proteins(df):
    df["peptideGroup"] = group_nodes_by_identical_edges(df, isPeptideNodes=True)
    df["proteinGroup"] = group_nodes_by_identical_edges(df, isPeptideNodes=False)

    df = df[["peptideGroup","proteinGroup"]]
    df = df.drop_duplicates(["peptideGroup","proteinGroup"]).reset_index(drop=True)
    df.columns = ["peptide","protein"]

    return df

def group_nodes_by_identical_edges(df, isPeptideNodes):
    if isPeptideNodes:
        mainNode = "peptide"
        connectedNode = "protein"
    else:
        mainNode = "protein"
        connectedNode = "peptide"
    tempDf = df.copy()
    tempDf["allConnections"] = df.groupby(mainNode)[connectedNode].transform(lambda x: [tuple(x.tolist())] * len(x))
    return tempDf.groupby("allConnections")[mainNode].transform(lambda x: [tuple(sorted(set(x.tolist())))] * len(x))

def separate__identify_and_label_independent_clusters(df):
    clusterColumn = np.array([-1] * len(df.index))
    clusters = extract_clusters_from_dataframe(df)
    for clusterNum in range(len(clusters)):
        clusterIdx = clusters[clusterNum]
        clusterColumn[clusterIdx] = clusterNum
    return clusterColumn

def extract_clusters_from_dataframe(df):
    clusters = []
    while len(df.index) > 0:
        initialPeptideSet, initialProteinSet = initialize_peptide_protein_sets_of_next_cluster(df)
        peptideSet, _ = identify_next_cluster_in_dataframe_recursively(df, initialPeptideSet, initialProteinSet)
        clusterDf = df[df["peptide"].isin(peptideSet)]
        clusters.append(clusterDf.index)
        df = df[~df.index.isin(clusterDf.index)]
    return clusters

def initialize_peptide_protein_sets_of_next_cluster(df):
    peptideSet = set([df.iloc[0]["peptide"]])
    proteinSet = set([df.iloc[0]["protein"]])
    return peptideSet, proteinSet

def identify_next_cluster_in_dataframe_recursively(df, oldPeptideSet, oldProteinSet):
    newPeptideSet, newProteinSet = identify_all_matching_peptide_proteins_in_cluster_from_old_set(df, oldPeptideSet, oldProteinSet)
    if newPeptideSet == oldPeptideSet and newProteinSet == oldProteinSet:
        return oldPeptideSet, newPeptideSet
    return identify_next_cluster_in_dataframe_recursively(df, newPeptideSet, newProteinSet)

def identify_all_matching_peptide_proteins_in_cluster_from_old_set(df, peptideSet, proteinSet):
    df = df[df["peptide"].isin(peptideSet) | df["protein"].isin(proteinSet)]
    return set(df["peptide"]), set(df["protein"])

def reduce__identify_minimum_number_of_most_connected_proteins(df):
    allAcceptedProteins = set()
    for _, clusterDf in df.groupby("cluster"):
        sortedClusterDf, initialAcceptedProteins = initialize_protein_identification_recursion_parameters(clusterDf)
        allAcceptedProteins.update(identify_acceptable_proteins_recursively(sortedClusterDf, initialAcceptedProteins))
    return allAcceptedProteins

def initialize_protein_identification_recursion_parameters(clusterDf):
    sortedClusterDf = sort_dataframe_by_descending_protein_count(clusterDf)
    initialAcceptedProteins = set([sortedClusterDf.iloc[0]["protein"]])
    return sortedClusterDf, initialAcceptedProteins

def sort_dataframe_by_descending_protein_count(df):
    return df.assign(freq=df \
                     .groupby("protein")["protein"] \
                     .transform('count'))\
             .sort_values(by=["freq","protein"],ascending=[False,True]) \
             .drop(["freq"], axis=1)

def identify_acceptable_proteins_recursively(sortedClusterDf, acceptedProteinSet):
    acceptedProteinDf = sortedClusterDf[sortedClusterDf["protein"].isin(acceptedProteinSet)]
    unclaimedPeptideDf = sortedClusterDf[~sortedClusterDf["peptide"].isin(acceptedProteinDf["peptide"])]
    if len(unclaimedPeptideDf.index) == 0:
        return set(acceptedProteinDf["protein"])
    nextProtein = unclaimedPeptideDf.iloc[0]["protein"]
    acceptedProteinSet.add(nextProtein)
    return identify_acceptable_proteins_recursively(sortedClusterDf, acceptedProteinSet)

def identify_high_confidence_proteins(df):
    peptideProteinEdgeDf = initialize__format_peptide_protein_connections(df)
    peptideProteinEdgeDf = collapse__group_identically_connected_peptides_and_proteins(peptideProteinEdgeDf)
    peptideProteinEdgeDf["cluster"] = separate__identify_and_label_independent_clusters(peptideProteinEdgeDf)
    return reduce__identify_minimum_number_of_most_connected_proteins(peptideProteinEdgeDf)

