import numpy as np

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
    clusters = extract_clusters_from_dataframe_recursion(df, clusters=[])
    for clusterNum in range(len(clusters)):
        clusterIdx = clusters[clusterNum]
        clusterColumn[clusterIdx] = clusterNum
    return clusterColumn

def extract_clusters_from_dataframe_recursion(df, clusters):
    if len(df.index) == 0:
        return clusters
    initialPeptideSet, initialProteinSet = initialize_peptide_protein_sets_of_next_cluster(df)
    peptideSet, _ = identify_next_cluster_in_dataframe_recursion(df, initialPeptideSet, initialProteinSet)
    clusterDf = df[df["peptide"].isin(peptideSet)]
    clusters.append(clusterDf.index)
    df = df[~df.index.isin(clusterDf.index)]
    return extract_clusters_from_dataframe_recursion(df, clusters)

def initialize_peptide_protein_sets_of_next_cluster(df):
    peptideSet = set([df.iloc[0]["peptide"]])
    proteinSet = set([df.iloc[0]["protein"]])
    return peptideSet, proteinSet

def identify_next_cluster_in_dataframe_recursion(df, oldPeptideSet, oldProteinSet):
    newPeptideSet, newProteinSet = identify_all_matching_peptide_proteins_in_cluster_from_old_set(df, oldPeptideSet, oldProteinSet)
    if newPeptideSet == oldPeptideSet and newProteinSet == oldProteinSet:
        return oldPeptideSet, newPeptideSet
    return identify_next_cluster_in_dataframe_recursion(df, newPeptideSet, newProteinSet)

def identify_all_matching_peptide_proteins_in_cluster_from_old_set(df, peptideSet, proteinSet):
    df = df[df["peptide"].isin(peptideSet) | df["protein"].isin(proteinSet)]
    return set(df["peptide"]), set(df["protein"])