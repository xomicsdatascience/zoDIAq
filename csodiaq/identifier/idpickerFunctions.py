

def collapse_identically_connected_peptides_and_proteins(df):

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

def separate_independent_clusters(df): pass
'''
    df["cluster"] = [-1] * len(df.index)
    while -1 in df["clusters"]:
        peptidesInCluster = set()
        proteinsInCluster = set()
        noClusterAssignedDf = df[df["clusters"] == -1]
        firstElementNotInClusterIdx = clusters.index(-1)
        peptidesInCluster.add(df.iloc[firstElementNotInClusterIdx]["peptide"])
        proteinsInCluster.add(df.iloc[firstElementNotInClusterIdx]["protein"])
        tempPeptideSet = set()
        tempProteinSet = set()
        while len(tempPeptideSet) != len(peptidesInCluster) and len(tempProteinSet) != len(proteinsInCluster):
            tempProteinSet = proteinsInCluster
            tempPeptideSet = peptidesInCluster
            tempDf = df[df["peptide"].isin(peptidesInCluster) & df["protein"].isin(proteinsInCluster)]
            peptidesInCluster = set(tempDf["peptide"])
            proteinsInCluster = set(tempDf["protein"])
        tempDf = df[df["peptide"].isin(peptidesInCluster) & df["protein"].isin(proteinsInCluster)]

        df["cluster"] = [-1] * len(df.index)
'''


