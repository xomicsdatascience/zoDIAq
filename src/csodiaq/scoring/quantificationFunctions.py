import pandas as pd
from csodiaq.utils import format_protein_string_to_list
from scipy import mean
from collections import defaultdict
from itertools import chain

def calculate_ion_count_from_peptides_of_protein(ionCountList):
    return mean(ionCountList)


def calculate_ion_count_for_each_protein_in_protein_fdr_df(proteinDf):
    separateProteinData = []
    for _, row in proteinDf.iterrows():
        proteins = format_protein_string_to_list(row["leadingProtein"])
        for protein in proteins:
            separateProteinData.append([protein, row["ionCount"]])
    separateProteinDf = pd.DataFrame(
        separateProteinData, columns=["protein", "ionCount"]
    )
    proteinIonCountDf = (
        separateProteinDf.groupby("protein")
        .apply(lambda x: calculate_ion_count_from_peptides_of_protein(x["ionCount"]))
        .reset_index(name="ionCount")
    )
    return proteinIonCountDf


def compile_ion_count_comparison_across_runs_df(inputDfs, columnName):
    allValuesToCompare = sorted(
        set(list(chain.from_iterable([df[columnName] for df in inputDfs.values()])))
    )
    comparisonData = []
    inputFileNames = []
    for name, df in inputDfs.items():
        comparisonData.append(
            extract_all_ion_counts_from_df(df, allValuesToCompare, columnName)
        )
        inputFileNames.append(name)
    outputDf = pd.DataFrame(
        comparisonData, columns=allValuesToCompare, index=inputFileNames
    )
    return outputDf


def extract_all_ion_counts_from_df(df, allValuesToCompare, columnName):
    valueDict = defaultdict(
        int, pd.Series(df["ionCount"].values, index=df[columnName]).to_dict()
    )
    return [valueDict[x] for x in allValuesToCompare]