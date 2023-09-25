import pandas as pd
import numpy as np
from zodiaq.utils import format_protein_string_to_list
from .idpickerFunctions import initialize__format_peptide_protein_connections
from scipy import mean
from collections import defaultdict
from itertools import chain
import scipy.linalg as linalg


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


def compile_common_protein_quantification_file(
    proteinDfs, commonPeptidesDf, proteinQuantificationMethod, minNumDifferences
):
    headerToProteinPresenceDict = {header: set(proteinDf["leadingProtein"]) for header, proteinDf in proteinDfs.items()}
    proteinPeptideDict = (
        pd.concat(list(proteinDfs.values()))
        .groupby("leadingProtein")["peptide"]
        .apply(set)
        .to_dict()
    )
    proteinQuantitiesDict = {}
    for protein, peptideSet in proteinPeptideDict.items():
        peptideQuantityDf = commonPeptidesDf[list(peptideSet)]
        peptideQuantityDf = set_non_present_protein_levels_to_zero(
            peptideQuantityDf, protein, headerToProteinPresenceDict
        )
        if proteinQuantificationMethod == "sum":
            proteinSampleQuantities = ion_count_sum(peptideQuantityDf)
        elif proteinQuantificationMethod == "maxlfq":
            proteinSampleQuantities = run_maxlfq_with_normalizations(peptideQuantityDf, minNumDifferences)
        proteinQuantitiesDict[protein] = proteinSampleQuantities
    commonProteinsDf = pd.DataFrame.from_dict(proteinQuantitiesDict)
    commonProteinsDf.index = commonPeptidesDf.index
    return commonProteinsDf

def ion_count_sum(sampleByPeptideMatrix):
    sampleByPeptideMatrix = sampleByPeptideMatrix.loc[:, (sampleByPeptideMatrix != 0).all(axis=0)]
    return np.array(sampleByPeptideMatrix.sum(axis=1))

def set_non_present_protein_levels_to_zero(
    peptideQuantityDf, protein, headerToProteinPresenceDict
):
    peptideQuantityDf = peptideQuantityDf.copy()
    indicesWithoutProtein = set(
        [
            key if protein not in value else None
            for key, value in headerToProteinPresenceDict.items()
        ]
    )
    peptideQuantityDf.loc[
        peptideQuantityDf.index.isin(indicesWithoutProtein), :
    ] = 0
    return peptideQuantityDf


def run_maxlfq_with_normalizations(peptideQuantityDf, minNumDifferences):
    proteinSampleQuantities = maxlfq(np.log(peptideQuantityDf).to_numpy(), minNumDifferences)
    proteinSampleQuantities = np.exp(proteinSampleQuantities)
    proteinSampleQuantities[proteinSampleQuantities == 1] = 0
    return proteinSampleQuantities


def prepare_matrices_for_cholesky_factorization(
    sampleByPeptideMatrix, minNumDifferences=2
):
    sampleNum = len(sampleByPeptideMatrix)
    A, B = np.zeros((sampleNum, sampleNum)), np.zeros(sampleNum)
    for i in range(sampleNum):
        for j in range(i + 1, sampleNum):
            sample1Peptides = sampleByPeptideMatrix[i]
            sample2Peptides = sampleByPeptideMatrix[j]
            matches = ~((sample1Peptides == 0) | (sample2Peptides == 0))
            numDifferences = np.count_nonzero(matches)
            if numDifferences < minNumDifferences:
                continue
            sample1Peptides, sample2Peptides = (
                sample1Peptides[matches],
                sample2Peptides[matches],
            )
            diff = np.median(sample1Peptides - sample2Peptides)
            A[i, i] += 1
            A[j, j] += 1
            A[i, j] = A[j, i] = -1
            B[i] += diff
            B[j] -= diff
    return A, B


def calculate_regularization_factor(A):
    regularization = np.array(np.diagonal(A))
    regularization[regularization < 2] = 1.0
    return regularization * 0.0001


def apply_regularization_to_matrices(sampleByPeptideMatrix, A, B):
    reg = calculate_regularization_factor(A)
    A[np.diag_indices_from(A)] = np.diagonal(A) + reg
    B += np.amax(sampleByPeptideMatrix, axis=1) * reg
    return A, B


def calculate_protein_intensities_across_samples_from_cholesky_factorization(A, B):
    a = linalg.cho_factor(A)
    return linalg.cho_solve(a, B)


def maxlfq(sampleByPeptideMatrix, minNumDifferences, tolerance=-10.0):
    sampleByPeptideMatrix[sampleByPeptideMatrix < tolerance] = 0
    A, B = prepare_matrices_for_cholesky_factorization(sampleByPeptideMatrix, minNumDifferences)
    unmatchedIdx = np.diagonal(A) == 0
    A, B = apply_regularization_to_matrices(sampleByPeptideMatrix, A, B)
    proteinScoreAcrossSamples = (
        calculate_protein_intensities_across_samples_from_cholesky_factorization(A, B)
    )
    proteinScoreAcrossSamples[unmatchedIdx] = 0
    return proteinScoreAcrossSamples
