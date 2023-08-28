import os
import subprocess
from tempfile import TemporaryDirectory
import pytest
import pandas as pd
import numpy as np
from . import create_input_template_for_scoring_module

@pytest.fixture(scope="module")
def systemTestFileDirectory():
    return TemporaryDirectory(prefix="csodiaq_scoring_system_test_files_")

@pytest.fixture(scope="module")
def inputFileDirectory(systemTestFileDirectory):
    inputDirectory = os.path.join(systemTestFileDirectory.name, "input_files")
    os.mkdir(inputDirectory)
    return inputDirectory

@pytest.fixture(scope="module")
def expectedOutputDirectory(systemTestFileDirectory):
    expectedOutputDirectory = os.path.join(systemTestFileDirectory.name, "expected_output_files")
    os.mkdir(expectedOutputDirectory)
    return expectedOutputDirectory

def assert_pandas_dataframes_are_equal(expectedDf, df):
    """
    NOTE: This function does NOT compare the specific score calculated (ie MaCC score). It is implied in the rank setup.
    """
    for columnName in expectedDf.columns:
        assert columnName in df.columns
        expectedColumn = np.array(expectedDf[columnName])
        column = np.array(df[columnName])
        if expectedColumn.dtype.kind in np.typecodes["AllFloat"]:
            np.testing.assert_array_almost_equal(expectedColumn, column)
        else:
            np.testing.assert_array_equal(expectedColumn, column)

def test__scoring__baseline_macc_score_run(inputFileDirectory, expectedOutputDirectory):
    inputTemplateDf = create_input_template_for_scoring_module()
    maccDf = add_macc_rank_values_to_input_template_and_remove_rank_label(inputTemplateDf)
    inputFilePath = os.path.join(inputFileDirectory, 'CsoDIAq-file_fullOutput.csv')
    maccDf.to_csv(inputFilePath, index=False)
    expectedSpectralOutput = make_expected_spectral_output(maccDf)
    expectedSpectralOutput.to_csv(os.path.join(expectedOutputDirectory, 'spectralFDR.csv'), index=False)
    expectedPeptideOutput = make_expected_peptide_output(maccDf)
    expectedPeptideOutput.to_csv(os.path.join(expectedOutputDirectory, 'peptideFDR.csv'), index=False)
    expectedProteinOutput = make_expected_protein_output(expectedPeptideOutput)
    expectedProteinOutput.to_csv(os.path.join(expectedOutputDirectory, 'proteinFDR.csv'), index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectory,
    ]
    subprocess.run(args, capture_output=True)
    assert 'fdrScores-macc' in os.listdir(inputFileDirectory)
    outputDirPath = os.path.join(inputFileDirectory, 'fdrScores-macc')
    outputDirContents = os.listdir(outputDirPath)
    assert 'CsoDIAq-file_fullOutput_spectralFDR.csv' in outputDirContents
    spectralOutputDf = pd.read_csv(os.path.join(outputDirPath, 'CsoDIAq-file_fullOutput_spectralFDR.csv'))
    assert_pandas_dataframes_are_equal(expectedSpectralOutput, spectralOutputDf)
    assert 'CsoDIAq-file_fullOutput_peptideFDR.csv' in outputDirContents
    peptideOutputDf = pd.read_csv(os.path.join(outputDirPath, 'CsoDIAq-file_fullOutput_peptideFDR.csv'))
    assert_pandas_dataframes_are_equal(expectedPeptideOutput, peptideOutputDf)
    assert 'CsoDIAq-file_fullOutput_proteinFDR.csv' in outputDirContents
    proteinOutputDf = pd.read_csv(os.path.join(outputDirPath, 'CsoDIAq-file_fullOutput_proteinFDR.csv'))
    assert_pandas_dataframes_are_equal(expectedProteinOutput, proteinOutputDf)


def add_macc_rank_values_to_input_template_and_remove_rank_label(inputTemplateDf):
    maccRanking = [
        (10, 0.92),
        (9, 0.93),
        (8, 0.94),
        (7, 0.95),
        (6, 0.96),
        (5, 0.97),
        (4, 0.98),
    ]
    maccRankValues = inputTemplateDf.groupby("rank")["rank"].transform(lambda x: [maccRanking[x.iloc[0]]]*len(x))
    inputTemplateDf["shared"], inputTemplateDf["cosine"] = zip(*maccRankValues)
    del inputTemplateDf["rank"]
    return inputTemplateDf

def make_expected_spectral_output(scoredDf):
    expectedSpectralOutput = scoredDf.copy()
    expectedSpectralOutput = expectedSpectralOutput.iloc[:414,:].reset_index(drop=True)
    spectralFdr = [0] * len(expectedSpectralOutput.index)
    lastTargetIdx = 410
    decoyFdrs = [i/(lastTargetIdx+i) for i in range(1,5)]
    spectralFdr[-4:] = decoyFdrs
    expectedSpectralOutput["spectralFDR"] = spectralFdr
    return expectedSpectralOutput

def make_expected_peptide_output(scoredDf):
    expectedPeptideOutput = scoredDf.copy()
    lastTargetIdx = 210
    expectedPeptideOutput = expectedPeptideOutput.iloc[:lastTargetIdx,:]
    decoysInPeptideOutput = scoredDf.iloc[410:412,:]
    expectedPeptideOutput = pd.concat([expectedPeptideOutput, decoysInPeptideOutput]).reset_index(drop=True)
    peptideFdr = [0] * len(expectedPeptideOutput.index)
    decoyFdrs = [i/(lastTargetIdx+i) for i in range(1,3)]
    peptideFdr[-2:] = decoyFdrs
    expectedPeptideOutput["peptideFDR"] = peptideFdr
    return expectedPeptideOutput

def make_expected_protein_output(peptideDf):
    proteinDf = peptideDf.copy()
    proteinDf = proteinDf.reindex(list(proteinDf.index)+[1]).sort_index()
    leadingProteins = list(proteinDf["protein"])
    idPickerLeadingProteins = [
        "1/protein7",
        "2/protein4/protein9",
        "1/protein6",
        "1/protein1",
        "1/protein1",
        "1/protein7",
        "1/protein6",
        "1/protein1",
        "1/protein1",
        "1/protein1",
        "2/protein4/protein9",
    ]
    leadingProteins[:len(idPickerLeadingProteins)] = idPickerLeadingProteins
    proteinDf["leadingProtein"] = leadingProteins
    proteinDf["proteinCosine"] = proteinDf["cosine"]
    numUniqueLeadingProteins = len(set(idPickerLeadingProteins))
    lastTargetIdx = len(set(proteinDf["leadingProtein"])) - 2
    proteinFdrs = list(proteinDf["peptideFDR"])
    decoyFdrs = [i/(lastTargetIdx+i) for i in range(1,3)]
    proteinFdrs[-2:] = decoyFdrs
    proteinDf["leadingProteinFDR"] = proteinFdrs
    proteinDf["uniquePeptide"] = [0] * len(idPickerLeadingProteins) + [1] * (len(proteinDf.index) - len(idPickerLeadingProteins))
    return proteinDf




