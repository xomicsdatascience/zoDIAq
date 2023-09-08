import os
import subprocess
from tempfile import TemporaryDirectory
import pytest
import pandas as pd
import numpy as np
from . import MaccScoresBreakdown, ProteinCosineEvalScoresBreakdown, StandardSample1And2Breakdown, StandardSample3Breakdown


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
    NOTE: This function does NOT compare the specific score calculated (ie MaCC score).
        It is implied in the rank setup of each breakdown.
    """
    assert list(expectedDf.index) == list(df.index)
    for columnName in expectedDf.columns:
        assert columnName in df.columns
        expectedColumn = np.array(expectedDf[columnName])
        column = np.array(df[columnName])
        if expectedColumn.dtype.kind in np.typecodes["AllFloat"]:
            np.testing.assert_array_almost_equal(expectedColumn, column)
        else:
            np.testing.assert_array_equal(expectedColumn, column)

def assert_all_fdr_outputs_are_correct(inputFileDirectory, inputHeader, breakdown):
    assert 'fdrScores-macc' in os.listdir(inputFileDirectory)
    outputDirPath = os.path.join(inputFileDirectory, 'fdrScores-macc')
    outputDirContents = os.listdir(outputDirPath)
    assert f'CsoDIAq-file_{inputHeader}_fullOutput_spectralFDR.csv' in outputDirContents
    spectralOutputDf = pd.read_csv(os.path.join(outputDirPath, f'CsoDIAq-file_{inputHeader}_fullOutput_spectralFDR.csv'))
    assert_pandas_dataframes_are_equal(breakdown.outputDict['spectralFDR'], spectralOutputDf)
    assert f'CsoDIAq-file_{inputHeader}_fullOutput_peptideFDR.csv' in outputDirContents
    peptideOutputDf = pd.read_csv(os.path.join(outputDirPath, f'CsoDIAq-file_{inputHeader}_fullOutput_peptideFDR.csv'))
    assert_pandas_dataframes_are_equal(breakdown.outputDict['peptideFDR'], peptideOutputDf)
    assert f'CsoDIAq-file_{inputHeader}_fullOutput_proteinFDR.csv' in outputDirContents
    proteinOutputDf = pd.read_csv(os.path.join(outputDirPath, f'CsoDIAq-file_{inputHeader}_fullOutput_proteinFDR.csv'))
    assert_pandas_dataframes_are_equal(breakdown.outputDict['proteinFDR'], proteinOutputDf)

def assert_common_peptide_outputs_are_correct(inputFileDirectory, headerToBreakdownDict, hasProteins=True):
    outputDirPath = os.path.join(inputFileDirectory, 'fdrScores-macc')
    outputDirContents = os.listdir(outputDirPath)
    assert 'commonPeptides.csv' in outputDirContents
    expectedCommonPeptidesDfs = []
    for header, breakdown in headerToBreakdownDict.items():
        commonPeptidesDict = breakdown.outputDict['peptideFDR'].set_index('peptide')['ionCount'].to_dict()
        singleExpectedCommonPeptidesDf = pd.DataFrame(commonPeptidesDict, index=[f'CsoDIAq-file_{header}_fullOutput'])
        expectedCommonPeptidesDfs.append(singleExpectedCommonPeptidesDf)

    expectedCommonPeptidesDf = pd.concat(expectedCommonPeptidesDfs).sort_index().fillna(0)
    commonPeptidesDf = pd.read_csv(os.path.join(outputDirPath, 'commonPeptides.csv'), index_col=0).sort_index()
    assert_pandas_dataframes_are_equal(expectedCommonPeptidesDf, commonPeptidesDf)
    if hasProteins:
        assert 'commonProteins.csv' in outputDirContents

def test__scoring__baseline_macc_score_run(inputFileDirectory, expectedOutputDirectory):
    inputHeader = 'macc_baseline'
    maccBreakdown = MaccScoresBreakdown(expectedOutputDirectory)
    inputFileDirectoryChild = os.path.join(inputFileDirectory, inputHeader)
    os.mkdir(inputFileDirectoryChild)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{inputHeader}_fullOutput.csv')
    maccBreakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectoryChild,
    ]
    subprocess.run(args, capture_output=True)
    assert_all_fdr_outputs_are_correct(inputFileDirectoryChild, inputHeader, maccBreakdown)
    assert_common_peptide_outputs_are_correct(inputFileDirectoryChild, {inputHeader:maccBreakdown})

def test__scoring__evaluate_protein_cosine_column_appears_correctly(inputFileDirectory, expectedOutputDirectory):
    inputHeader = 'protein_cosine_eval'
    proteinCosineBreakdown = ProteinCosineEvalScoresBreakdown(expectedOutputDirectory)
    inputFileDirectoryChild = os.path.join(inputFileDirectory, inputHeader)
    os.mkdir(inputFileDirectoryChild)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{inputHeader}_fullOutput.csv')
    proteinCosineBreakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectoryChild,
    ]
    subprocess.run(args, capture_output=True)
    assert_all_fdr_outputs_are_correct(inputFileDirectoryChild, inputHeader, proteinCosineBreakdown)
    assert_common_peptide_outputs_are_correct(inputFileDirectoryChild, {inputHeader: proteinCosineBreakdown})

def test__scoring__evaluate_common_peptides_dataframe_forms_correctly_with_multiple_files(inputFileDirectory, expectedOutputDirectory):
    inputHeader = 'common_peptides_multiple'
    inputFileDirectoryChild = os.path.join(inputFileDirectory, inputHeader)
    os.mkdir(inputFileDirectoryChild)

    maccHeader = 'macc_baseline'
    maccBreakdown = MaccScoresBreakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{maccHeader}_fullOutput.csv')
    maccBreakdown.inputDf.to_csv(inputFilePath, index=False)

    proteinCosineHeader = 'protein_cosine_eval'
    proteinCosineBreakdown = ProteinCosineEvalScoresBreakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{proteinCosineHeader}_fullOutput.csv')
    proteinCosineBreakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectoryChild,
    ]
    subprocess.run(args, capture_output=True)
    assert_common_peptide_outputs_are_correct(inputFileDirectoryChild, {
        maccHeader: maccBreakdown,
        proteinCosineHeader: proteinCosineBreakdown,
    })

def evaluate_standard_protein_quant_comparison(inputFileDirectory, expectedOutputDirectory, method):
    inputHeader = f'standard_common_proteins_{method}_method'

    inputFileDirectoryChild = os.path.join(inputFileDirectory, inputHeader)
    os.mkdir(inputFileDirectoryChild)

    sample1Header = 'standard_common_protein_eval_sample_1'
    sample1Breakdown = StandardSample1And2Breakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{sample1Header}_fullOutput.csv')
    sample1Breakdown.inputDf.to_csv(inputFilePath, index=False)
    sample2Header = 'standard_common_protein_eval_sample_2'
    sample2Breakdown = StandardSample1And2Breakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{sample2Header}_fullOutput.csv')
    sample2Breakdown.inputDf.to_csv(inputFilePath, index=False)
    sample3Header = 'standard_common_protein_eval_sample_3'
    sample3Breakdown = StandardSample3Breakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectoryChild, f'CsoDIAq-file_{sample3Header}_fullOutput.csv')
    sample3Breakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectoryChild,
        "-p",
        method,
    ]
    subprocess.run(args, capture_output=True)
    assert_common_peptide_outputs_are_correct(inputFileDirectoryChild, {
        sample1Header: sample1Breakdown,
        sample2Header: sample2Breakdown,
        sample3Header: sample3Breakdown,
    })

    outputDirPath = os.path.join(inputFileDirectoryChild, 'fdrScores-macc')
    outputDirContents = os.listdir(outputDirPath)
    assert 'commonProteins.csv' in outputDirContents
    expectedCommonProteinDf = pd.DataFrame([
        [100.0, 100.0],
        [100.0, 100.0],
        [0.0, 100.0],
    ], columns=['protein','proteinX'], index=[f'CsoDIAq-file_{sample1Header}_fullOutput',f'CsoDIAq-file_{sample2Header}_fullOutput',f'CsoDIAq-file_{sample3Header}_fullOutput'])
    commonProteinDf = pd.read_csv(os.path.join(outputDirPath, 'commonProteins.csv'), index_col=0).sort_index()
    assert_pandas_dataframes_are_equal(expectedCommonProteinDf, commonProteinDf)

def test__scoring__evaluate_common_proteins_has_zero_quantity_in_missing_sample__averaging_method(inputFileDirectory, expectedOutputDirectory):
    evaluate_standard_protein_quant_comparison(inputFileDirectory, expectedOutputDirectory, method='average')

def test__scoring__evaluate_common_proteins_has_zero_quantity_in_missing_sample__maxlfq_method(inputFileDirectory, expectedOutputDirectory):
    evaluate_standard_protein_quant_comparison(inputFileDirectory, expectedOutputDirectory, method='maxlfq')

