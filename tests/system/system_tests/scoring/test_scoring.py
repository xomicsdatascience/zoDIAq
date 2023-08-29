import os
import subprocess
from tempfile import TemporaryDirectory
import pytest
import pandas as pd
import numpy as np
from . import MaccScoresBreakdown, ProteinCosineEvalScoresBreakdown

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
    for columnName in expectedDf.columns:
        assert columnName in df.columns
        expectedColumn = np.array(expectedDf[columnName])
        column = np.array(df[columnName])
        if expectedColumn.dtype.kind in np.typecodes["AllFloat"]:
            np.testing.assert_array_almost_equal(expectedColumn, column)
        else:
            np.testing.assert_array_equal(expectedColumn, column)

def assert_all_outputs_are_correct(inputFileDirectory, inputHeader, breakdown):
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

def test__scoring__baseline_macc_score_run(inputFileDirectory, expectedOutputDirectory):
    inputHeader = 'macc_baseline'
    maccBreakdown = MaccScoresBreakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectory, f'CsoDIAq-file_{inputHeader}_fullOutput.csv')
    maccBreakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectory,
    ]
    subprocess.run(args, capture_output=True)
    assert_all_outputs_are_correct(inputFileDirectory, inputHeader, maccBreakdown)

def test__scoring__evaluate_protein_cosine_column_appears_correctly(inputFileDirectory, expectedOutputDirectory):
    inputHeader = 'protein_cosine_eval'
    proteinCosineBreakdown = ProteinCosineEvalScoresBreakdown(expectedOutputDirectory)
    inputFilePath = os.path.join(inputFileDirectory, f'CsoDIAq-file_{inputHeader}_fullOutput.csv')
    proteinCosineBreakdown.inputDf.to_csv(inputFilePath, index=False)

    args = [
        "csodiaq",
        "score",
        "-i",
        inputFileDirectory,
    ]
    subprocess.run(args, capture_output=True)
    assert_all_outputs_are_correct(inputFileDirectory, inputHeader, proteinCosineBreakdown)




