import os
import pandas as pd
import numpy as np
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
import numpy as np
import pytest
import re
from . import (
    create_template_library_dataframe,
    spectrastColumns,
    BaselineSpectraBreakdown,
    NoMatchSpectraBreakdown,
    NoCompensationVoltageSpectraBreakdown,
)


def get_parent_dir():
    return os.path.dirname(os.path.abspath(__file__))


def get_file_from_system_test_folder(file):
    return os.path.join(get_parent_dir(), "test_files", file)


def assert_pandas_dataframes_are_equal(expectedDf, df):
    for columnName in expectedDf.columns:
        assert columnName in df.columns
        expectedColumn = np.array(expectedDf[columnName])
        column = np.array(df[columnName])
        if expectedColumn.dtype.kind in np.typecodes["AllFloat"]:
            np.testing.assert_array_almost_equal(expectedColumn, column)
        else:
            np.testing.assert_array_equal(expectedColumn, column)


@pytest.fixture(scope="module")
def libraryTemplateDataFrame():
    return create_template_library_dataframe(
        minMz=400,
        maxMz=1600,
        precursorMzDiff=1,
        peakMzDiff=0.1,
        peakIntensityDiff=1000.0,
    )


@pytest.fixture(scope="module")
def systemTestFileDirectory():
    return TemporaryDirectory(prefix="csodiaq_system_test_files_")


@pytest.fixture(scope="module")
def libraryFileDirectory(libraryTemplateDataFrame, systemTestFileDirectory):
    libraryDirectory = os.path.join(systemTestFileDirectory.name, "library_files")
    os.mkdir(libraryDirectory)
    libraryTemplateDataFrame.to_csv(
        os.path.join(libraryDirectory, "template_test_library.csv"), index=False
    )
    spectrastDataFrame = libraryTemplateDataFrame.copy()
    spectrastDataFrame.columns = spectrastColumns
    spectrastDataFrame.to_csv(
        os.path.join(libraryDirectory, "spectrast_test_library.csv"), index=False
    )
    return libraryDirectory


@pytest.fixture(scope="module")
def inputFileDirectory(systemTestFileDirectory):
    inputDirectory = os.path.join(systemTestFileDirectory.name, "input_files")
    os.mkdir(inputDirectory)
    return inputDirectory


def test__identification__baseline_run(
    libraryTemplateDataFrame, libraryFileDirectory, inputFileDirectory
):
    baselineSpectraBreakdown = BaselineSpectraBreakdown(libraryTemplateDataFrame)
    inputFileHeader = "baseline_run"
    baselineSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, inputFileHeader
    )
    inputQueryFile = os.path.join(inputFileDirectory, f"{inputFileHeader}.mzXML")
    libraryFile = os.path.join(libraryFileDirectory, "spectrast_test_library.csv")
    outputDir = TemporaryDirectory(prefix="csodiaq_system_test")
    args = [
        "csodiaq",
        "id",
        "-i",
        inputQueryFile,
        "-l",
        libraryFile,
        "-o",
        outputDir.name,
        "-nc",
    ]
    subprocess.run(args, capture_output=True)
    outputDirContents = os.listdir(outputDir.name)
    assert len(outputDirContents) == 1
    csodiaqDir = os.path.join(outputDir.name, outputDirContents[0])
    csodiaqDirContents = os.listdir(csodiaqDir)
    assert len(csodiaqDirContents) == 1
    outputFile = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputDf = pd.read_csv(outputFile)
    outputDf = (
        outputDf.drop(["fileName", "MaCC_Score"], axis=1)
        .sort_values(["cosine", "MzLIB"], ascending=[False, True])
        .reset_index(drop=True)
    )
    assert_pandas_dataframes_are_equal(
        baselineSpectraBreakdown.expectedOutputDf, outputDf
    )


def test__identification__baseline__two_inputs_creates_two_outputs_independent_of_one_another(
    libraryTemplateDataFrame, libraryFileDirectory, inputFileDirectory
):
    baselineSpectraBreakdown = BaselineSpectraBreakdown(libraryTemplateDataFrame)
    inputFileHeader1 = "baseline_output_1"
    baselineSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, inputFileHeader1
    )
    inputFileHeader2 = "baseline_output_2"
    baselineSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, inputFileHeader2
    )
    libraryFile = os.path.join(libraryFileDirectory, "spectrast_test_library.csv")
    inputQueryFile1 = os.path.join(inputFileDirectory, f"{inputFileHeader1}.mzXML")
    inputQueryFile2 = os.path.join(inputFileDirectory, f"{inputFileHeader2}.mzXML")
    outputDir = TemporaryDirectory(prefix="csodiaq_system_test")
    args = [
        "csodiaq",
        "id",
        "-i",
        inputQueryFile1,
        "-i",
        inputQueryFile2,
        "-l",
        libraryFile,
        "-o",
        outputDir.name,
        "-nc",
    ]

    subprocess.run(args, capture_output=True)
    outputDirContents = os.listdir(outputDir.name)
    assert len(outputDirContents) == 1
    csodiaqDir = os.path.join(outputDir.name, outputDirContents[0])
    csodiaqDirContents = os.listdir(csodiaqDir)
    assert len(csodiaqDirContents) == 2
    outputFile1 = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputFile2 = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputDf1 = pd.read_csv(outputFile1)
    outputDf2 = pd.read_csv(outputFile2)
    outputDf1 = (
        outputDf1.drop(["fileName", "MaCC_Score"], axis=1)
        .sort_values(["cosine", "MzLIB"], ascending=[False, True])
        .reset_index(drop=True)
    )
    outputDf2 = (
        outputDf2.drop(["fileName", "MaCC_Score"], axis=1)
        .sort_values(["cosine", "MzLIB"], ascending=[False, True])
        .reset_index(drop=True)
    )
    assert_pandas_dataframes_are_equal(
        baselineSpectraBreakdown.expectedOutputDf, outputDf1
    )
    assert outputDf1.equals(outputDf2)


def test__identification__baseline_no_matches_raises_error(
    libraryTemplateDataFrame, libraryFileDirectory, inputFileDirectory
):
    noMatchSpectraBreakdown = NoMatchSpectraBreakdown(libraryTemplateDataFrame)
    noMatchFileHeader = "no_match_input"
    noMatchSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, noMatchFileHeader
    )
    baselineSpectraBreakdown = BaselineSpectraBreakdown(libraryTemplateDataFrame)
    inputFileHeader = "match_input"
    baselineSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, inputFileHeader
    )
    inputQueryFile = os.path.join(inputFileDirectory, f"{inputFileHeader}.mzXML")
    noMatchQueryFile = os.path.join(inputFileDirectory, f"{noMatchFileHeader}.mzXML")
    libraryFile = os.path.join(libraryFileDirectory, "spectrast_test_library.csv")
    outputDir = TemporaryDirectory(prefix="csodiaq_system_test")
    assert_no_matches_before_correction_raises_warning_and_skips_offending_file_only(
        inputQueryFile,
        noMatchQueryFile,
        libraryFile,
        outputDir,
        baselineSpectraBreakdown,
    )
    assert_no_matches_after_correction_raises_warning(
        inputQueryFile, libraryFile, outputDir, baselineSpectraBreakdown
    )


def assert_no_matches_before_correction_raises_warning_and_skips_offending_file_only(
    inputQueryFile, noMatchQueryFile, libraryFile, outputDir, baselineSpectraBreakdown
):
    args = [
        "csodiaq",
        "id",
        "-i",
        inputQueryFile,
        "-i",
        noMatchQueryFile,
        "-l",
        libraryFile,
        "-o",
        outputDir.name,
        "-nc",
    ]
    errorOutput = f"No matches found between library and query spectra. Skipping {noMatchQueryFile} file."
    output = subprocess.run(args, capture_output=True)
    assert errorOutput in str(output.stderr)
    outputDirContents = os.listdir(outputDir.name)
    assert len(outputDirContents) == 1
    csodiaqDir = os.path.join(outputDir.name, outputDirContents[0])
    csodiaqDirContents = os.listdir(csodiaqDir)
    assert len(csodiaqDirContents) == 1
    outputFile = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputDf = pd.read_csv(outputFile)
    outputDf = (
        outputDf.drop(["fileName", "MaCC_Score"], axis=1)
        .sort_values(["cosine", "MzLIB"], ascending=[False, True])
        .reset_index(drop=True)
    )
    assert_pandas_dataframes_are_equal(
        baselineSpectraBreakdown.expectedOutputDf, outputDf
    )


def assert_no_matches_after_correction_raises_warning(
    inputQueryFile, libraryFile, outputDir, baselineSpectraBreakdown
):
    args = [
        "csodiaq",
        "id",
        "-i",
        inputQueryFile,
        "-l",
        libraryFile,
        "-o",
        outputDir.name,
    ]
    errorOutput = f"No matches found between library and query spectra. Skipping {inputQueryFile} file."
    output = subprocess.run(args, capture_output=True)
    assert errorOutput in str(output.stderr)


def test__identification__no_compensation_voltage_succeeds(
    libraryTemplateDataFrame, libraryFileDirectory, inputFileDirectory
):
    noCompensationVoltageSpectraBreakdown = NoCompensationVoltageSpectraBreakdown(
        libraryTemplateDataFrame
    )
    inputFileHeader = "no_compensation_voltage"
    noCompensationVoltageSpectraBreakdown.write_query_scan_data_input_files(
        inputFileDirectory, inputFileHeader
    )
    inputQueryFile = os.path.join(inputFileDirectory, f"{inputFileHeader}.mzXML")
    libraryFile = os.path.join(libraryFileDirectory, "spectrast_test_library.csv")
    outputDir = TemporaryDirectory(prefix="csodiaq_system_test")
    args = [
        "csodiaq",
        "id",
        "-i",
        inputQueryFile,
        "-l",
        libraryFile,
        "-o",
        outputDir.name,
        "-nc",
    ]
    subprocess.run(args, capture_output=True)
    outputDirContents = os.listdir(outputDir.name)
    assert len(outputDirContents) == 1
    csodiaqDir = os.path.join(outputDir.name, outputDirContents[0])
    csodiaqDirContents = os.listdir(csodiaqDir)
    assert len(csodiaqDirContents) == 1
    outputFile = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputDf = pd.read_csv(outputFile)
    outputDf = (
        outputDf.drop(["fileName", "MaCC_Score"], axis=1)
        .sort_values(["cosine", "MzLIB"], ascending=[False, True])
        .reset_index(drop=True)
    )
    assert_pandas_dataframes_are_equal(
        noCompensationVoltageSpectraBreakdown.expectedOutputDf, outputDf
    )
