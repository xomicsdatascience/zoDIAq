import os
import pandas as pd
import numpy as np
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
import numpy as np
import pytest
from . import (
    create_template_library_dataframe,
    BaselineSpectraBreakdown,
    spectrastColumns,
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
    args = ["csodiaq", "id"]
    args += ["-i", inputQueryFile]
    args += ["-l", libraryFile]
    args += ["-o", outputDir.name]
    args += ["-nc"]
    subprocess.run(args)
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
