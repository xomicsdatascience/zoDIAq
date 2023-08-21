import os
import pandas as pd
import numpy as np
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
import numpy as np
import pytest
from . import (
    create_template_library_dataframe,
    separate_target_and_decoy_library_spectra,
    make_mz_window_spectra_summary_for_library_spectra,
    add_query_file_components_to_mz_window_spectra_summary,
    spectrastColumns,
    make_output_df_from_spectra_breakdown,
    organize_query_scan_data_for_writing,
    write_query_scan_data_to_mzml_file,
    convert_mzml_file_to_mzxml_file,
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
def librarySpectraBreakdown(libraryTemplateDataFrame):
    targetDecoySpectra = separate_target_and_decoy_library_spectra(
        libraryTemplateDataFrame
    )
    return make_mz_window_spectra_summary_for_library_spectra(targetDecoySpectra)


@pytest.fixture(scope="module")
def standardMzWindowSpectraBreakdown(librarySpectraBreakdown):
    return add_query_file_components_to_mz_window_spectra_summary(
        librarySpectraBreakdown
    )


@pytest.fixture(scope="module")
def inputFileDirectory(standardMzWindowSpectraBreakdown, systemTestFileDirectory):
    inputDirectory = os.path.join(systemTestFileDirectory.name, "input_files")
    os.mkdir(inputDirectory)
    mzmlFile = os.path.join(inputDirectory, "query.mzML")
    queryScanData = organize_query_scan_data_for_writing(
        standardMzWindowSpectraBreakdown
    )
    write_query_scan_data_to_mzml_file(queryScanData, mzmlFile)
    convert_mzml_file_to_mzxml_file(mzmlFile)
    return inputDirectory


@pytest.fixture(scope="module")
def standardExpectedOutputDataFrame(standardMzWindowSpectraBreakdown):
    return make_output_df_from_spectra_breakdown(standardMzWindowSpectraBreakdown)


def test__identification__standard_run(
    libraryFileDirectory, inputFileDirectory, standardExpectedOutputDataFrame
):
    inputQueryFile = os.path.join(inputFileDirectory, "query.mzXML")
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
    assert_pandas_dataframes_are_equal(standardExpectedOutputDataFrame, outputDf)
