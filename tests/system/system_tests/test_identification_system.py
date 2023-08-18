import os
import pandas as pd
import subprocess
from tempfile import TemporaryDirectory, NamedTemporaryFile
import numpy as np
from csodiaq.utils import Printer

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

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

def test__system__basic_run():
    expectedOutputDf = pd.read_csv(os.path.join(get_file_from_system_test_folder('outputs'), 'query.csv'))
    inputQueryFile = os.path.join(get_file_from_system_test_folder('inputs'), 'query.mzXML')
    libraryFile = os.path.join(get_file_from_system_test_folder('libraries'), 'spectrast_test_library.csv')
    outputDir = TemporaryDirectory(prefix='csodiaq_system_test')
    args = ['csodiaq', 'id']
    args += ['-i', inputQueryFile]
    args += ['-l', libraryFile]
    args += ['-o', outputDir.name]
    args += ['-nc']
    subprocess.run(args)
    outputDirContents = os.listdir(outputDir.name)
    assert len(outputDirContents) == 1
    csodiaqDir = os.path.join(outputDir.name, outputDirContents[0])
    csodiaqDirContents = os.listdir(csodiaqDir)
    assert len(csodiaqDirContents) == 1
    outputFile = os.path.join(csodiaqDir, csodiaqDirContents[0])
    outputDf = pd.read_csv(outputFile)
    outputDf = outputDf\
        .drop(
        ["fileName", "MaCC_Score"], axis=1
    )\
        .sort_values(['cosine','MzLIB'], ascending=[False, True])\
        .reset_index(drop=True)
    assert_pandas_dataframes_are_equal(expectedOutputDf, outputDf)

