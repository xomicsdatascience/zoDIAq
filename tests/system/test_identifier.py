import pandas as pd
import numpy as np
from csodiaq.identifier import Identifier
from csodiaq.loaders.query import QueryLoaderContext
import os
import pickle
import pytest

def get_parent_dir():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_file_from_system_test_folder(file):
    return os.path.join(
        get_parent_dir(), "test_files/system", file
    )

@pytest.fixture
def commandLineArgs():
    args = {
        'command': 'id',
        'files': [
            get_file_from_system_test_folder('20190412_DI2A_1to1_tMS2_n3.mzXML')
            ],
        'library': get_file_from_system_test_folder('system_test_library.csv'),
        'outDirectory': '',
        'fragmentMassTolerance': 30,
        'correction': 0,
        'histogram': False,
        'proteinTargets': 1,
        'query': 0,
        'heavyMz': True,
        'peaks': 0,
        'commonpeptide': False,
        'commonprotein': False
    }
    return args

def assert_numeric_pandas_dataframes_are_equal(expectedDf, df, type):
    assert len(expectedDf.index) == len(df.index)
    columns = get_columns_that_should_match(type)
    for columnName in columns:
        assert columnName in df.columns
        expectedColumn = np.array(expectedDf[columnName]).sort()
        column = np.array(df[columnName]).sort()
        np.testing.assert_array_equal(expectedColumn, column)

def get_columns_that_should_match(type):
    if type == "match":
        return [
            "libraryIntensity",
            "queryIntensity",
            "ppmDifference",
        ]
    else:
        return [
            "cosineScore",
            "maccScore",
        ]

def test__identifier__main_workflow(commandLineArgs):
    identifier = Identifier(commandLineArgs)
    queryFile = commandLineArgs['files'][0]
    identifier._queryContext = QueryLoaderContext(queryFile)

    expectedMatchDf = pd.read_csv(get_file_from_system_test_folder('matchDf_precorrected.csv.gz'), compression='gzip')
    matchDf = identifier._match_library_to_query_spectra()
    assert_numeric_pandas_dataframes_are_equal(expectedMatchDf, matchDf, "match")

    expectedScoreDf = pd.read_csv(get_file_from_system_test_folder('scoreDf_precorrected.csv.gz'), compression='gzip')
    scoreDf = identifier._score_spectra_matches(matchDf)
    assert_numeric_pandas_dataframes_are_equal(expectedScoreDf, scoreDf, "score")

    expectedCorrectedMatchDf = pd.read_csv(get_file_from_system_test_folder('matchDf_postcorrected.csv.gz'), compression='gzip')
    expectedCorrectedScoreDf = pd.read_csv(get_file_from_system_test_folder('scoreDf_postcorrected.csv.gz'), compression='gzip')
    matchDf, scoreDf = identifier._apply_correction_to_dataframes(matchDf, scoreDf)
    assert_numeric_pandas_dataframes_are_equal(expectedCorrectedMatchDf, matchDf, "match")
    assert_numeric_pandas_dataframes_are_equal(expectedCorrectedScoreDf, scoreDf, "score")






