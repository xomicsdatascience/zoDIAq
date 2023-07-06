from csodiaq.identifier.Pooler import Pooler
from csodiaq.loaders.query.QueryLoaderContext import QueryLoaderContext
from csodiaq.loaders.query.QueryLoaderStrategyMzxml import QueryLoaderStrategyMzxml

import os
import pytest
from unittest.mock import Mock


def get_parent_dir():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def libraryFile():
    return os.path.join(
        get_parent_dir(), "test_files", "sample_lib_table_spectrast.tsv"
    )


@pytest.fixture
def queryFile():
    return os.path.join(get_parent_dir(), "test_files", "sample_query_mzxml.mzXML")


def assert_final_dict_output_matches_expected(outputDict, expectedOutputDict):
    for csodiaqLibKey in expectedOutputDict:
        assert csodiaqLibKey in outputDict
        for libEntryKey in expectedOutputDict[csodiaqLibKey]:
            assert libEntryKey in outputDict[csodiaqLibKey]
            assert (
                outputDict[csodiaqLibKey][libEntryKey]
                == expectedOutputDict[csodiaqLibKey][libEntryKey]
            )


def test__pooler__initialization(libraryFile, queryFile):
    pooler = Pooler(libraryFile, queryFile)
    assert hasattr(pooler, "libraryDict")
    expectedOutputDict = {
        (516.801083027, "YRPGTVALR"): {
            "precursorCharge": 2,
            "identifier": "51327_YRPGTVALR_2",
            "proteinName": "5/sp|Q71DI3|H32_HUMAN/sp|Q6NXT2|H3C_HUMAN/sp|Q16695|H31T_HUMAN/sp|P84243|H33_HUMAN/sp|P68431|H31_HUMAN",
            "peaks": [
                (429.745245163, 39.7, 0),
                (435.269418758, 84.0, 0),
                (458.308543918, 61.0, 0),
                (559.3562223920001, 67.3, 0),
                (575.293622136, 57.3, 0),
                (616.3776861179999, 124.8, 0),
                (674.362036054, 143.6, 0),
                (713.4304499719999, 1795.5, 0),
                (745.399149844, 324.1, 0),
                (858.483213826, 271.9, 0),
            ],
            "csodiaqKeyIdx": 0,
            "isDecoy": 0,
        }
    }
    assert_final_dict_output_matches_expected(pooler.libraryDict, expectedOutputDict)
    assert hasattr(pooler, "queryContext")
    assert isinstance(pooler.queryContext, QueryLoaderContext)
    assert isinstance(pooler.queryContext._strategy, QueryLoaderStrategyMzxml)


@pytest.fixture
def pooler(libraryFile, queryFile):
    pooler = Pooler(libraryFile, queryFile)


def test__pooler__generate_pooled_library_and_query_spectra_by_mz_windows():
    pass
