from csodiaq.identifier.poolingFunctions import _pool_library_spectra_by_mz_window, generate_pooled_library_and_query_spectra_by_mz_windows
from csodiaq.loaders.library.LibraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.query.QueryLoaderContext import QueryLoaderContext
from csodiaq.loaders.query.QueryLoaderStrategyMzxml import QueryLoaderStrategyMzxml
from expectedPooledPeaks import expectedLibraryPeaks, expectedQueryPeaks
import os
import pytest
from unittest.mock import Mock
import re
import pytest


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

def test__pooler__pool_library_spectra_by_mz_window():
    aPeaks = [(1.0,1.0,0)]
    bPeaks = [(1.0,1.0,1)]
    cPeaks = [(1.0,1.0,2)]
    dPeaks = [(1.0,1.0,3)]
    ePeaks = [(1.0,1.0,4)]
    testLibDict = {
        (18.9, "A"): {'peaks':aPeaks},
        (19.0, "B"): {'peaks':bPeaks},
        (20.0, "C"): {'peaks':cPeaks},
        (21.0, "D"): {'peaks':dPeaks},
        (21.1, "E"): {'peaks':ePeaks},
    }
    mzWindow = (20.0, 2.0)
    expectedOutput = bPeaks + cPeaks + dPeaks
    output = _pool_library_spectra_by_mz_window(mzWindow, testLibDict)
    assert output == expectedOutput

def test__pooler__pool_library_spectra_by_mz_window__throws_warning_when_no_lib_peptides_found_in_window():
    aPeaks = [(1.0,1.0,0)]
    testLibDict = {
        (18.9, "A"): {'peaks':aPeaks},
    }
    mzWindow = (30.0, 2.0)
    expectedOutput = []
    errorOutput = f"No library spectra found in the {mzWindow} m/z window. Skipping"
    with pytest.warns(SyntaxWarning, match=re.escape(errorOutput)):
        output = _pool_library_spectra_by_mz_window(mzWindow, testLibDict)
        assert output == expectedOutput


def test__pooler__generate_pooled_library_and_query_spectra_by_mz_windows(libraryFile, queryFile):
    errorOutput = f"No library spectra found in the {(781.400024414063, 2.0)} m/z window. Skipping"
    with pytest.warns(SyntaxWarning, match=re.escape(errorOutput)):
        for libPeaks, queryPeaks in generate_pooled_library_and_query_spectra_by_mz_windows(libraryFile, queryFile):
            assert libPeaks == expectedLibraryPeaks
            assert queryPeaks == expectedQueryPeaks

