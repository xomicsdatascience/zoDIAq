from csodiaq.loaders.library.libraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.library.libraryLoaderStrategyTable import (
    LibraryLoaderStrategyTable,
)
from csodiaq.loaders.library.libraryLoaderStrategyMgf import LibraryLoaderStrategyMgf

import pytest
import os
import re


def get_parent_dir():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture
def tableContext():
    tableLibPath = os.path.join(
        get_parent_dir(), "test_files", "sample_lib_table_spectrast.tsv"
    )
    return LibraryLoaderContext(tableLibPath)


def test__library_loader_context__table__initialization(tableContext):
    assert isinstance(tableContext._strategy, LibraryLoaderStrategyTable)


def test__library_loader_context__table__load_csodiaq_library_dict(tableContext):
    expectedCsodiaqLibDict = {
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
    csodiaqLibDict = tableContext.load_csodiaq_library_dict()
    assert csodiaqLibDict == expectedCsodiaqLibDict


def test__library_loader_context__table__initialization_csv():
    tableLibPath = os.path.join(
        get_parent_dir(), "test_files", "sample_lib_table_prosit.csv"
    )
    context = LibraryLoaderContext(tableLibPath)
    assert isinstance(context._strategy, LibraryLoaderStrategyTable)


@pytest.fixture
def mgfContext():
    mgfLibPath = os.path.join(get_parent_dir(), "test_files", "sample_lib_mgf.mgf")
    return LibraryLoaderContext(mgfLibPath)


def test__library_loader_context__mgf__initialization(mgfContext):
    assert isinstance(mgfContext._strategy, LibraryLoaderStrategyMgf)


def test__library_loader_context__mgf__load_csodiaq_library_dict(mgfContext):
    expectedCsodiaqLibDict = {
        (798.93, "AAAAAAAAAAAAAAAGAGAGAK"): {
            "precursorCharge": 2,
            "identifier": 'human.faims.1.1. File:"", NativeID:"scan=1" Retention Time: 1910.238',
            "proteinName": "1/sp|P55011|S12A2_HUMAN",
            "peaks": [
                (356.193, 498.5, 0),
                (427.23, 356.3, 0),
                (498.267, 380.3, 0),
                (569.304, 415.2, 0),
                (640.341, 324.3, 0),
                (798.43, 475.3, 0),
                (799.414, 2485.1, 0),
                (815.437, 274.3, 0),
                (886.474, 448.2, 0),
                (1099.585, 366.1, 0),
            ],
            "csodiaqKeyIdx": 0,
            "isDecoy": 0,
        }
    }
    csodiaqLibDict = mgfContext.load_csodiaq_library_dict()
    assert csodiaqLibDict == expectedCsodiaqLibDict


def test__library_loader_context__initialization_msp_raises_specific_error():
    parentDir = os.path.dirname(os.getcwd())
    mspLibPath = os.path.join(parentDir, "test_files", "sample_lib_msp_prosit.msp")
    errorOutput = "The .msp library format is not currently supported. If the library file was generated via prosit, please reset the output into a tab-delimited (.tsv) format."
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        context = LibraryLoaderContext(mspLibPath)
