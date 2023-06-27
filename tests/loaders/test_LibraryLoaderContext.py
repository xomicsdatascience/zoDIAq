from csodiaq.loaders.LibraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.LibraryLoaderStrategyTable import LibraryLoaderStrategyTable
from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf

import pytest
import os


@pytest.fixture
def tableContext():
    parentDir = os.path.dirname(os.getcwd())
    tableLibPath = os.path.join(parentDir,  'test_files', 'sample_lib_table_spectrast.tsv')
    return LibraryLoaderContext(tableLibPath)

def test__library_loader_context__table__initialization(tableContext):
    assert isinstance(tableContext._strategy, LibraryLoaderStrategyTable)

def test__library_loader_context__table__load_csodiaq_library_dict(tableContext):
    expectedCsodiaqLibDict = {
        (375.87322429333335, 'FANYIDKVR'): {
            'precursorCharge': 3,
            'identifier': '1_FANYIDKVR_3',
            'proteinName': '1/sp|P08670|VIME_HUMAN',
            'peaks': [(333.15573166, 1072.6, 0), (397.23197061, 2082.7, 0), (402.2823291699999, 4930.4, 0), (445.738434336, 746.7, 0), (454.253434336, 1301.3, 0), (489.771991231, 1398.4, 0), (500.2792722, 863.7, 0), (517.3092722, 10000.0, 0), (630.393336182, 8235.6, 0), (793.4566647199999, 5098.5, 0)],
            'csodiaqKeyIdx': 0,
            'isDecoy': 0
        }
    }
    csodiaqLibDict = tableContext.load_csodiaq_library_dict()
    assert csodiaqLibDict == expectedCsodiaqLibDict

@pytest.fixture
def mgfContext():
    parentDir = os.path.dirname(os.getcwd())
    mgfLibPath = os.path.join(parentDir,  'test_files', 'sample_lib_mgf.mgf')
    return LibraryLoaderContext(mgfLibPath)

def test__library_loader_context__mgf__initialization(mgfContext):
    assert isinstance(mgfContext._strategy, LibraryLoaderStrategyMgf)

def test__library_loader_context__mgf__load_csodiaq_library_dict(mgfContext):
    expectedCsodiaqLibDict = {
        (798.93, 'AAAAAAAAAAAAAAAGAGAGAK'): {
            'precursorCharge': 2,
            'identifier': 'human.faims.1.1. File:"", NativeID:"scan=1" Retention Time: 1910.238',
            'proteinName': '1/sp|P55011|S12A2_HUMAN',
            'peaks': [(356.193, 498.5, 0), (427.23, 356.3, 0), (498.267, 380.3, 0), (569.304, 415.2, 0), (640.341, 324.3, 0), (798.43, 475.3, 0), (799.414, 2485.1, 0), (815.437, 274.3, 0), (886.474, 448.2, 0), (1099.585, 366.1, 0)],
            'csodiaqKeyIdx': 0,
            'isDecoy': 0
        }
    }
    csodiaqLibDict = mgfContext.load_csodiaq_library_dict()
    assert csodiaqLibDict == expectedCsodiaqLibDict
