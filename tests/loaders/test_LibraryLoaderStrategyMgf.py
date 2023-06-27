import pytest
import os
from pyteomics import mgf

from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf

@pytest.fixture
def loader():
    return LibraryLoaderStrategyMgf()

def test__library_loader_strategy_mgf__initialization(loader): pass

@pytest.fixture
def libFilePath():
    cwd = os.getcwd()
    parent = os.path.dirname(cwd)
    return os.path.join(parent,  'test_files', 'sample_lib_mgf.mgf')

def test__library_loader_strategy_mgf__load_raw_library_object_from_file(loader, libFilePath):
    loader._load_raw_library_object_from_file(libFilePath)
    assert hasattr(loader, 'rawLibMgf')
    assert isinstance(loader.rawLibMgf, mgf.IndexedMGF)

def test__library_loader_strategy_mgf__format_raw_library_object_into_csodiaq_library_dict(loader, libFilePath):
    loader._load_raw_library_object_from_file(libFilePath)
    outputDict = loader._format_raw_library_object_into_csodiaq_library_dict()
    expectedOutputDict = {
        (798.93, 'AAAAAAAAAAAAAAAGAGAGAK'): {
            'precursorCharge': 2,
            'identifier': 'human.faims.1.1. File:"", NativeID:"scan=1" Retention Time: 1910.238',
            'proteinName': '1/sp|P55011|S12A2_HUMAN',
            'peaks': [(356.193, 498.5, 0), (427.23, 356.3, 0), (498.267, 380.3, 0), (569.304, 415.2, 0), (640.341, 324.3, 0), (798.43, 475.3, 0), (799.414, 2485.1, 0), (815.437, 274.3, 0), (886.474, 448.2, 0), (1099.585, 366.1, 0)],
            'csodiaqKeyIdx': 0,
            'isDecoy': 0
        }
    }

    assert outputDict == expectedOutputDict

