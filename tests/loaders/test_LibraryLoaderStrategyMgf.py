import pytest
import os
from pyteomics import mgf

from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf

@pytest.fixture
def loader():
    return LibraryLoaderStrategyMgf()

def test__library_loader_strategy_mgf__initialization(loader): pass

def get_parent_dir(): return os.path.dirname(os.getcwd())

@pytest.fixture
def libFilePath():
    return os.path.join(get_parent_dir(),  'test_files', 'sample_lib_mgf.mgf')

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

def assert_final_dict_output_matches_expected(outputDict, expectedOutputDict):
    for csodiaqLibKey in expectedOutputDict:
        assert csodiaqLibKey in outputDict
        for libEntryKey in expectedOutputDict[csodiaqLibKey]:
            assert libEntryKey in outputDict[csodiaqLibKey]
            assert outputDict[csodiaqLibKey][libEntryKey] == expectedOutputDict[csodiaqLibKey][libEntryKey]

def test__library_loader_strategy_mgf__format_raw_library_object_into_csodiaq_library_dict(loader, libFilePath):
    loader._load_raw_library_object_from_file(os.path.join(get_parent_dir(),  'test_files', 'sample_lib_mgf_multiple.mgf'))
    outputDict = loader._format_raw_library_object_into_csodiaq_library_dict()
    expectedOutputDict = {
        (798.93, 'AAAAAAAAAAAAAAAGAGAGAK'): {
            'precursorCharge': 2,
            'identifier': 'human.faims.1.1. File:"", NativeID:"scan=1" Retention Time: 1910.238',
            'proteinName': '1/sp|P55011|S12A2_HUMAN',
            'peaks': [(356.193, 498.5, 0), (427.23, 356.3, 0), (498.267, 380.3, 0), (569.304, 415.2, 0), (640.341, 324.3, 0), (798.43, 475.3, 0), (799.414, 2485.1, 0), (815.437, 274.3, 0), (886.474, 448.2, 0), (1099.585, 366.1, 0)],
            'csodiaqKeyIdx': 0,
            'isDecoy': 0
        }, (798.93, 'AGAAGAAAAAAAAAAAAAAGAK'): {
            'precursorCharge': 2,
            'identifier': 'DECOY_human.faims.1.1. File:"", NativeID:"scan=1" PROTEIN: DECOY_null', 'proteinName': '1/DECOY_0_sp|P55011|S12A2_HUMAN',
            'peaks': [(328.162, 498.5, 1), (342.177, 498.5, 1), (470.236, 380.3, 1), (484.251, 380.3, 1), (541.273, 415.2, 1), (826.442, 475.3, 1), (826.461, 475.3, 1), (827.445, 2485.1, 1), (914.505, 448.2, 1), (1127.617, 366.1, 1)],
            'csodiaqKeyIdx': 1,
            'isDecoy': 1
        }
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)
