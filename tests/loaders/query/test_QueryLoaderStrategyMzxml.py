from csodiaq.loaders.query.QueryLoaderStrategy import QueryLoaderStrategy
from csodiaq.loaders.query.QueryLoaderStrategyMzxml import QueryLoaderStrategyMzxml
import pytest
import os
from expectedPooledQueryPeaks import (
    expectedOutputPoolSingleScan,
    expectedOutputPoolMultipleScans,
)


def get_parent_dir():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture
def testFile():
    return os.path.join(get_parent_dir(), "test_files", "sample_query_mzxml.mzXML")


@pytest.fixture
def loader(testFile):
    return QueryLoaderStrategyMzxml(testFile)


def test__library_loader_strategy_table__initialization(loader):
    assert hasattr(loader, "filePath")


def test__query_loader_strategy_mzxml__map_query_scan_ids_to_dia_mz_windows(loader):
    expectedOutputDict = {
        (781.400024414063, 2.0): ["1", "119"],
        (816.419982910156, 2.0): ["2"],
    }
    outputDict = loader.map_query_scan_ids_to_dia_mz_windows()
    assert isinstance(outputDict, dict)
    for key, value in expectedOutputDict.items():
        assert key in outputDict
        assert outputDict[key] == value


def test__query_loader_strategy_mzxml__extract_metadata_from_query_scans(loader):
    expectedOutputDict = {
        "1": {
            "precursorMz": 781.400024414063,
            "windowWidth": 2.0,
            "peaksCount": 266,
            "CV": -30.0,
        },
        "2": {
            "precursorMz": 816.419982910156,
            "windowWidth": 2.0,
            "peaksCount": 217,
            "CV": -30.0,
        },
        "119": {
            "precursorMz": 781.400024414063,
            "windowWidth": 2.0,
            "peaksCount": 660,
            "CV": -40.0,
        },
    }
    outputDict = loader.extract_metadata_from_query_scans()
    for key, value in expectedOutputDict.items():
        assert key in outputDict
        for metadataKey, metadataValue in expectedOutputDict[key].items():
            assert metadataKey in outputDict[key]
            assert outputDict[key][metadataKey] == metadataValue


def test__query_loader_strategy_mzxml__pool_query_scans__single_scan(loader):
    singleScan = ["2"]
    outputPool = loader.pool_peaks_of_query_scans(singleScan)
    for i in range(len(expectedOutputPoolSingleScan)):
        assert outputPool[i] == expectedOutputPoolSingleScan[i]


def test__query_loader_strategy_mzxml__pool_query_scans__multiple_scan(loader):
    multipleScans = ["1", "119"]
    outputPool = loader.pool_peaks_of_query_scans(multipleScans)
    for i in range(len(expectedOutputPoolMultipleScans)):
        assert outputPool[i] == expectedOutputPoolMultipleScans[i]
