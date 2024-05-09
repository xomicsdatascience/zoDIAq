from zodiaq.loaders.query.queryLoaderStrategy import (
    QueryLoaderStrategy,
    precursor_mz_missing_warning_text,
)
from zodiaq.loaders.library.libraryLoaderStrategy import (
    create_peaks_from_mz_intensity_lists_and_zodiaq_key_id,
)
from collections import defaultdict
from pyteomics import mzxml
import warnings
import numpy as np

class QueryLoaderStrategyMzxml(QueryLoaderStrategy):
    """
    Concrete strategy implementation of the QueryLoaderStrategy strategy class specific for
        loading mzXML query files.
    """

    def count_number_of_ms1_and_ms2_scans_per_cycle(self):
        mzWindowToScanIdDict = defaultdict(list)
        ms1Count = 0
        ms2Count = 0
        totalCount = 0
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                scan = spec["num"]
                if "precursorMz" not in spec:
                    if ms2Count == 0:
                        ms1Count += 1
                    totalCount += 1
                    continue
                if totalCount <= ms2Count+ms1Count:
                    ms2Count += 1
                totalCount += 1
        return ms1Count, ms2Count, totalCount

    def extract_metadata_from_query_scans(self) -> dict:
        scanMetadataDict = {}
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                if "precursorMz" not in spec:
                    continue
                metadataDict = {}
                metadataDict["precursorMz"] = spec["precursorMz"][0]["precursorMz"]
                metadataDict["windowWidth"] = spec["precursorMz"][0]["windowWideness"]
                metadataDict["peaksCount"] = spec["peaksCount"]
                metadataDict["retentionTime"] = spec["retentionTime"]
                if "nameValue" in spec:
                    for key, value in spec["nameValue"].items():
                        spec[key] = value
                if "compensationVoltage" in spec:
                    CV = spec["compensationVoltage"]
                else:
                    CV = ""
                metadataDict["CV"] = CV
                scanMetadataDict[spec["num"]] = metadataDict
        return scanMetadataDict

    def get_query_file_reader(self):
        return mzxml.read(self.filePath, use_index=True)

    def pool_peaks_of_query_scans(self, scans: list, reader) -> list:
        pooledQueryPeaks = []
        for scan in scans:
            spectrum = reader.get_by_id(scan)
            pooledQueryPeaks += create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
                spectrum["m/z array"], spectrum["intensity array"], int(scan)
            )
        return sorted(pooledQueryPeaks)

    def identify_mz_window_highest_and_lowest_bound(self, scans: list, reader) -> list:
        highestWindowBound = -np.inf
        lowestWindowBound = np.inf
        for scan in scans:
            spectrum = reader.get_by_id(scan)
            assert 'precursorMz' in spectrum
            precMz = spectrum["precursorMz"][0]["precursorMz"]
            windowWidth = spectrum["precursorMz"][0]["windowWideness"]
            highBound = precMz + windowWidth / 2
            lowBound = precMz - windowWidth / 2
            if highBound > highestWindowBound:
                highestWindowBound = highBound
            if lowBound < lowestWindowBound:
                lowestWindowBound = lowBound
        assert highestWindowBound > -np.inf
        assert lowestWindowBound < np.inf
        return lowestWindowBound, highestWindowBound
