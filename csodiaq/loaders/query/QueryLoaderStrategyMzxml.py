from csodiaq.loaders.query.QueryLoaderStrategy import QueryLoaderStrategy
from csodiaq.loaders.library.LibraryLoaderStrategy import create_peaks_from_mz_intensity_lists_and_csodiaq_key_id
from collections import defaultdict
from pyteomics import mzxml

class QueryLoaderStrategyMzxml(QueryLoaderStrategy):
    def map_query_scan_ids_to_dia_mz_windows(self) -> dict:
        mzWindowToScanIdDict = defaultdict(list)
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                if "precursorMz" not in spec:
                    continue
                scan = spec["num"]
                precMz = spec["precursorMz"][0]["precursorMz"]
                windowWidth = spec["precursorMz"][0]["windowWideness"]
                mzWindowToScanIdDict[precMz, windowWidth].append(scan)
        return dict(mzWindowToScanIdDict)

    def extract_metadata_from_query_scans(self) -> dict:
        scanMetadataDict = {}
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                if "precursorMz" not in spec:
                    continue
                metadataDict = {}
                metadataDict["precursorMz"] = spec["precursorMz"][0]["precursorMz"]
                metadataDict["windowWideness"] = spec["precursorMz"][0]["windowWideness"]
                metadataDict["peaksCount"] = spec["peaksCount"]
                if "compensationVoltage" in spec:
                    CV = spec["compensationVoltage"]
                else:
                    CV = ""
                metadataDict["CV"] = CV
                scanMetadataDict[spec["num"]] = metadataDict
        return scanMetadataDict

    def pool_peaks_of_query_scans(self, scans: list) -> list:
        pooledQueryPeaks = []
        with mzxml.read(self.filePath, use_index=True) as spectra:
            for scan in scans:
                spectrum = spectra.get_by_id(scan)
                pooledQueryPeaks += create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(spectrum['m/z array'],spectrum['intensity array'],int(scan))
        return sorted(pooledQueryPeaks)

