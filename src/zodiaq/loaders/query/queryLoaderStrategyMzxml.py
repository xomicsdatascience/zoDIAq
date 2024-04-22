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

    def map_query_scan_ids_to_dia_mz_windows(self) -> dict:
        mzWindowToScanIdDict = defaultdict(list)
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                scan = spec["num"]
                if "precursorMz" not in spec:
                    warnings.warn(
                        precursor_mz_missing_warning_text(scan),
                        SyntaxWarning,
                    )
                    continue
                precMz = spec["precursorMz"][0]["precursorMz"]
                windowWidth = spec["precursorMz"][0]["windowWideness"]
                mzWindowToScanIdDict[precMz, windowWidth].append(scan)
        return dict(mzWindowToScanIdDict)

    def find_all_ms1_mz_peaks_per_scan(self) -> dict:
        mzWindowToMs1MzDict = defaultdict(list)
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                scan = spec["num"]
                if "precursorMz" in spec:
                    #print()
                    continue

                #print(spec)
                midIntensity = sorted(spec['intensity array'])[-10]#[-len(spec['intensity array'])//5]
                #midIntensity = 0
                #print(len(spec['intensity array']))
                mzWindowToMs1MzDict[int(scan)] += list(spec['m/z array'][spec['intensity array'] >= midIntensity])

        for key, value in mzWindowToMs1MzDict.items():
            mzWindowToMs1MzDict[key] = np.array(sorted(value))
        return dict(mzWindowToMs1MzDict)

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
