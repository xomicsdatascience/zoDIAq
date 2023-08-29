from csodiaq.loaders.query.queryLoaderStrategy import (
    QueryLoaderStrategy,
    precursor_mz_missing_warning_text,
)
from csodiaq.loaders.library.libraryLoaderStrategy import (
    create_peaks_from_mz_intensity_lists_and_csodiaq_key_id,
)
from collections import defaultdict
from pyteomics import mzxml
import warnings

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

    def extract_metadata_from_query_scans(self) -> dict:
        scanMetadataDict = {}
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                metadataDict = {}
                metadataDict["precursorMz"] = spec["precursorMz"][0]["precursorMz"]
                metadataDict["windowWidth"] = spec["precursorMz"][0]["windowWideness"]
                metadataDict["peaksCount"] = spec["peaksCount"]
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

    def pool_peaks_of_query_scans(self, scans: list) -> list:
        pooledQueryPeaks = []
        with mzxml.read(self.filePath, use_index=True) as spectra:
            for scan in scans:
                spectrum = spectra.get_by_id(scan)
                pooledQueryPeaks += (
                    create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(
                        spectrum["m/z array"], spectrum["intensity array"], int(scan)
                    )
                )
        return sorted(pooledQueryPeaks)
