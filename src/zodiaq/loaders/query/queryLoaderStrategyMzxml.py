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
import pandas as pd
import numpy as np
import os

class QueryLoaderStrategyMzxml(QueryLoaderStrategy):
    """
    Concrete strategy implementation of the QueryLoaderStrategy strategy class specific for
        loading mzXML query files.
    """

    def map_query_scan_ids_to_dia_mz_windows(self) -> dict:
        mzWindowToScanIdDict = defaultdict(list)
        count = 0
        precMzs = set()
        with mzxml.read(self.filePath) as spectra:
            for spec in spectra:
                scan = spec["num"]
                if "precursorMz" not in spec:
                    count = 0
                    warnings.warn(
                        precursor_mz_missing_warning_text(scan),
                        SyntaxWarning,
                    )
                    continue
                else:
                    count += 1
                precMz = spec["precursorMz"][0]["precursorMz"]
                #if precMz in precMzs:
                #    return
                precMzs.add(precMz)
                windowWidth = spec["precursorMz"][0]["windowWideness"]
                print(spec["precursorMz"])
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
        #dir = '/Users/cranneyc/Desktop/summedSpectraAnalysis'
        pooledQueryPeaks = []
        for scan in scans:
            spectrum = reader.get_by_id(scan)
            queryPeaks = create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
                spectrum["m/z array"], spectrum["intensity array"], int(scan)
            )
            pooledQueryPeaks += queryPeaks
            #df = pd.DataFrame(sorted(queryPeaks), columns=['mz','intensity','scan'])
            #df.to_csv(os.path.join(dir, f'{scan}.csv'), index=False)

        sortedPooledQueryPeaks = sorted(pooledQueryPeaks)
        return sortedPooledQueryPeaks

