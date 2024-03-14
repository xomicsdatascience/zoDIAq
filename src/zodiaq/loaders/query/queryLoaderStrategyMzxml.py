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
                windowWidth = spec["precursorMz"][0]["windowWideness"]
                mzWindowToScanIdDict[precMz, windowWidth].append(scan)
        keys = sorted(dict(mzWindowToScanIdDict))
        d = dict(mzWindowToScanIdDict)
        data = []
        for i, window in enumerate(keys):
            scans = d[window]
            data.extend([(i,int(scan)) for scan in sorted(scans)])
        df = pd.DataFrame(data, columns=['cycle','scan'])
        df.to_csv('/Users/cranneyc/Desktop/cycle_to_scan_delete.csv', index=False)
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
        dir = '/Users/cranneyc/Desktop/summedSpectraAnalysis'
        pooledQueryPeaks = []
        for scan in scans:
            spectrum = reader.get_by_id(scan)
            queryPeaks = create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
                spectrum["m/z array"], spectrum["intensity array"], int(scan)
            )
            pooledQueryPeaks += queryPeaks
            df = pd.DataFrame(sorted(queryPeaks), columns=['mz','intensity','scan'])
            df.to_csv(os.path.join(dir, f'{scan}.csv'), index=False)

        sortedPooledQueryPeaks = sorted(pooledQueryPeaks)
        sortedPooledQueryPeaks = sum_pooled_query_peaks(sortedPooledQueryPeaks)
        df = pd.DataFrame(sortedPooledQueryPeaks, columns=['mz', 'intensity', 'scan'])
        df.to_csv(os.path.join(dir, f'{"_".join(scans[:10])}_etc.csv'), index=False)
        return 0
        return sortedPooledQueryPeaks

def sum_pooled_query_peaks(sortedPooledQueryPeaks):
    df = pd.DataFrame(sortedPooledQueryPeaks, columns=['mz','intensity','id'])
    minId = min(df['id'])
    df['sharedPeaks'] = identify_metabolites_using_threshold(df['mz'], threshold=0.07)
    value_counts = df['sharedPeaks'].value_counts()
    #mask = df['sharedPeaks'].map(value_counts) >= 2
    #df = df[mask]
    summedDf = df.groupby('sharedPeaks').agg({'mz': 'mean', 'intensity':'sum'})
    summedDf['id'] = [minId for i in range(len(summedDf))]
    min_value = summedDf['intensity'].min()
    max_value = summedDf['intensity'].max()
    summedDf['intensity'] = 1000000 * (summedDf['intensity'] - min_value) / (max_value - min_value)
    return sorted([tuple(x) for x in summedDf.values.tolist()])


def identify_metabolites_using_threshold(mz, threshold):
    mzDiff = calculate_differences_of_adjacent_mz_values(mz)
    split = (mzDiff >= threshold).astype(int)
    return np.cumsum(split)

def identify_MAD_from_mz_values(mz):
    mzDiff = calculate_differences_of_adjacent_mz_values(mz)
    return calculate_MAD(mzDiff)

def calculate_parts_per_million_relative_difference(referenceMz, targetMz):
    return (referenceMz - targetMz) * (1e6) / referenceMz

def calculate_differences_of_adjacent_mz_values(mz):
    mz = np.sort(mz)
    mzDiff = calculate_parts_per_million_relative_difference(mz[1:], mz[:-1])
    return np.append(0, mzDiff)

def calculate_MAD(values):
    return np.median(np.abs(values - np.median(values)))

