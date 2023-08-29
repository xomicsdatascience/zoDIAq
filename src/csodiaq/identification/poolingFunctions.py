import warnings
from bisect import bisect


def generate_pooled_library_and_query_spectra_by_mz_windows(libDict, queryContext):
    queDict = queryContext.map_query_scan_ids_to_dia_mz_windows()
    for mzWindow, scans in queDict.items():
        pooledLibraryPeaks = _pool_library_spectra_by_mz_window(mzWindow, libDict)
        if len(pooledLibraryPeaks) == 0:
            continue
        pooledQueryPeaks = queryContext.pool_peaks_of_query_scans(scans)
        yield pooledLibraryPeaks, pooledQueryPeaks


def _pool_library_spectra_by_mz_window(mzWindow, libDict):
    libKeys = _find_keys_of_library_spectra_in_mz_window(mzWindow, libDict.keys())
    pooledLibPeaks = []
    for key in libKeys:
        pooledLibPeaks.extend(libDict[key]["peaks"])
    return sorted(pooledLibPeaks)


def _find_keys_of_library_spectra_in_mz_window(mzWindow, libDictKeys):
    sortedLibDictKeys = sorted(libDictKeys)
    topMz = mzWindow[0] + mzWindow[1] / 2
    bottomMz = mzWindow[0] - mzWindow[1] / 2
    topIndex = bisect(sortedLibDictKeys, (topMz, "z"))
    bottomIndex = bisect(sortedLibDictKeys, (bottomMz, ""))
    if topIndex == bottomIndex:
        warnings.warn(
            f"No library spectra found in the {mzWindow} m/z window. Skipping",
            Warning,
        )
    return sortedLibDictKeys[bottomIndex:topIndex]
