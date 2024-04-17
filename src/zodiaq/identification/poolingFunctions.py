import warnings
from bisect import bisect
from zodiaq.utils import Printer
from zodiaq.identification.matchingFunctions import (
    calculate_parts_per_million_relative_difference,
    is_within_tolerance,
)
import numpy as np

def generate_pooled_library_and_query_spectra_by_mz_windows(libDict, queryContext):
    printer = Printer()
    queDict = queryContext.map_query_scan_ids_to_dia_mz_windows()
    ms1Dict = queryContext._strategy.find_all_ms1_mz_peaks_per_scan()
    numScans, numMs2Gaps = _find_ms1_scan_span_and_number_of_ms2_cycles(sorted(ms1Dict.keys()))
    printer(f"Total number of m/z windows: {len(queDict.keys())}")
    numWindowsTraversed = 0
    with queryContext.get_query_file_reader() as reader:
        for mzWindow, scans in queDict.items():
            scans = sorted([int(scan) for scan in scans])
            scans = [str(scan) for scan in scans]
            #print(scans)
            for i in range(0, len(scans)):
                ms2Scan = scans[i]
                ms1Scan = _find_corresponding_number(int(ms2Scan), numScans, numMs2Gaps)
                ms1MzValues = ms1Dict[ms1Scan] #np.array([1])
                numWindowsTraversed += 1
                printer(
                    f"Checkpoint: {numWindowsTraversed} / {len(queDict.keys())*len(scans)} windows traversed",
                    checkPoint=True,
                )
                pooledLibraryPeaks = _pool_library_spectra_by_mz_window(mzWindow, ms1MzValues, libDict)
                if len(pooledLibraryPeaks) == 0:
                    continue
                pooledQueryPeaks = queryContext.pool_peaks_of_query_scans([ms2Scan], reader)
                yield pooledLibraryPeaks, pooledQueryPeaks

def _find_corresponding_number(num, numScans=15, numMs2Gaps=3):
    baseline = ((num-1) // (numScans*numMs2Gaps)) * (numScans*numMs2Gaps) + 1
    diff = (num - baseline - numScans) // (numMs2Gaps - 1)
    return baseline + diff

def _find_ms1_scan_span_and_number_of_ms2_cycles(lst):
    for i in range(len(lst)):
        if lst[i] != lst[0] + i:
            return lst[i-1], (lst[i] - lst[0]) // lst[i-1]
    return None, None  # If no gap is found

def _pool_library_spectra_by_mz_window(mzWindow, ms1MzValues, libDict):
    #print('*'*20)
    libKeys = _find_keys_of_library_spectra_in_mz_window(mzWindow, libDict.keys())
    #print(len(libKeys))
    libKeys = _filter_lib_keys_by_ms1_values(libKeys, ms1MzValues)
    #print(len(libKeys))
    pooledLibPeaks = []
    for key in libKeys:
        pooledLibPeaks.extend(libDict[key]["peaks"])
    return sorted(pooledLibPeaks)


def _filter_lib_keys_by_ms1_values(libKeys, ms1MzValues):
    filteredKeys = [(mz, peptide) for mz, peptide in libKeys if _is_lib_mz_found_in_ms1(mz, ms1MzValues)]
    return filteredKeys

def _is_lib_mz_found_in_ms1(libMz, ms1MzValues, tolerance = 0.01):
    tempList = ms1MzValues[np.abs(libMz - ms1MzValues) <= tolerance]
    return any([is_within_tolerance(calculate_parts_per_million_relative_difference(libMz,ms1Mz),30) for ms1Mz in tempList])



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
