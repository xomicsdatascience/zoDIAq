import warnings
from bisect import bisect
from zodiaq.utils import Printer
from zodiaq.slim import sum_spectra_into_single_spectrum, filter_library_keys_by_ms1_mz_precursor_values

def generate_pooled_library_and_query_spectra_by_mz_windows(libDict, queryContext):
    printer = Printer()
    ms1Count, ms2Count, totalCount = queryContext.count_number_of_ms1_and_ms2_scans_per_cycle()
    printer(f"ms1 Count: {ms1Count}, ms2Count: {ms2Count}")
    allScans = [str(i+1) for i in range(totalCount)]
    number_of_ms2_scans_per_ms1_scan = ms2Count // ms1Count

    with queryContext.get_query_file_reader() as reader:
        numMs1ScansTraversed = 0

        for i in range(ms1Count):
            matchingMs1Scans = allScans[i::ms1Count + ms2Count]
            print(matchingMs1Scans)
            ms1MergedSpectra = queryContext.pool_peaks_of_query_scans(matchingMs1Scans, reader)
            summedMs1Spectrum = sum_spectra_into_single_spectrum(ms1MergedSpectra)   # sum ms1 scans HERE
            for j in range(ms1Count, ms1Count + number_of_ms2_scans_per_ms1_scan):
                startIndex = (i * (number_of_ms2_scans_per_ms1_scan)) + j
                matchingMs2Scans = allScans[startIndex::ms1Count + ms2Count]
                print(f'--->{matchingMs2Scans}')
                lowestWindowBound, highestWindowBound = queryContext.identify_mz_window_highest_and_lowest_bound(matchingMs2Scans, reader)
                libKeys = _find_keys_of_library_spectra_in_mz_window(lowestWindowBound, highestWindowBound,
                                                                     libDict.keys())
                libKeys = filter_library_keys_by_ms1_mz_precursor_values(summedMs1Spectrum, libKeys)  # filter by ms1 scan peaks HERE
                pooledLibraryPeaks = _pool_library_spectra_by_mz_window(libKeys, libDict)
                if len(pooledLibraryPeaks) == 0:
                    continue
                ms2MergedSpectra = queryContext.pool_peaks_of_query_scans(matchingMs2Scans, reader)
                summedMs2Spectrum = sum_spectra_into_single_spectrum(ms2MergedSpectra)  # sum ms2 scans HERE
                yield pooledLibraryPeaks, summedMs2Spectrum


def _pool_library_spectra_by_mz_window(libKeys, libDict):
    pooledLibPeaks = []
    for key in libKeys:
        pooledLibPeaks.extend(libDict[key]["peaks"])
    return sorted(pooledLibPeaks)


def _find_keys_of_library_spectra_in_mz_window(lowestWindowBound, highestWindowBound, libDictKeys):
    sortedLibDictKeys = sorted(libDictKeys)
    topIndex = bisect(sortedLibDictKeys, (highestWindowBound, "z"))
    bottomIndex = bisect(sortedLibDictKeys, (lowestWindowBound, ""))
    if topIndex == bottomIndex:
        warnings.warn(
            f"No library spectra found in the {mzWindow} m/z window. Skipping",
            Warning,
        )
    return sortedLibDictKeys[bottomIndex:topIndex]
