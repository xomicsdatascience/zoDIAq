import numpy as np
from zodiaq.identification.matchingFunctions import calculate_parts_per_million_relative_difference

def passthrough_sum_spectra_into_single_spectrum(mixedSpectra, type='ms2'):
    return mixedSpectra

def sum_spectra_into_single_spectrum(mixedSpectra, type='ms2'):
    mzs = np.array([x[0] for x in mixedSpectra])
    threshold = _identify_MAD_from_mz_values(mzs)
    clusters = _identify_clusters_using_threshold(mzs, threshold)
    mixedSpectraArray = np.array(mixedSpectra)
    summedMzs = np.array([np.average(mixedSpectraArray[clusters == cat][:, 0]) for cat in np.unique(clusters)])
    summedIntensities = np.array([np.sum(mixedSpectraArray[clusters == cat][:, 1]) for cat in np.unique(clusters)])
    minTag = min([x[-1] for x in mixedSpectra])
    return list(zip(summedMzs, summedIntensities, [minTag]*len(summedMzs)))

def _identify_MAD_from_mz_values(mz):
    mzDiff = _calculate_differences_of_adjacent_mz_values(mz)
    return _calculate_MAD(mzDiff)

def _calculate_differences_of_adjacent_mz_values(mz):
    mz = np.sort(mz)
    mzDiff = calculate_parts_per_million_relative_difference(mz[1:], mz[:-1])
    return np.append(0, mzDiff)

def _calculate_MAD(values):
    return np.median(np.abs(values - np.median(values)))

def _identify_clusters_using_threshold(mz, threshold):
    mzDiff = _calculate_differences_of_adjacent_mz_values(mz)
    split = (mzDiff >= threshold).astype(int)
    return np.cumsum(split)
