import numpy as np
from numba import njit
import matplotlib.pyplot as pyplot
from timeit import default_timer as timer
from datetime import timedelta

@njit
def calc_ppm_diff(ref_val: float, val: float):
    '''
    Calculates the difference and returns the value in parts per million (ppm) of ref_val.
    Parameters
    ----------
    ref_val : float
        Value used as reference for the difference and ppm.
    val : float
        Value to use for the difference.
    Returns
    -------
    float
        Difference between ref_val and val in ppm.
    '''
    return (ref_val-val)*(1e6)/ref_val

@njit
def is_within_ppm_tol(ref_val: float, val: float, ppm_tol: float):
    '''
    Verifies that the input values are within ppm_tol (±ppm_tol) of each other.
    Parameters
    ----------
    ref_val : float
        First value.
    val : float
        Second value.
    ppm_tol : float
        Tolerance value in ppm.

    Returns
    -------
    bool
        True if the input values are within tolerance. False otherwise
    '''
    ppm_diff = calc_ppm_diff(ref_val, val)
    return abs(ppm_diff) <= ppm_tol

@njit
def find_matching_peaks(libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, ppmTol):
    lenLib = len(libMzs)
    lenQue = len(queMzs)
    matchLibTags = []
    matchLibIntensities = []
    matchQueTags = []
    matchQueIntensities = []
    ppmMatches = []
    i, j = 0, 0
    while i < lenLib and j < lenQue:
        if not is_within_ppm_tol(libMzs[i], queMzs[j], ppmTol):
            if libMzs[i] > queMzs[j]:
                j += 1
                continue
            if libMzs[i] < queMzs[j]:
                i += 1
                continue
        p = i + 0
        while (p < lenLib):
            ppm = calc_ppm_diff(libMzs[p], queMzs[j])
            ppm_within_tol = is_within_ppm_tol(libMzs[p], queMzs[j], ppmTol)
            if not ppm_within_tol:  # if library and que aren't within tol, continue to next match
                break
            # Store matches
            matchLibTags.append(libTags[p])
            matchLibIntensities.append(libIntensities[p])
            matchQueTags.append(queTags[j])
            matchQueIntensities.append(queIntensities[j])
            ppmMatches.append(ppm)
            p += 1
        j += 1
    return matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches


def find_offset_tolerance(data, histFile, stdev, mean=True):
    if len(data) == 0:
        return 0, 10
    hist, bins = np.histogram(data, bins=200)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    # offset is calculated as the mean or median value of the provided data, though this is overwritten if no corrected standard deviation is provided
    if mean:
        offset = sum(data)/len(data)
    else:
        offset = data[len(data)//2]

    # If a corrected standard deviation is provided, it is used to set the tolerance. If not, see below conditional statement
    tolerance = np.std(data)*stdev

    # If no corrected standard deviation is provided, one is customized.
    #  Offset is considered the highest point of the histogram,
    #  tolerance the range necessary before the peaks are the same size as the noise on either end.
    if not stdev:
        index_max = max(range(len(hist)), key=hist.__getitem__)
        histList = list(hist)
        noise = np.mean(hist[:10] + hist[-10:])
        min_i = 0
        max_i = len(hist)-1
        for i in range(index_max, 0, -1):
            if hist[i] < noise:
                min_i = i
                break
        for i in range(index_max, len(hist)):
            if hist[i] < noise:
                max_i = i
                break

        offset = center[index_max]
        if index_max - min_i >= max_i - index_max:
            tolerance = offset - center[min_i]
        else:
            tolerance = center[max_i] - offset

    # if a histogram file is provided, it is created, with offset (black) and tolerance (red) lines drawn for reference
    if histFile:
        pyplot.clf()
        pyplot.bar(center, hist, align='center', width=width)
        pyplot.axvline(x=offset, color='black',
                       linestyle='dashed', linewidth=4)
        pyplot.axvline(x=offset-tolerance, color='red',
                       linestyle='dashed', linewidth=4)
        pyplot.axvline(x=offset+tolerance, color='red',
                       linestyle='dashed', linewidth=4)
        pyplot.suptitle('offset: '+str(offset) +
                        ', tolerance: '+str(tolerance))
        pyplot.savefig(histFile)
    return offset, tolerance


def print_milestone(text):
    print(text, flush=True)
    # print(str(timedelta(seconds=timer())), flush=True)
    return

def format_spectra_for_pooling(spectrum, scanNumber, sqrt=True):
    scanNumber = int(scanNumber)
    if sqrt:
        intensity = [x**0.5 for x in spectrum['intensity array']]
    else:
        intensity = [x for x in spectrum['intensity array']]
    peakIDs = [scanNumber for x in range(spectrum['peaksCount'])]
    return list(zip(spectrum['m/z array'], intensity, peakIDs))


def calculate_heavy_mz(seq, mz, z):
    hK = 8.014199  # mass of heavy lysine
    hR = 10.00827  # mass of heavy arg

    nK = seq.count('K')
    nR = seq.count('R')
    heavyMz = mz + (nK*hK)/z + (nR*hR)/z
    return heavyMz


def approx_list(x, l, ppmTol=10):
    for i in range(len(l)):
        if is_within_ppm_tol(x, l[i], ppm_tol=ppmTol):
            return i
    return -1
