import numpy as np
from numba import njit
import matplotlib.pyplot as pyplot
from timeit import default_timer as timer
from datetime import timedelta


@njit
def approx(x, y, ppmTol):
    if x==y: return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)

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
        if not approx(libMzs[i],queMzs[j], ppmTol):
            if libMzs[i] > queMzs[j]:
                j += 1
                continue
            if libMzs[i] < queMzs[j]: i += 1; continue
        p = i + 0
        while (p < lenLib):
            ppm = approx(libMzs[p], queMzs[j], ppmTol)
            if p==lenLib or not ppm: break
            matchLibTags.append(libTags[p])
            matchLibIntensities.append(libIntensities[p])
            matchQueTags.append(queTags[j])
            matchQueIntensities.append(queIntensities[j])
            ppmMatches.append(ppm)
            p += 1
        j += 1
    return matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches

def find_offset_tolerance(data, histFile, stdev, mean=True):
    if len(data)==0: return 0, 10
    hist, bins = np.histogram(data, bins=200)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    # offset is calculated as the mean or median value of the provided data, though this is overwritten if no corrected standard deviation is provided
    if mean: offset = sum(data)/len(data)
    else: offset = data[len(data)//2]

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
            if hist[i] < noise: min_i = i; break
        for i in range(index_max, len(hist)):
            if hist[i] < noise: max_i = i; break

        offset = center[index_max]
        if index_max - min_i >= max_i - index_max: tolerance = offset - center[min_i]
        else: tolerance = center[max_i] - offset

    # if a histogram file is provided, it is created, with offset (black) and tolerance (red) lines drawn for reference
    if histFile:
        pyplot.clf()
        pyplot.bar(center, hist, align='center', width=width)
        pyplot.axvline(x=offset, color='black', linestyle = 'dashed', linewidth=4)
        pyplot.axvline(x=offset-tolerance, color='red', linestyle = 'dashed', linewidth=4)
        pyplot.axvline(x=offset+tolerance, color='red', linestyle = 'dashed', linewidth=4)
        pyplot.suptitle('offset: '+str(offset) + ', tolerance: '+str(tolerance))
        pyplot.savefig(histFile)
    return offset, tolerance

def print_milestone(text):
    print(text)
    print(str(timedelta(seconds=timer())),flush=True)

def format_spectra_for_pooling(spectrum, scanNumber, sqrt=True):
    scanNumber = int(scanNumber)
    if sqrt: intensity = [x**0.5 for x in spectrum['intensity array']]
    else: intensity = [x for x in spectrum['intensity array']]
    peakIDs = [scanNumber for x in range(spectrum['peaksCount'])]
    return list(zip(spectrum['m/z array'],intensity,peakIDs))

def calculate_heavy_mz(seq, mz, z):
    hK=8.014199 ## mass of heavy lysine
    hR=10.00827 ## mass of heavy arg

    nK = seq.count('K')
    nR = seq.count('R')
    heavyMz = mz + (nK*hK)/z + (nR*hR)/z
    return heavyMz

def approx_list(x, l, ppmTol=10):
    for i in range(len(l)):
        if approx(x, l[i], ppmTol): return i
    return -1
