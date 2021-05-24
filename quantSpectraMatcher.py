import numpy as np
from numba import njit, typeof
import matplotlib.pyplot as pyplot
import csv


@njit
def approx(x, y, ppmTol):
    if x==y: return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)

@njit
def ppm_offset(mz, ppm):
    return (mz*-ppm)/1000000

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

@njit
def too_few_counts(minMatch, lightCount, heavyCount, minRank):
    if minMatch: return (lightCount < minMatch and heavyCount < minMatch)
    else: return (minRank > 2)

@njit
def tag_sparse_peak_matches(matchLibTags, lightHeavyMarks, libraryIntensityRank, matchQueTags, minMatch):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    ranks = [libraryIntensityRank[0]]
    if lightHeavyMarks[0]: heavyCount = 1; lightCount = 0
    else: heavyCount = 0; lightCount = 1
    count = 1
    remove = []
    length = len(matchLibTags)
    for i in range(1,length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            if too_few_counts(minMatch, lightCount, heavyCount, min(ranks)): remove.extend([i-j for j in range(1,count+1)])
            if lightHeavyMarks[i]: heavyCount = 1; lightCount = 0
            else: heavyCount = 0; lightCount = 1
            count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
            ranks = [libraryIntensityRank[i]]
        else:
            if lightHeavyMarks[i]: heavyCount += 1
            else: lightCount += 1
            count += 1
            ranks.append(libraryIntensityRank[i])
    if too_few_counts(minMatch, lightCount, heavyCount, min(ranks)): remove.extend([length-j for j in range(1,count+1)])
    return remove

def insufficient_matches(currentRankTable, minMatch, smallestIntensityValue):
    if minMatch:
        smallestIntensityOccurrences = 30 - np.count_nonzero(currentRankTable == smallestIntensityValue, axis = 0)
        return max(smallestIntensityOccurrences[0]) < minMatch
    return np.all(currentRankTable[:3]==currentRankTable[0][0])


def calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue, printer):
    if insufficient_matches(currentRankTable, minMatch, smallestIntensityValue): return np.nan
    log2Ratios = []
    for row in currentRankTable:
        if row[0] != row[1]:
            log2Ratios.append(np.log2(row[1]/row[0]))
            if printer: print('heavy: ' + str(row[1]) + ', light: '+ str(row[0]))
    if printer: print('\n\n')
    if ratioType=='median': return np.median(log2Ratios)
    elif ratioType=='mean': return np.mean(log2Ratios)

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


class QuantSpectraMatcher:
    def __init__(self):
        self.libraryLightHeavyMark = np.array([],dtype=int)
        self.libraryIntensityRank = np.array([],dtype=int)
        self.libraryPeptides = np.array([],dtype=np.str_)
        self.queryTags = np.array([],dtype=int)
        self.libraryIntensities = np.array([],dtype=float)
        self.queryIntensities = np.array([],dtype=float)
        self.ppmMatches = np.array([],dtype=float)

    def compare_spectra(self, pooledLibSpectra, pooledQueSpectra, tolerance, minMatch):
        libMzs, libIntensities, libFullTags = list(map(list, zip(*pooledLibSpectra)))
        queMzs, queIntensities, queTags = list(map(list, zip(*pooledQueSpectra)))

        libraryTagIndices, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = find_matching_peaks(libMzs, libIntensities, list(range(len(libFullTags))), queMzs, queIntensities, queTags, tolerance)

        if len(libraryTagIndices) != 0:
            self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides = list(map(list, zip(*[(libFullTags[x][0],libFullTags[x][1],libFullTags[x][2]) for x in libraryTagIndices])))
            self.sort_matches_by_tags()
            self.remove_sparse_matches(minMatch)

    def sort_matches_by_tags(self):
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [np.array(x) for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i1 = self.queryTags.argsort(kind='mergesort')
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i1] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i2 = self.libraryPeptides.argsort(kind='mergesort')
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i2] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def remove_sparse_matches(self, minMatch):
        remove_indices = tag_sparse_peak_matches(
            self.libraryPeptides,
            self.libraryLightHeavyMark,
            self.libraryIntensityRank,
            self.queryTags,
            minMatch)
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [np.delete(x, remove_indices) for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def extend_all_spectra(self, spectraMatch):
        self.libraryLightHeavyMark = np.append(self.libraryLightHeavyMark, spectraMatch.libraryLightHeavyMark)
        self.libraryIntensityRank = np.append(self.libraryIntensityRank, spectraMatch.libraryIntensityRank)
        self.libraryPeptides = np.append(self.libraryPeptides, spectraMatch.libraryPeptides)
        self.queryTags = np.append(self.queryTags, spectraMatch.queryTags)
        self.libraryIntensities = np.append(self.libraryIntensities, spectraMatch.libraryIntensities)
        self.queryIntensities = np.append(self.queryIntensities, spectraMatch.queryIntensities)
        self.ppmMatches = np.append(self.ppmMatches, spectraMatch.ppmMatches)

    def filter_by_corrected_ppm_window(self, corrected, histFile, minMatch):
        offset, tolerance = find_offset_tolerance(self.ppmMatches, histFile, corrected)
        lowend = offset-tolerance
        highend = offset+tolerance
        ppmIndices = np.where((self.ppmMatches>lowend)*(self.ppmMatches<highend))[0]
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[ppmIndices] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        self.remove_sparse_matches(minMatch)

    def determine_ratios(self, ratioDict, scanToNoiseIntensityCutoffDict, ratioType, minMatch):

        curLibTag = self.libraryPeptides[0]
        curQueTag = self.queryTags[0]
        smallestIntensityValue = scanToNoiseIntensityCutoffDict[curQueTag]

        currentRankTable = np.full((30,2),smallestIntensityValue) # 30 should be a different variable
        if self.libraryLightHeavyMark[0]: currentRankTable[self.libraryIntensityRank[0]][1] = self.queryIntensities[0]
        else: currentRankTable[self.libraryIntensityRank[0]][0] = self.queryIntensities[0]


        length = len(self.libraryPeptides)
        for i in range(1,length):
            if self.libraryPeptides[i] != curLibTag or self.queryTags[i] != curQueTag:
                smallestIntensityValue = scanToNoiseIntensityCutoffDict[curQueTag]
                if curLibTag == 'C+57.0215DSSPDSAEDVRK': printer=1
                else: printer=0
                ratio = calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue, printer) # needs ratio type, probably minMatch too
                ratioDict[curQueTag, curLibTag] = ratio
                currentRankTable = np.full((30,2),smallestIntensityValue)
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]
                curLibTag = self.libraryPeptides[i]
                curQueTag = self.queryTags[i]
            else:
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]

        smallestIntensityValue = scanToNoiseIntensityCutoffDict[self.queryTags[-1]]
        ratio = calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue, printer) # needs ratio type, probably minMatch too
        ratioDict[curQueTag, curLibTag] = ratio

        return ratioDict
