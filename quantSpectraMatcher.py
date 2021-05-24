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
def tag_sparse_peak_matches(matchLibTags, lightHeavyMarks, matchQueTags):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    if lightHeavyMarks[0]: heavyCount = 1; lightCount = 0
    else: heavyCount = 0; lightCount = 1
    count = 1
    remove = []
    length = len(matchLibTags)
    for i in range(1,length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            if lightCount < 3 and heavyCount < 3: remove.extend([i-j for j in range(1,count+1)])
            if lightHeavyMarks[i]: heavyCount = 1; lightCount = 0
            else: heavyCount = 0; lightCount = 1
            count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
        else:
            if lightHeavyMarks[i]: heavyCount += 1
            else: lightCount += 1
            count += 1
    if lightCount < 3 and heavyCount < 3: remove.extend([length-j for j in range(1,count+1)])
    return remove

def calculate_ratio(currentRankTable):
    if np.all(currentRankTable[:3]==currentRankTable[0][0]): return np.nan
    log2Ratios = []
    for row in currentRankTable:
        if row[0] != row[1]: log2Ratios.append(np.log2(row[1]/row[0]))
    return np.median(log2Ratios)

class QuantSpectraMatcher:
    def __init__(self):
        self.libraryLightHeavyMark = np.array([],dtype=int)
        self.libraryIntensityRank = np.array([],dtype=int)
        self.libraryPeptides = np.array([],dtype=np.str_)
        self.queryTags = np.array([],dtype=int)
        self.libraryIntensities = np.array([],dtype=float)
        self.queryIntensities = np.array([],dtype=float)
        self.ppmMatches = np.array([],dtype=float)

    def compare_spectra(self, pooledLibSpectra, pooledQueSpectra, tolerance):
        for x in pooledLibSpectra: print(x)
        print(len(pooledQueSpectra))
        libMzs, libIntensities, libFullTags = list(map(list, zip(*pooledLibSpectra)))
        queMzs, queIntensities, queTags = list(map(list, zip(*pooledQueSpectra)))

        libraryTagIndices, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = find_matching_peaks(libMzs, libIntensities, list(range(len(libFullTags))), queMzs, queIntensities, queTags, tolerance)
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides = list(map(list, zip(*[(libFullTags[x][0],libFullTags[x][1],libFullTags[x][2]) for x in libraryTagIndices])))
        self.sort_matches_by_tags()
        self.remove_sparse_matches()

    def sort_matches_by_tags(self):
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [np.array(x) for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i1 = self.queryTags.argsort(kind='mergesort')
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i1] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i2 = self.libraryPeptides.argsort(kind='mergesort')
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i2] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def remove_sparse_matches(self):
        remove_indices = tag_sparse_peak_matches(
            self.libraryPeptides,
            self.libraryLightHeavyMark,
            self.queryTags)
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [np.delete(x, remove_indices) for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def extend_all_spectra(self, spectraMatch):
        self.libraryLightHeavyMark = np.append(self.libraryLightHeavyMark, spectraMatch.libraryLightHeavyMark)
        self.libraryIntensityRank = np.append(self.libraryIntensityRank, spectraMatch.libraryIntensityRank)
        self.libraryPeptides = np.append(self.libraryPeptides, spectraMatch.libraryPeptides)
        self.queryTags = np.append(self.queryTags, spectraMatch.queryTags)
        self.libraryIntensities = np.append(self.libraryIntensities, spectraMatch.libraryIntensities)
        self.queryIntensities = np.append(self.queryIntensities, spectraMatch.queryIntensities)
        self.ppmMatches = np.append(self.ppmMatches, spectraMatch.ppmMatches)

    def filter_by_corrected_ppm_window(self, corrected, histFile):
        offset, tolerance = find_offset_tolerance(self.ppmMatches, histFile, corrected)
        lowend = offset-tolerance
        highend = offset+tolerance
        ppmIndices = np.where((self.ppmMatches>lowend)*(self.ppmMatches<highend))[0]
        self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[ppmIndices] for x in [self.libraryLightHeavyMark, self.libraryIntensityRank, self.libraryPeptides, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        self.remove_sparse_matches()

    def determine_ratios(self, scanToNoiseIntensityCutoffDict):

        ratioDict = {}

        curLibTag = self.libraryPeptides[0]
        curQueTag = self.queryTags[0]

        currentRankTable = np.full((30,2),smallestIntensityValue) # 30 should be a different variable
        if self.libraryLightHeavyMark[0]: currentRankTable[self.libraryIntensityRank[0]][1] = self.queryIntensities[0]
        else: currentRankTable[self.libraryIntensityRank[0]][0] = self.queryIntensities[0]


        length = len(self.libraryPeptides)
        for i in range(1,length):
            if self.libraryPeptides[i] != curLibTag or self.queryTags[i] != curQueTag:
                smallestIntensityValue = scanToNoiseIntensityCutoffDict[curQueTag]
                ratio = calculate_ratio(currentRankTable, smallestIntensityValue) # needs ratio type, probably minMatch too
                ratioDict[curQueTag, curLibTag] = ratio
                currentRankTable = np.full((30,2),smallestIntensityValue)
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]
            else:
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]

        smallestIntensityValue = scanToNoiseIntensityCutoffDict[-1]
        ratio = calculate_ratio(currentRankTable, smallestIntensityValue) # needs ratio type, probably minMatch too
        ratioDict[curQueTag, curLibTag] = ratio

        return ratioDict
