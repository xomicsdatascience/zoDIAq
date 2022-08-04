import numpy as np
from numba import njit
from . import spectra_matcher_functions as smf

@njit
def too_few_counts(minMatch, lightCount, heavyCount, ranks):
    if minMatch: return len(set(ranks)) < minMatch
    else: return (min(ranks) > 2)

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
            if too_few_counts(minMatch, lightCount, heavyCount, ranks): remove.extend([i-j for j in range(1,count+1)])
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
    if too_few_counts(minMatch, lightCount, heavyCount, ranks): remove.extend([length-j for j in range(1,count+1)])
    return remove

def insufficient_matches(currentRankTable, minMatch, smallestIntensityValue):
    if minMatch:
        abbreviatedRankTable = currentRankTable[~np.all(currentRankTable == smallestIntensityValue, axis=1)]
        return len(abbreviatedRankTable) < minMatch
    return np.all(currentRankTable[:3]==smallestIntensityValue)


def calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue):
    if insufficient_matches(currentRankTable, minMatch, smallestIntensityValue): return np.nan
    log2Ratios = []
    for row in currentRankTable:
        if row[0] != row[1]:
            log2Ratios.append(np.log2(row[1]/row[0]))
    if ratioType=='median': return np.median(log2Ratios)
    elif ratioType=='mean': return np.mean(log2Ratios)


class QuantificationSpectraMatcher:
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

        libraryTagIndices, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = smf.find_matching_peaks(libMzs, libIntensities, list(range(len(libFullTags))), queMzs, queIntensities, queTags, tolerance)

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
        offset, tolerance = smf.find_offset_tolerance(self.ppmMatches, histFile, corrected)
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
                ratio = calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue) # needs ratio type, probably minMatch too
                ratioDict[curQueTag, curLibTag] = ratio
                curLibTag = self.libraryPeptides[i]
                curQueTag = self.queryTags[i]
                smallestIntensityValue = scanToNoiseIntensityCutoffDict[curQueTag]
                currentRankTable = np.full((30,2),smallestIntensityValue)
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]

            else:
                if self.libraryLightHeavyMark[i]: currentRankTable[self.libraryIntensityRank[i]][1] = self.queryIntensities[i]
                else: currentRankTable[self.libraryIntensityRank[i]][0] = self.queryIntensities[i]

        smallestIntensityValue = scanToNoiseIntensityCutoffDict[self.queryTags[-1]]
        ratio = calculate_ratio(currentRankTable, ratioType, minMatch, smallestIntensityValue) # needs ratio type, probably minMatch too
        ratioDict[curQueTag, curLibTag] = ratio

        return ratioDict
