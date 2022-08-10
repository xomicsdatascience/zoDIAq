import numpy as np
from numba import njit
import csv
from . import spectra_matcher_functions as smf


@njit
def cosine_similarity(AB, A, B):
    magnitude = (A**0.5) * (B**0.5)  # (sqrt(sum(A^2))*sqrt(sum(B^2)))
    return (AB / magnitude if magnitude else 0)


@njit
def tag_sparse_peak_matches_and_generate_scores(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities):

    # Because the tags are sorted together, each spectrum-spectrum tag pair is grouped together. As soon as the pair changes, that means you've gone to to a new pair.
    #  We keep track of the current tags to determine if they have changed.
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]

    # count keeps track of the number of peak matches in each spectra match
    count = 1

    # keeping track of variables necessary to calculate cosine score and, consequentally, the macc score
    AB = matchLibIntensities[0]*matchQueIntensities[0]
    A = matchLibIntensities[0]**2
    B = matchQueIntensities[0]**2

    # initializing the return values
    remove = []
    returnMaccScores = []

    # looping through every peak match
    length = len(matchLibTags)
    for i in range(1, length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:

            # if the count is 3 or higher, calculate the macc score and add it
            if count > 2:
                test = 0
                cosine = cosine_similarity(AB, A, B)
                macc = (count**(1/5))*cosine
                # (sqrt(sum(A^2))*sqrt(sum(B^2)))
                magnitude = (A**0.5) * (B**0.5)
                returnMaccScores.extend([macc]*count)

            # if the count is 1 or 2, remove the peak match(es) represented by the indices
            else:
                remove.extend([i-j for j in range(1, count+1)])

            # reset the variables
            count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
            AB = matchLibIntensities[i]*matchQueIntensities[i]
            A = matchLibIntensities[i]**2
            B = matchQueIntensities[i]**2

        # if its the same spectra tag pair, contribute to the variables used to calculate the cosine score
        else:
            AB += matchLibIntensities[i]*matchQueIntensities[i]
            A += matchLibIntensities[i]**2
            B += matchQueIntensities[i]**2
            count += 1

    # the middle of the loop is repeated for the last spectra key pairing as well
    if count > 2:
        test = 0
        cosine = cosine_similarity(AB, A, B)
        macc = (count**(1/5))*cosine
        returnMaccScores.extend([macc]*count)
    else:
        remove.extend([length-j for j in range(1, count+1)])

    return remove, returnMaccScores


@njit
def reduce_duplicate_rows(matchLibTags, matchQueTags, maccScores, decoys):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    returnMaccs = [maccScores[0]]
    returnDecoys = [decoys[0]]
    length = len(matchLibTags)
    for i in range(1, length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            returnMaccs.append(maccScores[i])
            returnDecoys.append(decoys[i])
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
    return returnMaccs, returnDecoys


@njit
def calculate_score_fdr_cutoff(scores, decoys):  # ***NOTE***make return type
    # initializing the two return values at 0
    numDecoys = 0

    # for every row in the dataframe
    count = 0
    for i in range(len(scores)):
        if decoys[i]:
            numDecoys += 1
        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(count+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if count < 1/0.01:
                return -1
            return scores[i-1]
        count += 1

    return scores[-1]


@njit
def collect_relevant_ppm_values(ppmMatches, maccScores, maccCutoff):
    ppms = []
    length = len(ppmMatches)
    for i in range(1, length):
        if maccScores[i] >= maccCutoff:
            ppms.append(ppmMatches[i])
    return ppms


class IdentificationSpectraMatcher:
    def __init__(self):
        self.libraryTags = np.array([], dtype=int)
        self.queryTags = np.array([], dtype=int)
        self.libraryIntensities = np.array([], dtype=float)
        self.queryIntensities = np.array([], dtype=float)
        self.ppmMatches = np.array([], dtype=float)
        self.scores = np.array([], dtype=float)
        self.decoys = np.array([], dtype=float)

    def compare_spectra(self, pooledLibSpectra, pooledQueSpectra, tolerance, idToDecoyDict):
        # NOTE: Numba is remarkably slow with 2d numpy arrays, so data is tracked in 1d arrays of the same length until that is updated
        libMzs, libIntensities, libTags = list(
            map(list, zip(*pooledLibSpectra)))
        queMzs, queIntensities, queTags = list(
            map(list, zip(*pooledQueSpectra)))
        libraryTags, libraryIntensities, queryTags, queryIntensities, ppmMatches = smf.find_matching_peaks(
            libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, tolerance)
        if len(libraryTags) == 0:
            return None
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = np.array(libraryTags, dtype=int), np.array(
            libraryIntensities, dtype=float), np.array(queryTags, dtype=int), np.array(queryIntensities, dtype=float), np.array(ppmMatches, dtype=float)
        self.sort_matches_by_tags()
        self.remove_sparse_matches_and_generate_scores()
        self.decoys = [idToDecoyDict[x] for x in self.libraryTags]

    def sort_matches_by_tags(self):
        self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.ppmMatches = [np.array(
            x) for x in [self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.ppmMatches]]
        i1 = self.queryTags.argsort(kind='mergesort')
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [
            x[i1] for x in [self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i2 = self.libraryTags.argsort(kind='mergesort')
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [
            x[i2] for x in [self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def remove_sparse_matches_and_generate_scores(self):
        remove_indices, self.scores = tag_sparse_peak_matches_and_generate_scores(
            self.libraryTags,
            self.libraryIntensities,
            self.queryTags,
            self.queryIntensities)
        self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.ppmMatches = [np.delete(
            x, remove_indices) for x in [self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.ppmMatches]]
        if len(self.decoys) > 0:
            self.decoys = np.delete(self.decoys, remove_indices)

    def filter_by_corrected_ppm_window(self, corrected, scoreCutoff, histFile):

        ppms = collect_relevant_ppm_values(
            self.ppmMatches, self.scores, scoreCutoff)  # requires loop
        offset, tolerance = smf.find_offset_tolerance(ppms, histFile, stdev=0)

        lowend = offset-tolerance
        highend = offset+tolerance
        ppmIndices = np.where((self.ppmMatches > lowend)
                              * (self.ppmMatches < highend))[0]
        self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.decoys = [
            x[ppmIndices] for x in [self.libraryTags, self.queryTags, self.libraryIntensities, self.queryIntensities, self.decoys]]
        self.remove_sparse_matches_and_generate_scores()

    def find_score_fdr_cutoff(self):
        scores, decoys = reduce_duplicate_rows(
            self.libraryTags, self.queryTags, self.scores, self.decoys)
        scores = np.array(scores)
        decoys = np.array(decoys)
        i1 = (-scores).argsort()
        scores = scores[i1]
        decoys = decoys[i1]
        return calculate_score_fdr_cutoff(scores, decoys)

    def extend_all_spectra(self, spectraMatch):
        self.libraryTags = np.append(
            self.libraryTags, spectraMatch.libraryTags)
        self.libraryIntensities = np.append(
            self.libraryIntensities, spectraMatch.libraryIntensities)
        self.queryTags = np.append(self.queryTags, spectraMatch.queryTags)
        self.queryIntensities = np.append(
            self.queryIntensities, spectraMatch.queryIntensities)
        self.ppmMatches = np.append(self.ppmMatches, spectraMatch.ppmMatches)
        self.scores = np.append(self.scores, spectraMatch.scores)
        self.decoys = np.append(self.decoys, spectraMatch.decoys)

    def write_output(self, outFile, expSpectraFile, scoreCutoff, queValDict, idToKeyDict, lib):
        columns = [
            'fileName',  # Name of the query spectra file.
            # Scan number, corresponding to scans in the query spectra file.
            'scan',
            # precursor m/z for query spectrum. Column 'windowWideness' corresponds to this value.
            'MzEXP',
            # Peptide sequence for the library spectrum corresponding to this row.
            'peptide',
            # Protein name the peptide corresponds to, also derived from the library spectrum corresponding to this row.
            'protein',
            # precursor m/z for the library spectrum corresponding to this row.
            'MzLIB',
            # precursor charge for the library spectrum corresponding to this row.
            'zLIB',
            # Cosine score comparing the library spectrum corresponding to this row with the query spectrum.
            'cosine',
            # Title - corresponds to the column "transition_group_id," a library spectrum identifier.
            'name',
            'Peak(Query)',  # The number of peaks in the query spectrum.
            'Peaks(Library)',  # The number of peaks in the library spectrum.
            # The number of peaks that matched between query spectrum/library spectrum.
            'shared',
            'ionCount',  # Sum of query spectrum intensities, excluding possible duplicates
            # The compensation voltage of the query spectrum.
            'CompensationVoltage',
            # width of m/z that was captured in the query spectrum around MzEXP.
            'totalWindowWidth',
            # score unique to CsoDIAq, the fifth root of the number of matches ('shared') multiplied by the cosine score ('cosine')
            'MaCC_Score',
        ]
        with open(outFile, 'w', newline='') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerow(columns)
            if len(self.libraryTags) == 0:
                return None

            # Looping format is similar to functions such as reduce_final_df(). I'd consolidate them, but it was getting tricky to use numba.
            curLibTag = self.libraryTags[0]
            curQueTag = self.queryTags[0]
            curIonCount = self.queryIntensities[0]
            count = 1
            length = len(self.libraryTags)
            for i in range(1, length):
                if self.libraryTags[i] != curLibTag or self.queryTags[i] != curQueTag:
                    curScore = self.scores[i-1]
                    if curScore >= scoreCutoff:
                        libKey = idToKeyDict[curLibTag]
                        scan = str(curQueTag)
                        temp = [
                            expSpectraFile,  # fileName
                            scan,  # scan
                            queValDict[scan]['precursorMz'],  # MzEXP
                            libKey[1],  # peptide
                            lib[libKey]['ProteinName'],  # protein
                            libKey[0],  # MzLIB
                            lib[libKey]['PrecursorCharge'],  # zLIB
                            curScore/(count**(1/5)),  # cosine
                            lib[libKey]['transition_group_id'],  # name
                            queValDict[scan]['peaksCount'],  # Peaks(Query)
                            len(lib[libKey]['Peaks']),  # Peaks(Library)
                            count,  # shared
                            curIonCount,  # ionCount
                            queValDict[scan]['CV'],  # compensationVoltage
                            # totalWindowWidth
                            queValDict[scan]['windowWideness'],
                            curScore
                        ]
                        writer.writerow(temp)
                    curLibTag = self.libraryTags[i]
                    curQueTag = self.queryTags[i]
                    curIonCount = self.queryIntensities[i]
                    count = 1
                else:
                    curIonCount += self.queryIntensities[i]
                    count += 1
            curScore = self.scores[-1]
            if curScore >= scoreCutoff:
                libKey = idToKeyDict[self.libraryTags[-1]]
                scan = str(self.queryTags[-1])
                temp = [
                    expSpectraFile,  # fileName
                    scan,  # scan
                    queValDict[scan]['precursorMz'],  # MzEXP
                    libKey[1],  # peptide
                    lib[libKey]['ProteinName'],  # protein
                    libKey[0],  # MzLIB
                    lib[libKey]['PrecursorCharge'],  # zLIB
                    curScore/(count**(1/5)),  # cosine
                    lib[libKey]['transition_group_id'],  # name
                    queValDict[scan]['peaksCount'],  # Peaks(Query)
                    len(lib[libKey]['Peaks']),  # Peaks(Library)
                    count,  # shared
                    curIonCount,  # ionCount
                    queValDict[scan]['CV'],  # compensationVoltage
                    queValDict[scan]['windowWideness'],  # totalWindowWidth
                    curScore
                ]
                writer.writerow(temp)
