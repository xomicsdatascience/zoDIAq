import csodiaq_base_functions as cbf
import numpy as np
from numba import njit

@njit
def find_ppm_offset(mz, ppm):
    return (mz*-ppm)/1000000


@njit
def is_approx(x, y, ppmTol):
    if x==y: return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)

@njit
def cosine_similarity(AB, A, B):
    magnitude = (A**0.5) * (B**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
    return (AB / magnitude if magnitude else 0)

class PooledSpectraMatcher:
    def __init__(self):
        self.libraryTags = np.array([],dtype=int)
        self.queryTags = np.array([],dtype=int)
        self.libraryIntensities = np.array([],dtype=float)
        self.queryIntensities = np.array([],dtype=float)
        self.ppmMatches = np.array([],dtype=float)
        self.decoys = np.array([],dtype=float)
        self.scores = np.array([],dtype=float)

    def compare_spectra(self, librarySpectrum, querySpectrum, tolerance, offset):

            libMzs, libIntensities, libTags = list(map(list, zip(*librarySpectrum)))
            queMzs, queIntensities, queTags = list(map(list, zip(*querySpectrum)))
            self.find_matching_peaks(libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, offset, tolerance)
            self.remove_matches_with_too_few_peaks()


    @njit
    def find_matching_peaks(self, libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, offset, tolerance):

        lenLib = len(libMzs)
        lenQue = len(queMzs)

        # By tracking the indices of the current library/query peaks we reduce the time complexity of the algorithm
        i, j = 0, 0
        expPeakMz = queMzs[j] - find_ppm_offset(queMzs[j], ppmOffset)
        while i < lenLib and j < lenQue:

            # If the m/z of the peaks are not within the given ppm tolerance, the indices of the smaller of the two is incremented
            #   and we are returned to the top of the while loop.
            if not is_approx(libMzs[i],expPeakMz, ppmTol):
                if libMzs[i] > expPeakMz:
                    j += 1
                    if j < lenQue: expPeakMz = queMzs[j] - find_ppm_offset(queMzs[j], ppmOffset)
                    continue
                if libMzs[i] < expPeakMz: i += 1; continue

            # To account for the matching of one query peak to multiple library peaks, library peaks are looped over
            #   after the initial match. Every matching peak contributes to the various variables in the final returned dictionary.
            p = i + 0
            while (p < lenLib):
                ppm = is_approx(libMzs[p], expPeakMz, ppmTol)
                if p==lenLib or not ppm: break
                self.libraryTags.append(libTags[p])
                self.libraryIntensities.append(libIntensities[p])
                self.queryTags.append(queTags[j])
                self.queryIntensities.append(queIntensities[j])
                self.ppmMatches.append(ppm)

                p += 1

            # Note that the possibility of one library peak matching to multiple query peaks is automatically accounted for
            #   by the fact that the query peak is the next default increment after all match calculations have been made.
            j += 1
            if j < lenQue: expPeakMz = queMzs[j] - find_ppm_offset(queMzs[j], ppmOffset)

    def remove_matches_with_too_few_peaks_and_generate_scores(self):
        self.sort_matches_by_tags()
        remove = self.find_indices_to_remove_and_generate_scores()
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [np.delete(x,remove) for x in [self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def sort_matches_by_tags(self):
        i1 = self.queryTags.argsort(kind='mergesort')
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i1] for x in [self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]
        i2 = self.libraryTags.argsort(kind='mergesort')
        self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches = [x[i2] for x in [self.libraryTags, self.libraryIntensities, self.queryTags, self.queryIntensities, self.ppmMatches]]

    def find_indices_to_remove_and_generate_scores(self):

        # Because the tags are sorted together, each spectrum-spectrum tag pair is grouped together. As soon as the pair changes, that means you've gone to to a new pair.
        #  We keep track of the current tags to determine if they have changed.
        curLibTag = self.libraryTags[0]
        curQueTag = self.queryTags[0]

        # count keeps track of the number of peak matches in each spectra match
        count = 1

        # keeping track of variables necessary to calculate cosine score and, consequentally, the macc score
        AB = self.libraryIntensities[0]*self.queryIntensities[0]
        A = self.libraryIntensities[0]**2
        B = self.queryIntensities[0]**2

        # initializing the return values
        remove = []

        # looping through every peak match
        length = len(self.libraryTags)
        for i in range(1,length):
            if self.libraryTags[i] != curLibTag or self.queryTags[i] != curQueTag:

                # if the count is 3 or higher, calculate the macc score and add it
                if count > 2:
                    test = 0
                    cosine = cosine_similarity(AB, A, B)
                    macc = (count**(1/5))*cosine
                    magnitude = (A**0.5) * (B**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
                    self.scores.extend([macc]*count)

                # if the count is 1 or 2, remove the peak match(es) represented by the indices
                else: remove.extend([i-j for j in range(1,count+1)])

                # reset the variables
                count = 1
                curLibTag = self.libraryTags[i]
                curQueTag = self.queryTags[i]
                AB = self.libraryIntensities[i]*self.queryIntensities[i]
                A = self.libraryIntensities[i]**2
                B = self.queryIntensities[i]**2

            # if its the same spectra tag pair, contribute to the variables used to calculate the cosine score
            else:
                AB += self.libraryIntensities[i]*self.queryIntensities[i]
                A += self.libraryIntensities[i]**2
                B += self.queryIntensities[i]**2
                count += 1

        # the middle of the loop is repeated for the last spectra key pairing as well
        if count > 2:
            test = 0
            cosine = cosine_similarity(AB, A, B)
            macc = (count**(1/5))*cosine
            self.scores.extend([macc]*count)
        else: remove.extend([length-j for j in range(1,count+1)])
    return remove

    def mark_libraries_as_targets_vs_decoys(self, decoys):
        self.decoys = np.append(self.decoys, decoys)

    def extend_by_other_spectra_match(self, spectraMatch): pass
