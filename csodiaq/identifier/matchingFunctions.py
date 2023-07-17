import pandas as pd

def calculate_parts_per_million_relative_difference(referenceMz, targetMz):
    return (referenceMz - targetMz) * (1e6) / referenceMz

def is_within_tolerance(ppm, tolerance):
    return abs(ppm) <= tolerance

def match_library_to_query_pooled_spectra(libraryPeaks, queryPeaks, ppmTolerance):
    """
    Peaks between library and query spectra are compared, and peaks with an m/z ppm tolerance
        are returned as matches.


    Parameters
    ----------
    libraryPeaks : list
        Each value in the list is a 3-length tuple with the following characteristics:
            m/z - the m/z value of the peak.
            intensity - the intensity value of the peak.
            tag - an int representing the library spectrum to which the peak belongs.
                This characteristics allows for pooling multiple library spectra
                together and being able to parse them out from each other after the
                analysis.

    queryPeaks : list
        Uses the same format as libraryPeaks, but the tag identifies query spectra
            instead.

    ppmTolerance : int
        The ppm tolerance allowed between library and query spectra peaks m/z values
            to be considered a match.

    Returns
    -------
    matchDf : pandas DataFrame
        A dataframe containing the library peak to query peak matches. Data includes
            the library/query tags of the peaks, their intensities, and the ppm
            difference between the two.
    """

    i, j = 0, 0
    mzIdx, intensityIdx, tagIdx = 0, 1, 2
    data = []
    while i < len(libraryPeaks) and j < len(queryPeaks):
        ppm = calculate_parts_per_million_relative_difference(libraryPeaks[i][mzIdx], queryPeaks[j][mzIdx])
        if not is_within_tolerance(ppm, ppmTolerance):
            if libraryPeaks[i][mzIdx] > queryPeaks[j][mzIdx]:
                j += 1
                continue
            if libraryPeaks[i][mzIdx] < queryPeaks[j][mzIdx]:
                i += 1
                continue
        p = i + 0
        while (p < len(libraryPeaks)):
            ppm = calculate_parts_per_million_relative_difference(libraryPeaks[p][mzIdx], queryPeaks[j][mzIdx])
            if not is_within_tolerance(ppm, ppmTolerance):
                break
            data.append(
                [libraryPeaks[p][tagIdx], libraryPeaks[p][intensityIdx], queryPeaks[j][tagIdx], queryPeaks[j][intensityIdx], ppm]
            )
            p += 1
        j += 1
    return pd.DataFrame(data, columns=["libraryIdx","libraryIntensity","queryIdx","queryIntensity","ppmDifference"])

def eliminate_low_count_matches(matches, minNumMatches = 3):
    return matches.groupby(["libraryIdx","queryIdx"]).filter(lambda x: x.shape[0] >= minNumMatches).reset_index(drop=True)

def eliminate_matches_below_fdr_cutoff(matches, groupsAboveCutoff):
    return matches.groupby(["libraryIdx","queryIdx"]).filter(lambda x: x.name in groupsAboveCutoff)


