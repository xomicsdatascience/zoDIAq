import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def calculate_parts_per_million_relative_difference(referenceMz, targetMz):
    return (referenceMz - targetMz) * (1e6) / referenceMz

def is_within_tolerance(ppm, tolerance):
    return abs(ppm) <= tolerance

def match_library_to_query_pooled_spectra(libraryPeaks, queryPeaks, ppmTolerance):
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

def eliminate_low_count_matches(matches):
    minNumMatches = 3
    return matches.groupby(["libraryIdx","queryIdx"]).filter(lambda x: x.shape[0] >= minNumMatches).reset_index(drop=True)

def eliminate_matches_below_fdr_cutoff(matches, groupsAboveCutoff):
    return matches.groupby(["libraryIdx","queryIdx"]).filter(lambda x: x.name in groupsAboveCutoff)


