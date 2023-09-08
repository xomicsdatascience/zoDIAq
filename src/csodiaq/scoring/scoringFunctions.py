import numpy as np
from scipy.spatial.distance import cosine
import matplotlib.pyplot as pyplot


def score_library_to_query_matches(matches):
    scoreDf = (
        matches.groupby(["libraryIdx", "queryIdx"])
        .apply(
            lambda x: calculate_cosine_similarity_score(
                x["libraryIntensity"], x["queryIntensity"]
            )
        )
        .reset_index(name="cosineScore")
    )
    return scoreDf.sort_values("cosineScore", ascending=False).reset_index(drop=True)


def calculate_cosine_similarity_score(vectorA, vectorB):
    return 1 - cosine(vectorA, vectorB)


def calculate_macc_score(numMatchedPeaks, cosineScore):
    return numMatchedPeaks ** (1 / 5) * cosineScore


def determine_index_of_fdr_cutoff(isDecoyArray, fdrCutoff=1e-2):
    if isDecoyArray[0]:
        raise ValueError(
            "None of the library peptides were identified in the query spectra (highest score was a decoy)."
        )
    if 1 not in isDecoyArray:
        return len(isDecoyArray)
    fdrs = calculate_fdr_rates_of_decoy_array(isDecoyArray)
    lastDecoyIdx = np.argmax(fdrs > fdrCutoff)
    return lastDecoyIdx


def calculate_fdr_rates_of_decoy_array(isDecoyArray):
    isDecoyArrayIdx = np.arange(1, len(isDecoyArray) + 1)
    decoyCumSum = np.array(isDecoyArray).cumsum()
    return decoyCumSum / isDecoyArrayIdx
