import numpy as np
from scipy.spatial.distance import cosine

def score_library_to_query_matches(matches):
    output = matches\
                .groupby(["libraryIdx", "queryIdx"])\
                .apply(lambda x: calculate_macc_score(x["libraryIntensity"], x["queryIntensity"]))\
                .reset_index(name='score')
    return output

def calculate_cosine_similarity_score(vectorA, vectorB):
    return 1 - cosine(vectorA, vectorB)

def calculate_macc_score(vectorA, vectorB):
    return len(vectorA.index)**(1/5) * calculate_cosine_similarity_score(vectorA, vectorB)

def identify_all_decoys(decoySet, scoreDf):
    return np.where(scoreDf["libraryIdx"].isin(decoySet), 1, 0)

def determine_index_of_fdr_cutoff(isDecoyArray):
    if isDecoyArray[0]:
        raise ValueError('None of the library peptides were identified in the query spectra (highest score was a decoy).')
    fdrCutoff = 1e-2
    decoyIndices = np.flatnonzero(isDecoyArray)
    fdrs = (np.arange(1, len(decoyIndices)+1))/(decoyIndices + 1)
    lastDecoy = np.argmax(fdrs > fdrCutoff)
    return decoyIndices[lastDecoy]

