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

def add_decoy_label_to_score_df(decoySet, scoreDf):
    scoreDf["isDecoy"] = np.where(scoreDf["libraryIdx"].isin(decoySet), 1, 0)
    return scoreDf

def determine_index_of_fdr_cutoff(isDecoySeries):
    if isDecoySeries[0]:
        raise ValueError('None of the library peptides were identified in the query spectra (highest score was a decoy).')
    fdrCutoff = 1e-2
    decoyIndices = np.flatnonzero(isDecoySeries)
    fdrs = (np.arange(1, len(decoyIndices)+1))/(decoyIndices + 1)
    lastDecoyIdx = np.argmax(fdrs > fdrCutoff)
    return decoyIndices[lastDecoyIdx] - 1

