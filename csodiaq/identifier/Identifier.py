from csodiaq.loaders.library.LibraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.query.QueryLoaderContext import QueryLoaderContext
from csodiaq.identifier.poolingFunctions import generate_pooled_library_and_query_spectra_by_mz_windows
from csodiaq.identifier.matchingFunctions import match_library_to_query_pooled_spectra, eliminate_low_count_matches, eliminate_matches_below_fdr_cutoff
from csodiaq.identifier.fdrFunctions import score_library_to_query_matches, identify_all_decoys, determine_index_of_fdr_cutoff, filter_matches_by_ppm_offset_and_tolerance
class Identifier():
    def __init__(self, args):
        self.args = args
        self.libraryDict = LibraryLoaderContext(self.args["library"]).load_csodiaq_library_dict()
        self.decoySet = set([key for key,value in self.libraryDict.items() if value["isDecoy"]])

    def identify_library_spectra_in_queries(self):
        outputs = []
        for queryFile in self.args["files"]:
            self.queryContext = QueryLoaderContext(queryFile)
            for pooledLibPeaks, pooledQueryPeaks in generate_pooled_library_and_query_spectra_by_mz_windows(libDict, queryContext):
                matchDf = match_library_to_query_pooled_spectra(pooledLibPeaks, pooledQueryPeaks, ppmTolerance)
                matchDf = eliminate_low_count_matches(matchDf)
                scoreDf = score_library_to_query_matches(matchDf)
                isDecoyArray = identify_all_decoys(self.decoySet, scoreDf)
                scoreDfCutoffIdx = determine_index_of_fdr_cutoff(isDecoyArray)
                scoreDf = scoreDf.iloc[:scoreDfCutoffIdx,:]
                if self.args["corrected"] != -1:
                    aboveCutoffGroups = set(scoreDf.groupby("libraryIndex", "queryIndex").groups())
                    matchDf = eliminate_matches_below_fdr_cutoff(matchDf, aboveCutoffGroups)
                    offset, tolerance = calculate_ppm_offset_tolerance(matchDf["ppmDifference"], args["corrected"])
                    matchDf = filter_matches_by_ppm_offset_and_tolerance(matchDf, offset, tolerance)
                    matchDf = eliminate_low_count_matches(matchDf)
                    scoreDf = score_library_to_query_matches(matchDf)
                    isDecoyArray = identify_all_decoys(self.decoySet, scoreDf)
                    scoreDfCutoffIdx = determine_index_of_fdr_cutoff(isDecoyArray)
                    scoreDf = scoreDf.iloc[:scoreDfCutoffIdx, :]

