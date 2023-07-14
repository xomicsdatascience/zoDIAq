from csodiaq.loaders.library.LibraryLoaderContext import LibraryLoaderContext
from csodiaq.loaders.query.QueryLoaderContext import QueryLoaderContext
from csodiaq.identifier.poolingFunctions import generate_pooled_library_and_query_spectra_by_mz_windows
from csodiaq.identifier.matchingFunctions import match_library_to_query_pooled_spectra, eliminate_low_count_matches, eliminate_matches_below_fdr_cutoff
from csodiaq.identifier.scoringFunctions import score_library_to_query_matches, identify_all_decoys, determine_index_of_fdr_cutoff, filter_matches_by_ppm_offset_and_tolerance, calculate_ppm_offset_tolerance
from csodiaq.identifier.outputWritingFunctions import extract_metadata_from_match_and_score_dataframes, format_output_line, format_output_as_pandas_dataframe
import pandas as pd

class Identifier():
    def __init__(self, args):
        self.args = args
        self.libraryDict = LibraryLoaderContext(self.args["library"]).load_csodiaq_library_dict()
        self.decoySet = set([key for key,value in self.libraryDict.items() if value["isDecoy"]])

    def identify_library_spectra_in_queries(self):
        for queryFile in self.args["files"]:
            self.queryContext = QueryLoaderContext(queryFile)
            matchDf = self.match_library_to_query_spectra()
            scoreDf = self.score_spectra_matches(matchDf)
            matchDf, scoreDf = self.apply_correction_to_dataframes(matchDf, scoreDf)
            self.write_identifications_to_dataframe(matchDf, scoreDf)

    def match_library_to_query_spectra(self):
        matchDfs = []
        for pooledLibPeaks, pooledQueryPeaks in generate_pooled_library_and_query_spectra_by_mz_windows(
                self.libraryDict, self.queryContext):
            matchDf = match_library_to_query_pooled_spectra(pooledLibPeaks, pooledQueryPeaks,
                                                            self.args["fragmentMassTolerance"])
            matchDf = eliminate_low_count_matches(matchDf)
            matchDfs.append(matchDf)
        return pd.concat(matchDfs)

    def score_spectra_matches(self, matchDf):
        scoreDf = score_library_to_query_matches(matchDf)
        isDecoyArray = identify_all_decoys(self.decoySet, scoreDf)
        scoreDfCutoffIdx = determine_index_of_fdr_cutoff(isDecoyArray)
        return scoreDf.iloc[:scoreDfCutoffIdx, :]

    def apply_correction_to_dataframes(self, matchDf, scoreDf):
        if self.args["correction"] == -1:
            return matchDf, scoreDf
        aboveCutoffGroups = set(scoreDf.groupby(["libraryIdx", "queryIdx"]).groups)
        matchDf = eliminate_matches_below_fdr_cutoff(matchDf, aboveCutoffGroups)
        offset, tolerance = calculate_ppm_offset_tolerance(matchDf["ppmDifference"], self.args["correction"])
        matchDf = filter_matches_by_ppm_offset_and_tolerance(matchDf, offset, tolerance)
        matchDf = eliminate_low_count_matches(matchDf)
        scoreDf = score_library_to_query_matches(matchDf)
        isDecoyArray = identify_all_decoys(self.decoySet, scoreDf)
        scoreDfCutoffIdx = determine_index_of_fdr_cutoff(isDecoyArray)
        scoreDf = scoreDf.iloc[:scoreDfCutoffIdx, :]
        return matchDf, scoreDf

    def write_identifications_to_dataframe(self, matchDf, scoreDf):
        matchDict = extract_metadata_from_match_and_score_dataframes(matchDf, scoreDf)
        queryDict = self.queryContext.extract_metadata_from_query_scans()
        sortedLibKeys = sorted(self.libraryDict.keys())
        outputs = []
        for key, matchMetadata in matchDict.items():
            libKeyIdx, queryScan = key
            libKey = sortedLibKeys[libKeyIdx]
            libraryMetadata = self.libraryDict[libKey]
            libraryMetadata["precursorMz"] = libKey[0]
            libraryMetadata["peptide"] = libKey[1]
            queryMetadata = queryDict[str(queryScan)]
            queryMetadata["scan"] = queryScan
            outputLine = format_output_line(libraryMetadata, queryMetadata, matchMetadata)
            outputs.append(outputLine)
        outputDf = format_output_as_pandas_dataframe(self.queryContext.filePath, outputs)
        outFileHeader = self.args['outDirectory'] + 'CsoDIAq-file' + '_' + '.'.join(
            self.queryContext.filePath.split('/')[-1].split('.')[:-1])
        if self.args['correction'] != -1:
            outFileHeader += '_corrected_new2'
        outFile = outFileHeader + '.csv'
        outputDf.to_csv(outFile, index=False)