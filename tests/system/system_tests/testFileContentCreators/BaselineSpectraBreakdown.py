from .SpectraBreakdown import (
    SpectraBreakdown,
    create_vector_that_can_be_used_to_create_cosine_score,
    find_mz_window_of_query_scan_from_spectra,
    find_number_of_total_query_peaks_in_scan,
)
import numpy as np
import pandas as pd


class BaselineSpectraBreakdown(SpectraBreakdown):
    def make_library_spectra_breakdown_summary(self, targetDecoyLibrarySpectra):
        return [
            {
                "title": "highCosineTargets",
                "scan": 1,
                "cosine": 0.99,
                "libSpectra": targetDecoyLibrarySpectra["targets"][:200],
            },
            {
                "title": "highCosineDecoys",
                "scan": 2,
                "cosine": 0.98,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][:2],
            },
            {
                "title": "midCosineTargets",
                "scan": 3,
                "cosine": 0.97,
                "libSpectra": targetDecoyLibrarySpectra["targets"][200:400],
            },
            {
                "title": "midCosineDecoys",
                "scan": 4,
                "cosine": 0.96,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][2:5],
            },
            {
                "title": "lowCosineTargets",
                "scan": 5,
                "cosine": 0.95,
                "libSpectra": targetDecoyLibrarySpectra["targets"][400:],
            },
            {
                "title": "lowCosineDecoys",
                "scan": 6,
                "cosine": 0.94,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][5:],
            },
        ]

    def add_query_file_components_to_library_spectra_summary(
        self, librarySpectraBreakdown
    ):
        for mzWindowGroup in librarySpectraBreakdown:
            (
                mzWindowGroup["queryScanMz"],
                mzWindowGroup["windowWidth"],
            ) = find_mz_window_of_query_scan_from_spectra(mzWindowGroup["libSpectra"])
            mzWindowGroup["querySpectra"] = [
                create_query_spectrum_from_library_spectrum_and_cosine_score(
                    spectrum, mzWindowGroup["cosine"]
                )
                for spectrum in mzWindowGroup["libSpectra"]
            ]
        return librarySpectraBreakdown

    def make_expected_output_df(self, mzWindowSpectraBreakdown):
        data = []
        for cosineGroup in mzWindowSpectraBreakdown:
            cosine = cosineGroup["cosine"]
            numQueryPeaks = find_number_of_total_query_peaks_in_scan(
                cosineGroup["querySpectra"]
            )
            for i in range(len(cosineGroup["libSpectra"])):
                libSpectra = cosineGroup["libSpectra"][i]
                querySpectra = cosineGroup["querySpectra"][i]
                queryMzLargerThanPrecursorIdx = np.argwhere(
                    np.array(querySpectra["mzs"]) > cosineGroup["queryScanMz"]
                )
                queryMzLargerThanPrecursorIntensities = np.array(
                    querySpectra["intensities"]
                )[queryMzLargerThanPrecursorIdx]
                ionCount = 0
                if len(queryMzLargerThanPrecursorIntensities):
                    ionCount = sum(queryMzLargerThanPrecursorIntensities)[0]
                if "DECOY" in libSpectra["protein"]:
                    decoy = 1
                else:
                    decoy = 0
                data.append(
                    [
                        cosineGroup["scan"],
                        cosineGroup["queryScanMz"],
                        libSpectra["peptide"],
                        libSpectra["protein"],
                        decoy,
                        libSpectra["precursorMz"],
                        libSpectra["precursorCharge"],
                        cosine,
                        libSpectra["id"],
                        numQueryPeaks,
                        len(libSpectra["mzs"]),
                        len(libSpectra["mzs"]),
                        ionCount,
                        -cosineGroup["scan"],
                        cosineGroup["windowWidth"],
                        len(libSpectra["mzs"])
                        - len(queryMzLargerThanPrecursorIntensities),
                    ]
                )
        return pd.DataFrame(
            data,
            columns=[
                "scan",
                "MzEXP",
                "peptide",
                "protein",
                "isDecoy",
                "MzLIB",
                "zLIB",
                "cosine",
                "name",
                "Peak(Query)",
                "Peaks(Library)",
                "shared",
                "ionCount",
                "CompensationVoltage",
                "totalWindowWidth",
                "exclude_num",
            ],
        )


def create_query_spectrum_from_library_spectrum_and_cosine_score(spectrum, cosineScore):
    return {
        "mzs": spectrum["mzs"],
        "intensities": create_vector_that_can_be_used_to_create_cosine_score(
            spectrum["intensities"], cosineScore
        ),
    }
