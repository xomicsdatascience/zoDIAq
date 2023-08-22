from .BaselineSpectraBreakdown import BaselineSpectraBreakdown
import numpy as np

class StressTestSpectraBreakdown(BaselineSpectraBreakdown):
    def make_library_spectra_breakdown_summary(self, targetDecoyLibrarySpectra):
        return [
            {
                "title": "highCosineTargets",
                "scan": 1,
                "cosine": 0.99,
                "libSpectra": targetDecoyLibrarySpectra["targets"],
            },
            {
                "title": "highCosineDecoys",
                "scan": 4,
                "cosine": 0.98,
                "libSpectra": targetDecoyLibrarySpectra["decoys"],
            }
        ]

    def add_query_file_components_to_library_spectra_summary(
        self, librarySpectraBreakdown
    ):
        for mzWindowGroup in librarySpectraBreakdown:
            (
                mzWindowGroup["queryScanMz"],
                mzWindowGroup["windowWidth"],
            ) = find_mz_window_of_query_scan_from_spectra(mzWindowGroup["libSpectra"])
            mzWindowGroup["querySpectra"] = generate_random_massive_query_spectra()
        return librarySpectraBreakdown

    def make_expected_output_df(self, mzWindowSpectraBreakdown):
        return pd.DataFrame(
            [],
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

def generate_random_massive_query_spectra():
    peakNum = 1000000000
    mzs = np.random.uniform(100.0, 2000.0, peakNum)
    intensities = np.random.uniform(1000.0, 10000.0, peakNum)
    return {
        "mzs": mzs,
        "intensities": intensities,
    }

