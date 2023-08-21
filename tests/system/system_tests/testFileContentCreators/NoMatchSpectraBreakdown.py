from .BaselineSpectraBreakdown import BaselineSpectraBreakdown
from .SpectraBreakdown import (
    create_vector_that_can_be_used_to_create_cosine_score,
    find_mz_window_of_query_scan_from_spectra,
)


class NoMatchSpectraBreakdown(BaselineSpectraBreakdown):
    def add_query_file_components_to_library_spectra_summary(
        self, librarySpectraBreakdown
    ):
        for mzWindowGroup in librarySpectraBreakdown:
            (
                mzWindowGroup["queryScanMz"],
                mzWindowGroup["windowWidth"],
            ) = find_mz_window_of_query_scan_from_spectra(mzWindowGroup["libSpectra"])
            topLibraryMz = (
                mzWindowGroup["queryScanMz"] + mzWindowGroup["windowWidth"] / 2
            )
            mzWindowGroup["querySpectra"] = [
                create_no_match_query_spectrum_from_library_spectrum_and_cosine_score(
                    spectrum, mzWindowGroup["cosine"], topLibraryMz
                )
                for spectrum in mzWindowGroup["libSpectra"]
            ]
        return librarySpectraBreakdown


def create_no_match_query_spectrum_from_library_spectrum_and_cosine_score(
    spectrum, cosineScore, topLibraryMz
):
    queryMz = spectrum["mzs"] + topLibraryMz + 1

    return {
        "mzs": queryMz,
        "intensities": create_vector_that_can_be_used_to_create_cosine_score(
            spectrum["intensities"], cosineScore
        ),
    }
