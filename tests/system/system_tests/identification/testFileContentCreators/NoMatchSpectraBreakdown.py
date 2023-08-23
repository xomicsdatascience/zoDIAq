from .BaselineSpectraBreakdown import BaselineSpectraBreakdown


class NoMatchSpectraBreakdown(BaselineSpectraBreakdown):
    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        topLibraryMz = mzWindowGroup["queryScanMz"] + mzWindowGroup["windowWidth"] / 2
        return libraryMzs + topLibraryMz + 1
