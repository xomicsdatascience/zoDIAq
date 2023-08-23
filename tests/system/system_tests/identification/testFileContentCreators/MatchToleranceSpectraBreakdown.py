from .BaselineSpectraBreakdown import (
    BaselineSpectraBreakdown,
    calculate_mz_value_from_ppm,
)


class MatchToleranceSpectraBreakdown(BaselineSpectraBreakdown):
    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        firstPpm = 29.5
        first = libraryMzs[0]
        libraryMzs[0] = calculate_mz_value_from_ppm(libraryMzs[0], firstPpm)
        secondPpm = 28.5
        libraryMzs[1] = calculate_mz_value_from_ppm(libraryMzs[1], secondPpm)
        return libraryMzs
