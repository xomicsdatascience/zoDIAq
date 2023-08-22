from .BaselineSpectraBreakdown import BaselineSpectraBreakdown

class MatchToleranceSpectraBreakdown(BaselineSpectraBreakdown):
    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        firstPpm = 29.5
        first = libraryMzs[0]
        libraryMzs[0] = calculate_mz_value_from_ppm(libraryMzs[0], firstPpm)
        secondPpm = 28.5
        libraryMzs[1] = calculate_mz_value_from_ppm(libraryMzs[1], secondPpm)
        return libraryMzs

def calculate_mz_value_from_ppm(mz, ppm):
    return mz * (1e6-ppm) / 1e6

def calculate_parts_per_million_relative_difference(referenceMz, targetMz):
    return (referenceMz - targetMz) * (1e6) / referenceMz