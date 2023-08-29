from .BaselineSpectraBreakdown import (
    BaselineSpectraBreakdown,
    calculate_mz_value_from_ppm,
)
import numpy as np


class StDevCorrectionSpectraBreakdown(BaselineSpectraBreakdown):
    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        ppmCenter = 10.0
        stDev = 1.0
        ppmDistribution = np.random.normal(ppmCenter, stDev, len(libraryMzs))
        return [
            calculate_mz_value_from_ppm(libraryMzs[i], ppmDistribution[i])
            for i in range(len(libraryMzs))
        ]
