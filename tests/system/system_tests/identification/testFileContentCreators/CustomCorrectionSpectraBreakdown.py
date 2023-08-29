from .BaselineSpectraBreakdown import (
    BaselineSpectraBreakdown,
    calculate_mz_value_from_ppm,
)
import numpy as np


class CustomCorrectionSpectraBreakdown(BaselineSpectraBreakdown):
    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        numBins = 200
        ppmCenterBin = 10.0
        binWidth = 0.1
        emptyBinNum = 4
        centerBinDistribution = np.random.uniform(
            ppmCenterBin - (binWidth / 2), ppmCenterBin + (binWidth / 2), size=(4,)
        )
        ppmLeftBin = ppmCenterBin - binWidth
        leftBinDistribution = np.random.uniform(
            ppmLeftBin - (binWidth / 2), ppmLeftBin + (binWidth / 2), size=(2,)
        )
        ppmRightBin = ppmCenterBin + binWidth
        rightBinDistribution = np.random.uniform(
            ppmRightBin - (binWidth / 2), ppmRightBin + (binWidth / 2), size=(2,)
        )
        farLeftBin = 10 - ((numBins - 3) // 2) * binWidth
        throwawayLeftDistribution = np.random.uniform(
            farLeftBin - (binWidth / 2),
            ppmLeftBin - emptyBinNum * binWidth - (binWidth / 2),
            size=(1,),
        )
        farRightBin = 10 + ((numBins - 3) // 2 + 1) * binWidth
        throwawayRightDistribution = np.random.uniform(
            ppmRightBin + emptyBinNum * binWidth + (binWidth / 2),
            farRightBin + (binWidth / 2),
            size=(1,),
        )
        ppmDistribution = np.concatenate(
            [
                throwawayLeftDistribution,
                leftBinDistribution,
                centerBinDistribution,
                rightBinDistribution,
                throwawayRightDistribution,
            ]
        )
        return [
            calculate_mz_value_from_ppm(libraryMzs[i], ppmDistribution[i])
            for i in range(len(libraryMzs))
        ]
