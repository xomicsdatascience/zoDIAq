from .BaselineSpectraBreakdown import BaselineSpectraBreakdown
import numpy as np


class NoCompensationVoltageSpectraBreakdown(BaselineSpectraBreakdown):
    def make_expected_output_df(self, mzWindowSpectraBreakdown):
        outputDf = super().make_expected_output_df(mzWindowSpectraBreakdown)
        outputDf["CompensationVoltage"] = [np.nan] * len(outputDf)
        return outputDf

    def make_input_query_spectra(self, mzWindowSpectraBreakdown):
        spectra = super().make_input_query_spectra(mzWindowSpectraBreakdown)
        for spectrum in spectra:
            spectrum.compensationVoltage = None
        return spectra
