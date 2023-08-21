from .BaselineSpectraBreakdown import BaselineSpectraBreakdown
from .SpectraBreakdown import (
    create_vector_that_can_be_used_to_create_cosine_score,
    find_mz_window_of_query_scan_from_spectra,
)
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
