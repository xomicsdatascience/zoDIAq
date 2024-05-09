from zodiaq.identification.matchingFunctions import (
    calculate_parts_per_million_relative_difference,
    is_within_tolerance,
)
import numpy as np

def placehoder_filter_library_keys_by_ms1_mz_precursor_values(ms1Spectrum, libKeys):
    return libKeys

def filter_library_keys_by_ms1_mz_precursor_values(ms1Spectrum, libKeys):
    ms1MzValues = [x[0] for x in ms1Spectrum]
    return _filter_lib_keys_by_ms1_values(libKeys, ms1MzValues)

def _filter_lib_keys_by_ms1_values(libKeys, ms1MzValues):
    filteredKeys = [(mz, peptide) for mz, peptide in libKeys if _is_lib_mz_found_in_ppm_tolerance_of_an_ms1_peak(mz, np.array(ms1MzValues))]
    return filteredKeys

def _is_lib_mz_found_in_ppm_tolerance_of_an_ms1_peak(libMz, ms1MzValues, ppmTolerance=20): #NOTE
    tempList = ms1MzValues[np.abs(libMz - ms1MzValues) <= 0.1]
    return any([is_within_tolerance(calculate_parts_per_million_relative_difference(libMz,ms1Mz),ppmTolerance) for ms1Mz in tempList])