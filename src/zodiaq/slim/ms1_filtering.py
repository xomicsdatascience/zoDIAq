from zodiaq.identification.matchingFunctions import (
    calculate_parts_per_million_relative_difference,
    is_within_tolerance,
)

def placehoder_filter_library_keys_by_ms1_mz_precursor_values(ms1Spectrum, libKeys):
    return libKeys

def filter_library_keys_by_ms1_mz_precursor_values(ms1Spectrum, libKeys):
    ms1MzValues = [x[1] for x in ms1Spectrum]
    return _filter_lib_keys_by_ms1_values_all(libKeys, ms1MzValues)

def _filter_lib_keys_by_ms1_values_all(libKeys, ms1MzValues):
    filteredKeys = [(mz, peptide) for mz, peptide in libKeys if _is_lib_mz_found_in_ms1(mz, np.array(ms1MzValues))]
    peptides = [key[1] for key in filteredKeys]
    target = 'VPTANVSVVDLTC(UniMod:4)R'
    return filteredKeys

def _is_lib_mz_found_in_ms1_float_tolerance(libMz, ms1MzValues, floatTolerance = 0.1): #NOTE
    tempList = ms1MzValues[np.abs(libMz - ms1MzValues) <= tolerance]
    return any([is_within_tolerance(calculate_parts_per_million_relative_difference(libMz,ms1Mz),20) for ms1Mz in tempList])