from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy, create_peaks_from_mz_intensity_lists_and_csodiaq_key_id, remove_low_intensity_peaks_below_max_peak_num
import pandas as pd
import os
from pyteomics import mgf
import re

class LibraryLoaderStrategyMgf(LibraryLoaderStrategy):

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawLibMgf = mgf.read(libraryFilePath)
        # NOTE: check if required values exist??

    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict:
        maxPeakNum = 10
        csodiaqLibDict = {}
        for csodiaqKeyIdx in range(len(self.rawLibMgf)):
            librarySpectrum = self.rawLibMgf[csodiaqKeyIdx]
            csodiaqKey, csodiaqValue = create_csodiaq_library_entry(librarySpectrum, maxPeakNum, csodiaqKeyIdx)
            csodiaqLibDict[csodiaqKey] = csodiaqValue
        return csodiaqLibDict

def create_csodiaq_library_entry(librarySpectrum: dict, maxPeakNum: int, csodiaqKeyIdx: int) -> dict:
    csodiaqKey, precursorCharge, identifier, proteinName, isDecoy = extract_variables_from_spectrum_metadata(librarySpectrum['params'])
    peaks = create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(librarySpectrum['m/z array'],
                                                                    librarySpectrum['intensity array'],
                                                                    csodiaqKeyIdx)
    reducedPeaks = remove_low_intensity_peaks_below_max_peak_num(peaks, maxPeakNum)
    return csodiaqKey, {
        'precursorCharge': precursorCharge,
        'identifier': identifier,
        'proteinName': proteinName,
        'peaks': sorted(reducedPeaks),
        'csodiaqKeyIdx': csodiaqKeyIdx,
        'isDecoy': isDecoy,
    }

def extract_variables_from_spectrum_metadata(metadata):
    csodiaqKey = (metadata['pepmass'][0], metadata['seq'])
    charge = int(re.sub('[+-]', '', str(metadata['charge'][0])))
    name = metadata['title']
    if 'protein' in metadata:
        protein = metadata['protein']
    else:
        protein = ''
    if 'DECOY' in name:
        decoy = 1
    else:
        decoy = 0
    return csodiaqKey, charge, name, protein, decoy