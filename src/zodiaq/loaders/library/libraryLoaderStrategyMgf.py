from zodiaq.loaders.library.libraryLoaderStrategy import (
    LibraryLoaderStrategy,
    create_peaks_from_mz_intensity_lists_and_zodiaq_key_id,
    remove_low_intensity_peaks_below_max_peak_num,
    finalVariableNames,
)
import os
from pyteomics import mgf
import re


class LibraryLoaderStrategyMgf(LibraryLoaderStrategy):
    """
    Concrete strategy implementation of the LibraryLoaderStrategy strategy class specific for
        loading mgf library files into a standardized dictionary for use in zoDIAq.

    Attributes
    ----------
    rawUploadedLibraryObject : pyteomics.mgf.IndexedMGF
        The original mgf library as loaded using the pyteomics mgf class.
    """

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawUploadedLibraryObject = mgf.read(libraryFilePath)

    def _format_raw_library_object_into_zodiaq_library_dict(self) -> dict:
        maxPeakNum = 10
        zodiaqLibDict = {}
        for zodiaqKeyIdx in range(len(self.rawUploadedLibraryObject)):
            librarySpectrum = self.rawUploadedLibraryObject[zodiaqKeyIdx]
            zodiaqKey, zodiaqValue = _create_zodiaq_library_entry(
                librarySpectrum, maxPeakNum, zodiaqKeyIdx
            )
            zodiaqLibDict[zodiaqKey] = zodiaqValue
        zodiaqLibDict = _rewrite_zodiaq_key_idx_to_correspond_to_sorted_zodiaq_key(
            zodiaqLibDict
        )
        return zodiaqLibDict


def _rewrite_zodiaq_key_idx_to_correspond_to_sorted_zodiaq_key(zodiaqLibDict):
    sortedLibKeys = sorted(zodiaqLibDict.keys())
    for i, key in enumerate(sortedLibKeys):
        oldPeaks = zodiaqLibDict[key][finalVariableNames["peaks"]]
        newPeaks = [(peak[0], peak[1], i) for peak in oldPeaks]
        zodiaqLibDict[key][finalVariableNames["peaks"]] = newPeaks
    return zodiaqLibDict


def _create_zodiaq_library_entry(
    librarySpectrum: dict, maxPeakNum: int, zodiaqKeyIdx: int
) -> dict:
    (
        zodiaqKey,
        precursorCharge,
        identifier,
        proteinName,
        isDecoy,
    ) = _extract_variables_from_spectrum_metadata(librarySpectrum["params"])
    peaks = create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
        librarySpectrum["m/z array"], librarySpectrum["intensity array"], zodiaqKeyIdx
    )
    reducedPeaks = remove_low_intensity_peaks_below_max_peak_num(peaks, maxPeakNum)
    return zodiaqKey, {
        finalVariableNames["precursorCharge"]: precursorCharge,
        finalVariableNames["identification"]: identifier,
        finalVariableNames["proteinName"]: proteinName,
        finalVariableNames["peaks"]: sorted(reducedPeaks),
        finalVariableNames["zodiaqKeyIdx"]: zodiaqKeyIdx,
        finalVariableNames["isDecoy"]: isDecoy,
    }


def _extract_variables_from_spectrum_metadata(metadata):
    zodiaqKey = (metadata["pepmass"][0], metadata["seq"])
    charge = int(re.sub("[+-]", "", str(metadata["charge"][0])))
    name = metadata["title"]
    if "protein" in metadata:
        protein = metadata["protein"]
    else:
        protein = ""
    if "DECOY" in name:
        decoy = 1
    else:
        decoy = 0
    return zodiaqKey, charge, name, protein, decoy
