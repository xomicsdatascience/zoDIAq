from abc import ABC, abstractmethod
import os
from zodiaq.loaders.library.decoyGenerationFunctions import (
    determine_if_decoys_should_be_generated,
    add_decoys_to_zodiaq_library,
)
import random


class LibraryLoaderStrategy(ABC):
    def __init__(self):
        self.rawUploadedLibraryObject = None

    @abstractmethod
    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        """
        Loads the library file in a format appropriate for the concrete strategy implementation.
            The object is saved as an attribute (rawLib__) for use by other class functions.
            This is an abstract method and must be implemented in child classes.

        Parameters
        ----------
        libraryFilePath : string (os.PathLike format)
            Path to the library file.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the library file don't match any scoped library format.
        """
        pass

    @abstractmethod
    def _format_raw_library_object_into_zodiaq_library_dict(self) -> dict:
        """
        Creates a standardized dictionary object comprised of library spectra.
            Relies on a raw library object created as an attribute 'rawLib__' in the
                _load_raw_library_object_from_file function.
            This is an abstract method and must be implemented in child classes.

        Parameters
        ----------
        None.

        Returns
        -------
        zodiaqLibDict : dict
            key: (precursorMz, peptideName) tuple
            value: dictionary
                contains the following `key - value` entries:
                    precursorCharge - int of the precursor charge.
                    identification - string of the unique title of the library spectrum.
                    proteinName - string of protein(s) this peptide is found in.
                    peaks - list of tuples repres. peaks in the spectrum. Tuples are:
                        (peak mz value, peak intensity value, zodiaqLibDict key idx)
                    zodiaqKeyIdx - int representing the zodiaqLibDict key idx,
                        matching the last element of peak tuples.
                    isDecoy - int, indicating if the spectrum represents a decoy
                        0 - is not a decoy
                        1 - is a decoy
        """
        pass

    def load_zodiaq_library_dict_from_file(
        self, libraryFilePath: os.PathLike, isTest
    ) -> dict:
        """
        Outputs a standardized dictionary object from a library file.
            The LibraryLoaderContext class calls this function, which uses
            the concrete-strategy-implemented abstract methods of the appropriate
            library format.

        Parameters
        ----------
        libraryFilePath : string (os.PathLike format)
            Path to the library file.

        Returns
        -------
        zodiaqLibDict : dict
            see _format_raw_library_object_into_zodiaq_library_dict return value.
        """
        random.seed(0)
        self._load_raw_library_object_from_file(libraryFilePath)
        zodiaqLibDict = self._format_raw_library_object_into_zodiaq_library_dict()
        format_proteins_into_list_format(zodiaqLibDict)
        if determine_if_decoys_should_be_generated(zodiaqLibDict) and not isTest:
            zodiaqLibDict = add_decoys_to_zodiaq_library(zodiaqLibDict)
        return zodiaqLibDict


def format_proteins_into_list_format(zodiaqLibDict):
    for key, value in zodiaqLibDict.items():
        if value["proteinName"] and not value["proteinName"][:1].isdigit():
            value["proteinName"] = f'1/{value["proteinName"]}'


def create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
    mz: list, intensities: list, id: int
) -> list:
    """
    Combines matching mz and intensity lists into a list of tuples.
        Includes a spectrum-identififying index as a third value of each tuple.
        (peak mz value, peak intensity value, zodiaqLibDict key idx)

    Parameters
    ----------
    mz : list
        List of floats corresponding to peak mz values.

    intensities : list
        List of floats corresponding to peak intensity values.

    id : int
        A spectrum-identififying index value.

    Returns
    -------
    peaks : list
        see _format_raw_library_object_into_zodiaq_library_dict return value for 'peak'.
    """
    idList = [id] * len(mz)
    return list(zip(mz, intensities, idList))


def remove_low_intensity_peaks_below_max_peak_num(peaks: list, maxPeakNum: int) -> list:
    """
    Reduces a pre-existing list of peaks to a maximum length,
        prioritizing peak intensity values.

    Parameters
    ----------
    peaks : list
        see _format_raw_library_object_into_zodiaq_library_dict return value for 'peak'.

    maxPeakNum : int
        Maximum number of peaks to be returned.

    Returns
    -------
    peaks : list
        see _format_raw_library_object_into_zodiaq_library_dict return value for 'peak',
        but reduced to a smaller length if the original length was greater than the max length.
    """
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks[:maxPeakNum]


finalVariableNames = {
    "precursorCharge": "precursorCharge",
    "identification": "identification",
    "proteinName": "proteinName",
    "peaks": "peaks",
    "zodiaqKeyIdx": "zodiaqKeyIdx",
    "isDecoy": "isDecoy",
    "fragmentTypes": "fragmentTypes",
}
