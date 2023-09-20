from abc import ABC, abstractmethod
import os


class QueryLoaderStrategy(ABC):
    """
    Strategy implementation for loading query files.

    Attributes
    ----------
    filePath : string (os.PathLike format)
        The file path to the query data file.
    """

    def __init__(self, queryFilePath: os.PathLike):
        """
        Initialization of this strategy saves the path to the query file for use in the abstract methods.

        Parameters
        ----------
        queryFilePath : string (os.PathLike format)
            Path to the query file.
        """
        self.filePath = queryFilePath

    @abstractmethod
    def map_query_scan_ids_to_dia_mz_windows(self) -> dict:
        """
        Identifies m/z windows from each scan, then groups scans by identical m/z windows.


        Extended Summary
        ---------------
        M/z windows are a feature of mass spectrometry DIA methods.
            An m/z window is represented by a precursor m/z and total area around m/z.
            The upper/lower bound of the window is found by adding/subtracting HALF
            of the window size to the precursor m/z, respectively.

        Parameters
        ----------
        None

        Returns
        -------
        mzWindowToScanIdDict : dict
            key: (precursorMz, windowSize) tuple

            value: list
            List contains strings representing the scan number. Pyteomics allows for access
                to query spectrum directly from the file given the scan number of the spectrum.
        """
        pass

    @abstractmethod
    def extract_metadata_from_query_scans(self) -> dict:
        """
        Identifies metadata for each query spectrum/scan for reference in the final csodiaq output.

        Parameters
        ----------
        None

        Returns
        -------
        scanMetadataDict : dict
            key: string
            integer in string format representing the scan (spectrum) number.
            value: dictionary
                contains the following `key - value` entries:
                    precursorMz - float of the precursor m/z value.
                    windowWidth - float of the m/z window.
                    peaksCount - int of the number of peaks in the scan.
                    compensationVoltage - compensation voltage of the scan. Compensation voltage (CV)
                        is a characteristic of DISPA data. If data is not DISPA data, CV is rendered as
                        an empty string in the final metadata dictionary.
        """
        pass

    @abstractmethod
    def get_query_file_reader(self):
        pass


    @abstractmethod
    def pool_peaks_of_query_scans(self, scans: list, reader) -> list:
        """
        Given a list of scans, this function gathers their peaks and pools them into
            a single list.

        Parameters
        ----------
        scans : list
            List of scan identifiers. Scan identifiers are in string format for compatibility
                with pyteomics indexing.

        Returns
        -------
        pooledPeaks : list
            List of pooled peaks. See _format_raw_library_object_into_csodiaq_library_dict
                return value for 'peak' in libraryLoaderStrategy.py.
        """
        pass


def precursor_mz_missing_warning_text(scanNum):
    return f"scan number {scanNum} has no precursorMz and will be ignored. This may be because it is from an ms1 scan."
