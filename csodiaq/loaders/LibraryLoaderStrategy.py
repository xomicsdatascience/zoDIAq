from abc import ABC, abstractmethod
import os

class LibraryLoaderStrategy(ABC):

    @abstractmethod
    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None: pass

    @abstractmethod
    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict: pass

    #NOTE: test with each of your test files?
    def load_csodiaq_library_dict_from_file(self, libraryFilePath: os.PathLike) -> dict:
        self._load_raw_library_object_from_file(libraryFilePath)
        return self._format_raw_library_object_into_csodiaq_library_dict()

def create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(mz: list, intensities: list, id: int) -> list:
    idList = [id for i in range(len(mz))]
    return list(zip(mz, intensities, idList))

def remove_low_intensity_peaks_below_max_peak_num(peaks: list, maxPeakNum: int) -> list:
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks[:maxPeakNum]