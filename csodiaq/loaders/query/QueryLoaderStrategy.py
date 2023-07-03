from abc import ABC, abstractmethod
import os

class QueryLoaderStrategy(ABC):
    def __init__(self, queryFilePath: os.PathLike):
        self.filePath = queryFilePath

    @abstractmethod
    def map_query_scan_ids_to_dia_mz_windows(self) -> dict: pass

    @abstractmethod
    def extract_metadata_from_query_scans(self) -> dict: pass

    @abstractmethod
    def pool_peaks_of_query_scans(self, scans: list) -> list: pass

