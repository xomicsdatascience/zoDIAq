from zodiaq.loaders.query.queryLoaderStrategyMzxml import (
    QueryLoaderStrategyMzxml as Mzxml,
)
import os


class QueryLoaderContext:
    """
    Class for accessing query loading strategies (as used by the strategy design pattern).
    """

    def __init__(self, queryFilePath: os.PathLike):
        self.filePath = queryFilePath
        self._strategy = Mzxml(queryFilePath)

    def count_number_of_ms1_and_ms2_scans_per_cycle(self):
        return self._strategy.count_number_of_ms1_and_ms2_scans_per_cycle()

    def extract_metadata_from_query_scans(self) -> dict:
        return self._strategy.extract_metadata_from_query_scans()

    def get_query_file_reader(self):
        return self._strategy.get_query_file_reader()

    def pool_peaks_of_query_scans(self, scans: list, reader) -> list:
        return self._strategy.pool_peaks_of_query_scans(scans, reader)

    def identify_mz_window_highest_and_lowest_bound(self, scans: list, reader) -> list:
        return self._strategy.identify_mz_window_highest_and_lowest_bound(scans, reader)
