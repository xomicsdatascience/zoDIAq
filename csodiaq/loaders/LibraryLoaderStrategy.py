from abc import ABC, abstractmethod
import os

class LibraryLoaderStrategy(ABC):

    @abstractmethod
    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None: pass

    @abstractmethod
    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict: pass

    #NOTE: test with each of your test files?
    def load_csodiaq_library_dict_from_file(self, libraryFilePath: os.PathLike) -> dict:
        self.load_raw_library_object_from_file(libraryFilePath)
        return self.format_raw_library_object_into_csodiaq_dict()
