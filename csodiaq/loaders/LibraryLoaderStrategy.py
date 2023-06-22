from abc import ABC, abstractmethod

class LibraryLoaderStrategy(ABC):

    @abstractmethod
    def load_raw_library_file(self, libraryFilePath) -> None: pass

    @abstractmethod
    def format_raw_library_file() -> dict: pass

    def load_library_file(self, libraryFilePath) -> dict:
        self.load_raw_library_file(libraryFilePath)
        return self.format_raw_library_file()
