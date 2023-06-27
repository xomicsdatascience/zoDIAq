from csodiaq.loaders.LibraryLoaderStrategyTable import LibraryLoaderStrategyTable as Table
from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf as Mgf
from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import os

class LibraryLoaderContext:
    def __init__(self, libraryFilePath: os.PathLike):
        self._libraryFilePath = libraryFilePath
        if libraryFilePath.endswith('.tsv'):
            self._strategy = Table()
        if libraryFilePath.endswith('.mgf'):
            self._strategy = Mgf()

    def load_csodiaq_library_dict(self):
        #NOTE: include max peaks variable?
        return self._strategy.load_csodiaq_library_dict_from_file(self._libraryFilePath)


