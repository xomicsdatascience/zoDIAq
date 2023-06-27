from csodiaq.loaders.LibraryLoaderStrategyTable import LibraryLoaderStrategyTable as Table
from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf as Mgf
from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import os

class LibraryLoaderContext:
    def __init__(self, libraryFilePath: os.PathLike):
        self._libraryFilePath = libraryFilePath
        if libraryFilePath.endswith('.tsv') or libraryFilePath.endswith('.csv'):
            self._strategy = Table()
        elif libraryFilePath.endswith('.mgf'):
            self._strategy = Mgf()
        elif libraryFilePath.endswith('.msp'):
            raise ValueError('The .msp library format is not currently supported. If the library file was generated via prosit, please reset the output into a tab-delimited (.tsv) format.')

    def load_csodiaq_library_dict(self):
        #NOTE: include max peaks variable?
        return self._strategy.load_csodiaq_library_dict_from_file(self._libraryFilePath)


