from csodiaq.loaders.LibraryLoaderStrategyTraml import LibraryLoaderStrategyTraml as Traml
from csodiaq.loaders.LibraryLoaderStrategyMgf import LibraryLoaderStrategyMgf as Mgf
from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import os

class LibraryLoaderContext:
    def __init__(self, libraryFilePath: os.PathLike):
        self._libraryFilePath = libraryFilePath
        if libraryFilePath.endswith('.tsv'):
            self._strategy = Traml('spectrast')
        if libraryFilePath.endswith('.mgf'):
            self._strategy = Mgf()

    def load_csodiaq_library_dict(self):
        return self._strategy.load_csodiaq_library_dict_from_file(self._libraryFilePath)


