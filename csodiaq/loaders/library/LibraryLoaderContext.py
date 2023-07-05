from csodiaq.loaders.library.LibraryLoaderStrategyTable import (
    LibraryLoaderStrategyTable as Table,
)
from csodiaq.loaders.library.LibraryLoaderStrategyMgf import (
    LibraryLoaderStrategyMgf as Mgf,
)
from csodiaq.loaders.library.LibraryLoaderStrategy import LibraryLoaderStrategy
import os


class LibraryLoaderContext:
    """
    Class for accessing library loading strategies (as used by the strategy design pattern).

    Attributes
    ----------
    _libraryFilePath : string (os.PathLike format)
        The path to the library file. The file type will indicate which concrete library loading
        strategy class will be used.
    _strategy : LibraryLoaderStrategy
        The concrete library loading strategy class.
    """

    def __init__(self, libraryFilePath: os.PathLike):
        self._libraryFilePath = libraryFilePath
        if libraryFilePath.endswith(".tsv") or libraryFilePath.endswith(".csv"):
            self._strategy = Table()
        elif libraryFilePath.endswith(".mgf"):
            self._strategy = Mgf()
        elif libraryFilePath.endswith(".msp"):
            raise ValueError(
                "The .msp library format is not currently supported. If the library file was generated via prosit, please reset the output into a tab-delimited (.tsv) format."
            )

    def load_csodiaq_library_dict(self):
        """See 'load_csodiaq_library_dict_from_file' in 'LibraryLoaderStrategy'"""
        return self._strategy.load_csodiaq_library_dict_from_file(self._libraryFilePath)
