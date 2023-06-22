from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import pandas as pd
import os

class LibraryLoaderStrategyTraml(LibraryLoaderStrategy):

    def load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawLibDf = pd.read_csv(libraryFilePath, sep='\t')
        requiredSpectrastColumnNames = [
            'PrecursorMz',
            'FullUniModPeptideName',
            'PrecursorCharge',
            'ProductMz',
            'LibraryIntensity',
            'transition_group_id',
            'ProteinName',
        ]
        missingColumnValues = set(requiredSpectrastColumnNames) - set(self.rawLibDf.columns)
        if len(missingColumnValues):
            raise ValueError(f'traml spectrast library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])')
        self.rawLibDf = self.rawLibDf[requiredSpectrastColumnNames]

    def format_raw_library_object_into_csodiaq_dict() -> dict: pass