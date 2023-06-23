from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import pandas as pd
import os

class LibraryLoaderStrategyTraml(LibraryLoaderStrategy):

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawLibDf = pd.read_csv(libraryFilePath, sep='\t')
        oldToNewColumnDict = { #NOTE: set in __init__ depending on type?
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'FullUniModPeptideName',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'transitionGroupId',
            'ProteinName': 'ProteinName'
        }
        assert_there_are_no_missing_columns(oldToNewColumnDict.keys(), self.rawLibDf.columns)

    def _format_raw_library_object_into_csodiaq_dict() -> dict: pass
        #self.rawLibDf = self.rawLibDf[oldToNewColumnDict.keys()]
        #remap_table_columns

def assert_there_are_no_missing_columns(requiredColumns, presentColumns):
    missingColumnValues = set(requiredColumns) - set(presentColumns)
    if len(missingColumnValues):
        raise ValueError(
            f'traml spectrast library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
        )

def remap_table_columns(df, oldToNewColumnDict):
    df.rename(columns=oldToNewColumnDict, inplace=True)
