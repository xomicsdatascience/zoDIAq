from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import pandas as pd
import os

class LibraryLoaderStrategyTraml(LibraryLoaderStrategy):

    def __init__(self):
        self.oldToNewColumnDict = {
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'FullUniModPeptideName',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'transitionGroupId',
            'ProteinName': 'ProteinName'
        }

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawLibDf = pd.read_csv(libraryFilePath, sep='\t')
        assert_there_are_no_missing_columns(self.oldToNewColumnDict.keys(), self.rawLibDf.columns)

    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict:
        self.rawLibDf = self.rawLibDf[self.oldToNewColumnDict.keys()]
        remap_table_columns(self.rawLibDf, self.oldToNewColumnDict)
        tupleToListMzDict, tupleToListIntensityDict, tupleToDictDataDict = create_data_dicts_that_correspond_to_csodiaq_library_dict_keys(
            self.rawLibDf)
        sortedCsodiaqKeys = sorted(tupleToDictDataDict.keys())
        for csodiaqKeyIdx in range(len(sortedCsodiaqKeys)):
            csodiaqKey = sortedCsodiaqKeys[csodiaqKeyIdx]
            mzList = tupleToListMzDict[csodiaqKey]
            intensityList = tupleToListIntensityDict[csodiaqKey]
            peaks = create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(mzList, intensityList, csodiaqKeyIdx)
        pass

def assert_there_are_no_missing_columns(requiredColumns: list, presentColumns: list):
    missingColumnValues = set(requiredColumns) - set(presentColumns)
    if len(missingColumnValues):
        raise ValueError(
            f'traml spectrast library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
        )

def remap_table_columns(df: pd.DataFrame, oldToNewColumnDict: dict):
    df.rename(columns=oldToNewColumnDict, inplace=True)

def create_csodiaq_library_dict_keys_as_new_column(df: pd.DataFrame):
    df['precursorMzAndPeptide'] = list(zip(df['PrecursorMz'].tolist(),
                          df['FullUniModPeptideName'].tolist()))

def create_data_dicts_that_correspond_to_csodiaq_library_dict_keys(df: pd.DataFrame):
    create_csodiaq_library_dict_keys_as_new_column(df)
    mz = df.groupby('precursorMzAndPeptide')['ProductMz'].apply(list).to_dict()
    intensities = df.groupby('precursorMzAndPeptide')['LibraryIntensity'].apply(list).to_dict()
    df.drop_duplicates(subset='precursorMzAndPeptide', inplace=True)
    df.set_index('precursorMzAndPeptide', drop=True, inplace=True)
    data = df.to_dict(orient='index')
    return mz,intensities,data

def create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(mzList, intensityList, id):
    idList = [id for i in range(len(mzList))]
    return list(zip(mzList, intensityList, idList))
