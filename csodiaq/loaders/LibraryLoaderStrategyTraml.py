from csodiaq.loaders.LibraryLoaderStrategy import LibraryLoaderStrategy
import pandas as pd
import os

class LibraryLoaderStrategyTraml(LibraryLoaderStrategy):

    def __init__(self, librarySource: str):
        self.oldToNewColumnDict = set_old_to_new_column_dict(librarySource)

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        self.rawLibDf = pd.read_csv(libraryFilePath, sep='\t')
        assert_there_are_no_missing_columns(self.oldToNewColumnDict.keys(), self.rawLibDf.columns)

    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict:
        maxPeakNum = 10
        reformattedLibDf = reformat_raw_library_object_columns(self.rawLibDf, self.oldToNewColumnDict)
        organizedDataDict = organize_data_by_csodiaq_library_dict_keys(reformattedLibDf)
        csodiaqLibraryDict = {}
        for csodiaqKeyIdx in range(len(organizedDataDict['csodiaqKeys'])):
            csodiaqKey = organizedDataDict['csodiaqKeys'][csodiaqKeyIdx]
            csodiaqLibraryDict[csodiaqKey] = create_csodiaq_library_entry(organizedDataDict, maxPeakNum, csodiaqKeyIdx)
        return csodiaqLibraryDict

def assert_there_are_no_missing_columns(requiredColumns: list, presentColumns: list) -> None:
    missingColumnValues = set(requiredColumns) - set(presentColumns)
    if len(missingColumnValues):
        raise ValueError(
            f'traml spectrast library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
        )

def reformat_raw_library_object_columns(df: pd.DataFrame, oldToNewColumnDict: dict) -> pd.DataFrame:
    reformattedDf = df[oldToNewColumnDict.keys()]
    reformattedDf = reformattedDf.rename(columns=oldToNewColumnDict)
    reformattedDf['csodiaqLibKey'] = list(zip(reformattedDf['precursorMz'].tolist(),
                          reformattedDf['peptideName'].tolist()))
    return reformattedDf

def organize_data_by_csodiaq_library_dict_keys(df: pd.DataFrame) -> dict:
    keys = sorted(set(df['csodiaqLibKey']))
    mz = df.groupby('csodiaqLibKey')['peakMz'].apply(list).to_dict()
    intensities = df.groupby('csodiaqLibKey')['peakIntensity'].apply(list).to_dict()
    df.drop_duplicates(subset='csodiaqLibKey', inplace=True)
    df.set_index('csodiaqLibKey', drop=True, inplace=True)
    df.drop(['precursorMz','peakMz','peptideName','peakIntensity'], axis=1, inplace=True)
    metadata = df.to_dict(orient='index')
    return {
        'csodiaqKeys': keys,
        'mz':mz,
        'intensities':intensities,
        'metadata':metadata
    }

def create_csodiaq_library_entry(organizedDataDict: dict, maxPeakNum: int, csodiaqKeyIdx: int) -> dict:
    csodiaqKey = organizedDataDict['csodiaqKeys'][csodiaqKeyIdx]
    peaks = create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(organizedDataDict['mz'][csodiaqKey],
                                                                    organizedDataDict['intensities'][csodiaqKey],
                                                                    csodiaqKeyIdx)
    reducedPeaks = remove_low_intensity_peaks_below_max_peak_num(peaks, maxPeakNum)
    isDecoy = int('decoy' in organizedDataDict['metadata'][csodiaqKey]['proteinName'].lower())
    return {
                'precursorCharge': organizedDataDict['metadata'][csodiaqKey]['precursorCharge'],
                'identifier': organizedDataDict['metadata'][csodiaqKey]['identifier'],
                'proteinName': organizedDataDict['metadata'][csodiaqKey]['proteinName'],
                'peaks': sorted(reducedPeaks),
                'csodiaqKeyIdx': csodiaqKeyIdx,
                'isDecoy': isDecoy,
            }

def create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(mz: list, intensities: list, id: int) -> list:
    idList = [id for i in range(len(mz))]
    return list(zip(mz, intensities, idList))

def remove_low_intensity_peaks_below_max_peak_num(peaks: list, maxPeakNum: int) -> list:
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks[:maxPeakNum]

def set_old_to_new_column_dict(librarySource):
    newColumns = [
        'precursorMz',
        'peptideName',
        'peakMz',
        'peakIntensity',
        'precursorCharge',
        'identifier',
        'proteinName',
    ]
    if librarySource == 'fragpipe':
        oldColumns = [
            'PrecursorMz',
            'ModifiedPeptideSequence',
            'ProductMz',
            'LibraryIntensity',
            'PrecursorCharge',
            'PeptideSequence',
            'ProteinId',
        ]
        return dict(zip(oldColumns, newColumns))
    else:
        oldColumns = [
            'PrecursorMz',
            'FullUniModPeptideName',
            'ProductMz',
            'LibraryIntensity',
            'PrecursorCharge',
            'transition_group_id',
            'ProteinName',
        ]
        return dict(zip(oldColumns, newColumns))