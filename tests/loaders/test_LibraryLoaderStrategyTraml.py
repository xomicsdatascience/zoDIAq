import pytest
import os
import pandas as pd
import re
from tempfile import TemporaryDirectory, NamedTemporaryFile
from csodiaq.loaders.LibraryLoaderStrategyTraml import LibraryLoaderStrategyTraml, remap_table_columns, create_csodiaq_library_dict_keys_as_new_column, create_data_dicts_that_correspond_to_csodiaq_library_dict_keys, create_peaks_from_mz_intensity_lists_and_csodiaq_key_id

@pytest.fixture
def loader():
    return LibraryLoaderStrategyTraml()

def test__library_loader_strategy_traml__initialization(loader): pass

@pytest.fixture
def tramlTestLibFilePath():
    cwd = os.getcwd()
    parent = os.path.dirname(cwd)
    return os.path.join(parent,  'test_files', 'spectra_st.tsv')

def test__library_loader_strategy_traml__load_raw_library_object_from_file__spectrast_library(loader, tramlTestLibFilePath):
    loader._load_raw_library_object_from_file(tramlTestLibFilePath)
    assert hasattr(loader, 'rawLibDf')
    assert isinstance(loader.rawLibDf, pd.DataFrame)

def check_value_error_thrown_when_missing_columns(loader, dfPath, missingColumnValues):
    libDf = pd.read_csv(dfPath, sep='\t')
    libDf.drop(missingColumnValues, axis=1, inplace=True)
    invalidLibFile = NamedTemporaryFile(prefix=f'csodiaq_traml_loader_missing_column_{"_".join(missingColumnValues)}_', suffix=".tsv")
    libDf.to_csv(invalidLibFile.name, sep='\t', index=False)
    errorOutput = f'traml spectrast library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        loader._load_raw_library_object_from_file(invalidLibFile.name)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__PrecursorMz(loader, tramlTestLibFilePath):
    missingColumnValues = ['PrecursorMz']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__FullUniModPeptideName(loader, tramlTestLibFilePath):
    missingColumnValues = ['FullUniModPeptideName']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__PrecursorCharge(loader, tramlTestLibFilePath):
    missingColumnValues = ['PrecursorCharge']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__ProductMz(loader, tramlTestLibFilePath):
    missingColumnValues = ['ProductMz']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__LibraryIntensity(loader, tramlTestLibFilePath):
    missingColumnValues = ['LibraryIntensity']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__transition_group_id(loader, tramlTestLibFilePath):
    missingColumnValues = ['transition_group_id']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__ProteinName(loader, tramlTestLibFilePath):
    missingColumnValues = ['ProteinName']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__2_columns_missing(loader, tramlTestLibFilePath):
    missingColumnValues = ['ProteinName', 'LibraryIntensity']
    check_value_error_thrown_when_missing_columns(loader, tramlTestLibFilePath, missingColumnValues)

def test__library_loader_strategy_traml__format_raw_library_object_into_csodiaq_library_dict__spectrast_library(loader, tramlTestLibFilePath):
    loader._load_raw_library_object_from_file(tramlTestLibFilePath)
    outputDict = loader._format_raw_library_object_into_csodiaq_library_dict()
    expectedOutputDict = {(375.87322429333335, 'FANYIDKVR'): {'transition_group_id': '1_FANYIDKVR_3', 'proteinName': '1/sp|P08670|VIME_HUMAN', 'precursorCharge': 3, 'peaks': [(517.3092722, 10000.0, 1), (630.393336182, 8235.6, 1), (793.4566647199999, 5098.5, 1), (402.2823291699999, 4930.4, 1), (397.23197061, 2082.7, 1), (489.771991231, 1398.4, 1), (454.253434336, 1301.3, 1), (333.15573166, 1072.6, 1), (500.2792722, 863.7, 1), (445.738434336, 746.7, 1)], 'csodiaqKeyID': 1, 'isDecoy': 0}}


def test__library_loader_strategy_traml__remap_table_columns():
    numColumns = 10
    oldColumns = [str(i) for i in range(numColumns)]
    newColumns = [columnName + '_new' for columnName in oldColumns]
    data = [
        [0 for i in range(numColumns)]
    ]
    df = pd.DataFrame(data, columns = oldColumns)
    oldToNewColumnDict = dict(zip(oldColumns, newColumns))
    remap_table_columns(df, oldToNewColumnDict)
    assert set(df.columns) == set(newColumns)

def test__library_loader_strategy_traml__create_csodiaq_library_dict_keys_as_new_column(loader, tramlTestLibFilePath):
    loader._load_raw_library_object_from_file(tramlTestLibFilePath)
    create_csodiaq_library_dict_keys_as_new_column(loader.rawLibDf)
    expectedNewColumn = [(375.87322429333335, 'FANYIDKVR') for i in range(len(loader.rawLibDf))]
    assert list(loader.rawLibDf['precursorMzAndPeptide']) == expectedNewColumn

def test__library_loader_strategy_traml__create_data_dicts_that_correspond_to_csodiaq_library_dict_keys(loader, tramlTestLibFilePath):
    loader._load_raw_library_object_from_file(tramlTestLibFilePath)
    create_csodiaq_library_dict_keys_as_new_column(loader.rawLibDf)
    expectedTupleToListMzDict = {(375.87322429333335, 'FANYIDKVR'): [517.3092722, 630.393336182, 793.4566647199999, 402.2823291699999, 397.23197061, 489.771991231, 454.253434336, 333.15573166, 500.2792722, 445.738434336]}
    expectedTupleToListIntensityDict = {(375.87322429333335, 'FANYIDKVR'): [10000.0, 8235.6, 5098.5, 4930.4, 2082.7, 1398.4, 1301.3, 1072.6, 863.7, 746.7]}
    expectedTupleToDictMetadataDict = {(375.87322429333335, 'FANYIDKVR'): {'PrecursorMz': 375.87322429333335, 'ProductMz': 517.3092722, 'Tr_recalibrated': 1107.8, 'transition_name': '18_y4_1_FANYIDKVR_3', 'CE': -1, 'LibraryIntensity': 10000.0, 'transition_group_id': '1_FANYIDKVR_3', 'decoy': 0, 'PeptideSequence': 'FANYIDKVR', 'ProteinName': '1/sp|P08670|VIME_HUMAN', 'Annotation': 'y4/-0.000', 'FullUniModPeptideName': 'FANYIDKVR', 'PrecursorCharge': 3, 'PeptideGroupLabel': '1_FANYIDKVR_3', 'UniprotID': '1/sp|P08670|VIME_HUMAN', 'FragmentType': 'y', 'FragmentCharge': 1, 'FragmentSeriesNumber': 4, 'LabelType': 'light'}}
    tupleToListMzDict, tupleToListIntensityDict, tupleToDictMetadataDict = create_data_dicts_that_correspond_to_csodiaq_library_dict_keys(loader.rawLibDf)
    assert tupleToListMzDict == expectedTupleToListMzDict
    assert tupleToListIntensityDict == expectedTupleToListIntensityDict
    assert tupleToDictMetadataDict == expectedTupleToDictMetadataDict

def test__library_loader_strategy_traml__create_peaks_from_mz_intensity_lists_and_csodiaq_key_id():
    numPeaks = 15
    mzList = [i for i in range(numPeaks-1,-1,-1)]
    intensityList = [i for i in range(numPeaks)]
    id = 0
    expectedPeaks = [
        (14, 0, 0),
        (13, 1, 0),
        (12, 2, 0),
        (11, 3, 0),
        (10, 4, 0),
        (9, 5, 0),
        (8, 6, 0),
        (7, 7, 0),
        (6, 8, 0),
        (5, 9, 0),
        (4, 10, 0),
        (3, 11, 0),
        (2, 12, 0),
        (1, 13, 0),
        (0, 14, 0),
    ]

    peaks = create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(mzList, intensityList, id)
    assert peaks == expectedPeaks

