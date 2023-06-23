import pytest
import os
import pandas as pd
import re
from tempfile import TemporaryDirectory, NamedTemporaryFile
from csodiaq.loaders.LibraryLoaderStrategyTraml import LibraryLoaderStrategyTraml, remap_table_columns

@pytest.fixture
def loader():
    return LibraryLoaderStrategyTraml()

def test__library_loader_strategy_traml__initialization(loader): pass

@pytest.fixture
def tramlTestLibFilePath():
    cwd = os.getcwd()
    parent = os.path.dirname(cwd)
    return os.path.join(parent,  'test_files', 'spectra_st.tsv')

def test__library_loader_strategy_traml__load_raw_library_object_from_file(loader, tramlTestLibFilePath):
    loader._load_raw_library_object_from_file(tramlTestLibFilePath)
    assert hasattr(loader, 'rawLibDf')
    assert isinstance(loader.rawLibDf, pd.DataFrame)
    '''
    # I want to include this under the formatting function instead
    assert loader.rawLibDf.shape[1] == 7
    expectedColumnsSet = set([
        'PrecursorMz',
        'FullUniModPeptideName',
        'PrecursorCharge',
        'ProductMz',
        'LibraryIntensity',
        'transitionGroupId',
        'ProteinName',
    ])
    assert set(loader.rawLibDf.columns) == expectedColumnsSet
    '''

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
