import pytest
import os
import pandas as pd
import re
from tempfile import TemporaryDirectory, NamedTemporaryFile
from zodiaq.loaders.library.libraryLoaderStrategyTable import (
    LibraryLoaderStrategyTable,
    _reformat_raw_library_object_columns,
    _organize_data_by_zodiaq_library_dict_keys,
    _determine_library_source_from_file,
)


@pytest.fixture
def loader():
    return LibraryLoaderStrategyTable()


def test__library_loader_strategy_table__initialization(loader):
    pass


def get_parent_dir():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture
def libFilePath():
    return os.path.join(
        get_parent_dir(), "test_files", "sample_lib_table_spectrast.tsv"
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__spectrast_library(
    loader, libFilePath
):
    loader._load_raw_library_object_from_file(libFilePath)
    assert hasattr(loader, "rawUploadedLibraryObject")
    assert isinstance(loader.rawUploadedLibraryObject, pd.DataFrame)


def check_value_error_thrown_when_missing_columns(loader, dfPath, missingColumnValues):
    if dfPath.endswith(".tsv"):
        separator = "\t"
    else:
        separator = ","
    libDf = pd.read_csv(dfPath, sep=separator)
    libDf.drop(missingColumnValues, axis=1, inplace=True)
    invalidLibFile = NamedTemporaryFile(
        prefix=f'zodiaq_table_loader_missing_column_{"_".join(missingColumnValues)}_',
        suffix=".tsv",
    )
    libDf.to_csv(invalidLibFile.name, sep="\t", index=False)
    errorOutput = f'table library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        loader._load_raw_library_object_from_file(invalidLibFile.name)


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__PrecursorMz(
    loader, libFilePath
):
    missingColumnValues = ["PrecursorMz"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__FullUniModPeptideName(
    loader, libFilePath
):
    missingColumnValues = ["FullUniModPeptideName"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__PrecursorCharge(
    loader, libFilePath
):
    missingColumnValues = ["PrecursorCharge"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__ProductMz(
    loader, libFilePath
):
    missingColumnValues = ["ProductMz"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__LibraryIntensity(
    loader, libFilePath
):
    missingColumnValues = ["LibraryIntensity"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__ProteinName(
    loader, libFilePath
):
    missingColumnValues = ["ProteinName"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__spectrast_library__2_columns_missing(
    loader, libFilePath
):
    missingColumnValues = ["ProteinName", "LibraryIntensity"]
    check_value_error_thrown_when_missing_columns(
        loader, libFilePath, missingColumnValues
    )


def assert_final_dict_output_matches_expected(outputDict, expectedOutputDict):
    for zodiaqLibKey in expectedOutputDict:
        assert zodiaqLibKey in outputDict
        for libEntryKey in expectedOutputDict[zodiaqLibKey]:
            assert libEntryKey in outputDict[zodiaqLibKey]
            assert (
                outputDict[zodiaqLibKey][libEntryKey]
                == expectedOutputDict[zodiaqLibKey][libEntryKey]
            )


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__spectrast_library(
    loader, libFilePath
):
    loader._load_raw_library_object_from_file(libFilePath)
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (516.801083027, "YRPGTVALR"): {
            "precursorCharge": 2,
            "identification": "51327_YRPGTVALR_2",
            "proteinName": "5/sp|Q71DI3|H32_HUMAN/sp|Q6NXT2|H3C_HUMAN/sp|Q16695|H31T_HUMAN/sp|P84243|H33_HUMAN/sp|P68431|H31_HUMAN",
            "peaks": [
                (429.745245163, 39.7, 0),
                (435.269418758, 84.0, 0),
                (458.308543918, 61.0, 0),
                (559.3562223920001, 67.3, 0),
                (575.293622136, 57.3, 0),
                (616.3776861179999, 124.8, 0),
                (674.362036054, 143.6, 0),
                (713.4304499719999, 1795.5, 0),
                (745.399149844, 324.1, 0),
                (858.483213826, 271.9, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        }
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__multi_peptide_spectrast_library(
    loader,
):
    libFilePath = os.path.join(
        get_parent_dir(), "test_files", "sample_lib_table_spectrast_multiple.tsv"
    )
    loader._load_raw_library_object_from_file(libFilePath)
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (300.83985497, "FVVGSHVR"): {
            "precursorCharge": 3,
            "identification": "1_FVVGSHVR_3",
            "proteinName": "1/sp|P49736|MCM2_HUMAN",
            "peaks": [
                (411.24627802, 1499.7, 0),
                (490.26601039, 107.9, 0),
                (498.27830643, 1696.3, 0),
                (555.299770156, 10000.0, 0),
                (654.368184074, 1368.4, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        },
        (300.83985497, "SVHGVVFR"): {
            "precursorCharge": 3,
            "identification": "2_SVHGVVFR_3",
            "proteinName": "1/DECOY_0_sp|P49736|MCM2_HUMAN",
            "peaks": [
                (421.255780064, 1499.7, 1),
                (480.2565083459999, 107.9, 1),
                (520.324193982, 1696.3, 1),
                (577.345657708, 10000.0, 1),
                (714.404569582, 1368.4, 1),
            ],
            "zodiaqKeyIdx": 1,
            "isDecoy": 1,
        },
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


def test__library_loader_strategy_table__reformat_raw_library_object_columns():
    numColumns = 10
    zodiaqKeyColumns = ["precursorMz", "peptideName"]
    oldMappedColumns = [str(i) for i in range(numColumns)]
    newMappedColumns = [columnName + "_new" for columnName in oldMappedColumns]
    oldMappedColumns += zodiaqKeyColumns
    newMappedColumns += zodiaqKeyColumns
    oldToNewColumnDict = dict(zip(oldMappedColumns, newMappedColumns))
    superfluousColumns = ["random", "superfluous", "columns"]
    oldColumns = oldMappedColumns + superfluousColumns
    newColumns = newMappedColumns + ["zodiaqLibKey"]
    data = [[0 for i in range(len(oldColumns))]]
    df = pd.DataFrame(data, columns=oldColumns)
    newDf = _reformat_raw_library_object_columns(df, oldToNewColumnDict)
    assert set(newDf.columns) == set(newColumns)


def test__library_loader_strategy_table__organize_data_by_zodiaq_library_dict_keys(
    loader, libFilePath
):
    loader._load_raw_library_object_from_file(libFilePath)
    reformattedDf = _reformat_raw_library_object_columns(
        loader.rawUploadedLibraryObject, loader.oldToNewColumnDict
    )
    df = reformattedDf[['peakMz','peakIntensity','fragmentType','fragmentNumber']]
    expectedKeys = [(516.801083027, "YRPGTVALR")]

    expectedTupleToListMzDict = {
        (516.801083027, "YRPGTVALR"): [
            713.4304499719999,
            745.399149844,
            858.483213826,
            674.362036054,
            616.3776861179999,
            435.269418758,
            559.3562223920001,
            458.308543918,
            575.293622136,
            429.745245163,
        ]
    }

    expectedTupleToListIntensityDict = {
        (516.801083027, "YRPGTVALR"): [
            1795.5,
            324.1,
            271.9,
            143.6,
            124.8,
            84.0,
            67.3,
            61.0,
            57.3,
            39.7,
        ]
    }

    expectedTupleToDictMetadataDict = {
        (516.801083027, "YRPGTVALR"): {
            "precursorCharge": 2,
            "identification": "51327_YRPGTVALR_2",
            "proteinName": "5/sp|Q71DI3|H32_HUMAN/sp|Q6NXT2|H3C_HUMAN/sp|Q16695|H31T_HUMAN/sp|P84243|H33_HUMAN/sp|P68431|H31_HUMAN",
        }
    }

    expectedTupleToFragmentTypeDf = {
        (516.801083027, "YRPGTVALR"): [
            ('y', 7),
            ('b', 7),
            ('b', 8),
            ('b', 6),
            ('y', 6),
            ('y', 8),
            ('y', 5),
            ('y', 4),
            ('b', 5),
            ('b', 8),
        ]

    }

    dataDict = _organize_data_by_zodiaq_library_dict_keys(reformattedDf)
    assert dataDict["zodiaqKeys"] == expectedKeys
    assert dataDict["mz"] == expectedTupleToListMzDict
    assert dataDict["intensities"] == expectedTupleToListIntensityDict
    assert dataDict["metadata"] == expectedTupleToDictMetadataDict
    assert dataDict["fragmentTypes"] == expectedTupleToFragmentTypeDf


@pytest.fixture
def fragpipeLibFilePath():
    return os.path.join(get_parent_dir(), "test_files", "sample_lib_table_fragpipe.tsv")


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__PrecursorMz(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["PrecursorMz"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__ModifiedPeptideSequence(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["ModifiedPeptideSequence"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__PrecursorCharge(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["PrecursorCharge"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__ProductMz(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["ProductMz"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__LibraryIntensity(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["LibraryIntensity"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__fragpipe_library__PeptideSequence(
    loader, fragpipeLibFilePath
):
    missingColumnValues = ["PeptideSequence"]
    check_value_error_thrown_when_missing_columns(
        loader, fragpipeLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__fragpipe_library(
    loader, fragpipeLibFilePath
):
    loader._load_raw_library_object_from_file(fragpipeLibFilePath)
    df = loader.rawUploadedLibraryObject
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (375.873226, "FANYIDKVR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "P08670",
            "peaks": [
                (175.118953, 2926.18, 0),
                (274.187367, 1647.689, 0),
                (333.155733, 1071.177, 0),
                (397.231972, 2078.822, 0),
                (402.282331, 4932.288, 0),
                (454.253437, 1301.4617, 0),
                (489.771994, 1395.553, 0),
                (517.309275, 10000.0, 0),
                (630.393339, 8233.006, 0),
                (793.456668, 5096.472, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
            "fragmentTypes": [
                ('y', 1),
                ('y', 2),
                ('b', 3),
                ('y', 6),
                ('y', 3),
                ('y', 7),
                ('y', 8),
                ('y', 4),
                ('y', 5),
                ('y', 6),
            ],
        }
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__multi_peptides_fragpipe_library(
    loader, fragpipeLibFilePath
):
    loader._load_raw_library_object_from_file(
        os.path.join(
            get_parent_dir(), "test_files", "sample_lib_table_fragpipe_multiple.tsv"
        )
    )
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (375.873226, "FANYIDKVR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "P08670",
            "peaks": [
                (175.118953, 2926.18, 0),
                (274.187367, 1647.689, 0),
                (333.155733, 1071.177, 0),
                (397.231972, 2078.822, 0),
                (402.282331, 4932.288, 0),
                (454.253437, 1301.4617, 0),
                (489.771994, 1395.553, 0),
                (517.309275, 10000.0, 0),
                (630.393339, 8233.006, 0),
                (793.456668, 5096.472, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        },
        (375.885354, "FGTINIVHPK"): {
            "precursorCharge": 3,
            "identification": "FGTINIVHPK",
            "proteinName": "Q9Y617",
            "peaks": [
                (205.097155, 1987.8866, 1),
                (244.165569, 10000.0, 1),
                (306.144834, 1832.419, 1),
                (381.224481, 3466.5054, 1),
                (480.292896, 3109.8425, 1),
                (490.271629, 801.70935, 1),
                (533.271827, 358.89835, 1),
                (593.37696, 2114.4893, 1),
                (707.419888, 5087.7124, 1),
                (820.503953, 865.17474, 1),
            ],
            "zodiaqKeyIdx": 1,
            "isDecoy": 0,
        },
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


@pytest.fixture
def prositLibFilePath():
    return os.path.join(get_parent_dir(), "test_files", "sample_lib_table_prosit.csv")


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__PrecursorMz(
    loader, prositLibFilePath
):
    missingColumnValues = ["PrecursorMz"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__ModifiedPeptide(
    loader, prositLibFilePath
):
    missingColumnValues = ["ModifiedPeptide"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__PrecursorCharge(
    loader, prositLibFilePath
):
    missingColumnValues = ["PrecursorCharge"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__FragmentMz(
    loader, prositLibFilePath
):
    missingColumnValues = ["FragmentMz"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__StrippedPeptide(
    loader, prositLibFilePath
):
    missingColumnValues = ["StrippedPeptide"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__load_raw_library_object_from_file__fails_when_missing_required_columns__prosit_library__FragmentLossType(
    loader, prositLibFilePath
):
    missingColumnValues = ["FragmentLossType"]
    check_value_error_thrown_when_missing_columns(
        loader, prositLibFilePath, missingColumnValues
    )


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__prosit_library(
    loader, prositLibFilePath
):
    loader._load_raw_library_object_from_file(prositLibFilePath)
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (374.1867597566666, "_MMPAAALIM[Oxidation (O)]R_"): {
            "precursorCharge": 3,
            "identification": "MMPAAALIMR",
            "proteinName": "noloss",
            "peaks": [
                (175.11895751953125, 1.0, 0),
                (263.0882568359375, 0.1923068910837173, 0),
                (322.15435791015625, 0.412596195936203, 0),
                (360.1410217285156, 0.085973247885704, 0),
                (431.1781311035156, 0.1523399353027343, 0),
                (435.2384033203125, 0.7306222915649414, 0),
                (502.2152404785156, 0.0825881585478782, 0),
                (548.322509765625, 0.3042449355125427, 0),
                (619.359619140625, 0.1164016127586364, 0),
                (690.396728515625, 0.0937163606286048, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        }
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


def test__library_loader_strategy_table__format_raw_library_object_into_zodiaq_library_dict__multi_peptide_prosit_library(
    loader, prositLibFilePath
):
    loader._load_raw_library_object_from_file(
        os.path.join(
            get_parent_dir(), "test_files", "sample_lib_table_prosit_multiple.csv"
        )
    )
    outputDict = loader._format_raw_library_object_into_zodiaq_library_dict()
    expectedOutputDict = {
        (254.3121828783333, "_MRALLLIPPPPM[Oxidation (O)]R_"): {
            "precursorCharge": 6,
            "identification": "MRALLLIPPPPMR",
            "proteinName": "noloss",
            "peaks": [
                (175.11895751953125, 0.6478143930435181, 0),
                (258.6335754394531, 0.3157764971256256, 0),
                (288.1488647460937, 0.985026240348816, 0),
                (307.15997314453125, 0.6946653723716736, 0),
                (355.68634033203125, 0.4260823428630829, 0),
                (359.1859741210937, 0.4684059023857116, 0),
                (419.2071228027344, 0.7623715996742249, 0),
                (516.2598876953125, 1.0, 0),
                (585.3541259765625, 0.5721365809440613, 0),
                (613.3126220703125, 0.7858051657676697, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        },
        (374.1867597566666, "_MMPAAALIM[Oxidation (O)]R_"): {
            "precursorCharge": 3,
            "identification": "MMPAAALIMR",
            "proteinName": "noloss",
            "peaks": [
                (175.11895751953125, 1.0, 1),
                (263.0882568359375, 0.1923068910837173, 1),
                (322.15435791015625, 0.412596195936203, 1),
                (360.1410217285156, 0.085973247885704, 1),
                (431.1781311035156, 0.1523399353027343, 1),
                (435.2384033203125, 0.7306222915649414, 1),
                (502.2152404785156, 0.0825881585478782, 1),
                (548.322509765625, 0.3042449355125427, 1),
                (619.359619140625, 0.1164016127586364, 1),
                (690.396728515625, 0.0937163606286048, 1),
            ],
            "zodiaqKeyIdx": 1,
            "isDecoy": 0,
        },
        (507.272473135, "_MLAPPPIM[Oxidation (O)]K_"): {
            "precursorCharge": 2,
            "identification": "MLAPPPIMK",
            "proteinName": "noloss",
            "peaks": [
                (245.1318206787109, 0.8553705215454102, 2),
                (294.148193359375, 0.0534038245677948, 2),
                (301.17254638671875, 0.1253907978534698, 2),
                (316.1689453125, 0.1860230714082718, 2),
                (349.69891357421875, 0.6669084429740906, 2),
                (385.2174682617188, 0.0563227161765098, 2),
                (504.2850341796875, 0.1659155339002609, 2),
                (601.3377685546875, 0.6199196577072144, 2),
                (698.3905639648438, 1.0, 2),
                (769.4276733398438, 0.3755797147750854, 2),
            ],
            "zodiaqKeyIdx": 2,
            "isDecoy": 0,
        },
    }
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)


def test__library_loader_strategy_table__determine_library_source_from_file(
    libFilePath, fragpipeLibFilePath, prositLibFilePath
):
    assert _determine_library_source_from_file(libFilePath) == "spectrast"
    assert _determine_library_source_from_file(fragpipeLibFilePath) == "fragpipe"
    assert _determine_library_source_from_file(prositLibFilePath) == "prosit"
    errorOutput = "The library table file provided does not match spectrast, fragpipe, or prosit formats."
    with pytest.raises(ValueError, match=re.escape(errorOutput)):
        _determine_library_source_from_file(
            os.path.join(get_parent_dir(), "test_files", "sample_lib_mgf.mgf")
        )
