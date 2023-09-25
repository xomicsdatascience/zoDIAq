from zodiaq.utils import (
    create_outfile_header,
    confirm_proteins_in_list_are_in_appropriate_format,
)
import pandas as pd
import pytest


@pytest.fixture
def outputDirectory():
    return "test/output/dir"


@pytest.fixture
def inputFileName():
    return "mzxml_test"


@pytest.fixture
def inputFile(inputFileName):
    return inputFileName + ".mzxml"


@pytest.fixture
def inputFilePath(inputFile):
    return "mzxml/directory/" + inputFile


@pytest.fixture
def outputZodiaqTag():
    return "zoDIAq-file_"


def test__utils__create_outfile_header__no_correction(
    outputDirectory, inputFileName, inputFilePath, outputZodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputZodiaqTag}{inputFileName}"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=-1)
    assert expectedOutput == output


def test__utils__create_outfile_header__no_correction__output_directory_ends_in_slash(
    outputDirectory, inputFileName, inputFilePath, outputZodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputZodiaqTag}{inputFileName}"
    output = create_outfile_header(outputDirectory + "/", inputFilePath, correction=-1)
    assert expectedOutput == output


def test__utils__create_outfile_header__no_correction__includes_non_file_type_dots(
    outputDirectory, inputFileName, inputFilePath, outputZodiaqTag
):
    inputFileNameWithPeriods = inputFileName + ".dots.added"
    inputFilePathWithPeriods = inputFileNameWithPeriods + ".mzxml"
    expectedOutput = f"{outputDirectory}/{outputZodiaqTag}{inputFileNameWithPeriods}"
    output = create_outfile_header(
        outputDirectory, inputFilePathWithPeriods, correction=-1
    )
    assert expectedOutput == output


def test__utils__create_outfile_header__custom_correction(
    outputDirectory, inputFileName, inputFilePath, outputZodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputZodiaqTag}{inputFileName}_corrected"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=0)
    assert expectedOutput == output


def test__utils__create_outfile_header__stdev_correction(
    outputDirectory, inputFileName, inputFilePath, outputZodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputZodiaqTag}{inputFileName}_corrected"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=1)
    assert expectedOutput == output


@pytest.fixture
def validProtein():
    return "2/protein1/protein2"


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__multiple_digits_returns_true(
    validProtein,
):
    proteinList = [
        validProtein,
        "10/protein1/protein2/protein3/protein4/protein5/protein6/protein7/protein8/protein9/protein10",
    ]
    assert confirm_proteins_in_list_are_in_appropriate_format(proteinList)


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__extra_prefix_returns_false(
    validProtein,
):
    proteinList = [
        validProtein,
        "extraStuff" + validProtein,
    ]
    assert not confirm_proteins_in_list_are_in_appropriate_format(proteinList)


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__ends_in_slash_returns_false(
    validProtein,
):
    proteinList = [
        validProtein,
        validProtein + "/",
    ]
    assert not confirm_proteins_in_list_are_in_appropriate_format(proteinList)


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__gibberish_returns_false(
    validProtein,
):
    proteinList = [
        validProtein,
        "gibberish",
    ]
    assert not confirm_proteins_in_list_are_in_appropriate_format(proteinList)


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__wrong_count_returns_false(
    validProtein,
):
    proteinList = [
        validProtein,
        "1/protein1/protein2",
    ]
    assert not confirm_proteins_in_list_are_in_appropriate_format(proteinList)


def test__utils__confirm_proteins_in_list_are_in_appropriate_format__blank_value_returns_false(
    validProtein,
):
    proteinList = [
        validProtein,
        "",
    ]
    assert not confirm_proteins_in_list_are_in_appropriate_format(proteinList)
