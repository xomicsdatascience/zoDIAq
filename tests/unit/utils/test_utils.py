from csodiaq.utils import (
    create_outfile_header,
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
def outputCsodiaqTag():
    return "CsoDIAq-file_"


def test__output_writing_functions__create_outfile_header__no_correction(
    outputDirectory, inputFileName, inputFilePath, outputCsodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputCsodiaqTag}{inputFileName}"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=-1)
    assert expectedOutput == output


def test__output_writing_functions__create_outfile_header__no_correction__output_directory_ends_in_slash(
    outputDirectory, inputFileName, inputFilePath, outputCsodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputCsodiaqTag}{inputFileName}"
    output = create_outfile_header(outputDirectory + "/", inputFilePath, correction=-1)
    assert expectedOutput == output


def test__output_writing_functions__create_outfile_header__no_correction__includes_non_file_type_dots(
    outputDirectory, inputFileName, inputFilePath, outputCsodiaqTag
):
    inputFileNameWithPeriods = inputFileName + ".dots.added"
    inputFilePathWithPeriods = inputFileNameWithPeriods + ".mzxml"
    expectedOutput = f"{outputDirectory}/{outputCsodiaqTag}{inputFileNameWithPeriods}"
    output = create_outfile_header(
        outputDirectory, inputFilePathWithPeriods, correction=-1
    )
    assert expectedOutput == output


def test__output_writing_functions__create_outfile_header__custom_correction(
    outputDirectory, inputFileName, inputFilePath, outputCsodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputCsodiaqTag}{inputFileName}_corrected"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=0)
    assert expectedOutput == output


def test__output_writing_functions__create_outfile_header__stdev_correction(
    outputDirectory, inputFileName, inputFilePath, outputCsodiaqTag
):
    expectedOutput = f"{outputDirectory}/{outputCsodiaqTag}{inputFileName}_corrected"
    output = create_outfile_header(outputDirectory, inputFilePath, correction=1)
    assert expectedOutput == output
