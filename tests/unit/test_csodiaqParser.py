import pytest
import argparse
import re
import os
from csodiaq import set_command_line_settings
from csodiaq.csodiaqParser import OutputDirectory, InputQueryFile, LibraryFile, RestrictedInt, RestrictedFloat
from unittest.mock import Mock
from tempfile import TemporaryDirectory, NamedTemporaryFile

@pytest.fixture
def parser():
    def mock_error(message):
        raise argparse.ArgumentTypeError(message)
    argparse.ArgumentParser.error = Mock(side_effect=mock_error)
    parser = set_command_line_settings()
    return parser

def test__csodiaq__set_command_line_settings__gui__explicitly_flagged_succeeds(parser):
    args = vars(parser.parse_args(["gui"]))
    assert args["command"] == "gui"

def test__csodiaq__set_command_line_settings__gui__succeeds_with_no_command_provided(parser):
    args = vars(parser.parse_args([]))
    assert not args["command"]

@pytest.fixture
def testOutputDir():
    return OutputDirectory("test")

def test__csodiaq__output_directory_parsing_class__rejects_file_as_input(testOutputDir):
    testFile = NamedTemporaryFile(prefix="csodiaq_test_file_", suffix=".txt")
    errorOutput = "The -o or --output argument must be a directory or to-be-created directory header, not an existing file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        testOutputDir(testFile.name)

def test__csodiaq__output_directory_parsing_class__rejects_input_that_has_no_existing_parent_directory(testOutputDir):
    testPath = 'this/path/does/not/exist'
    errorOutput = "The -o or --output argument directory requires an existing parent directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        testOutputDir(testPath)

@pytest.fixture
def outputDirectoryPattern():
    return r"\d{8}-\d{6}"

def test__csodiaq__output_directory_parsing_class__if_directory_exists_new_directory_made_in_directory(testOutputDir, outputDirectoryPattern):
    testDirectory = TemporaryDirectory(prefix="csodiaq_input_test_directory_")
    errorOutput = "The -o or --output argument directory requires an existing parent directory."
    newDirectoryPath = testOutputDir(testDirectory.name)
    newDirectoryPathName = newDirectoryPath.split('/')[-1]
    expectedDirectoryNamePattern = re.compile(rf'csodiaq-test-{outputDirectoryPattern}')
    assert expectedDirectoryNamePattern.search(newDirectoryPathName)

def test__csodiaq__output_directory_parsing_class__if_directory_does_not_exist_new_directory_made_with_header_at_end_of_given_path(testOutputDir, outputDirectoryPattern):
    testDirectory = TemporaryDirectory(prefix="csodiaq_input_test_directory_")
    newDirectoryHeader = 'newDirectoryHeader'
    inputPath = os.path.join(testDirectory.name, newDirectoryHeader)
    errorOutput = "The -o or --output argument directory requires an existing parent directory."
    newDirectoryPath = testOutputDir(inputPath)
    newDirectoryPathName = newDirectoryPath.split('/')[-1]
    expectedDirectoryNamePattern = re.compile(rf'{newDirectoryHeader}-csodiaq-test-{outputDirectoryPattern}')
    assert expectedDirectoryNamePattern.search(newDirectoryPathName)

@pytest.fixture
def testInputFile():
    return InputQueryFile()

def test__csodiaq__input_query_file_parsing_class__rejects_non_existing_file_values(testInputFile):
    testFile = 'this/path/does/not/exist'
    errorOutput = "The -i or --input argument must be an existing file (and not a directory)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        testInputFile(testFile)

def test__csodiaq__input_query_file_parsing_class__rejects_files_of_wrong_type(testInputFile):
    testFile = NamedTemporaryFile(prefix='csodiaq_input_test_file_', suffix=".txt")
    errorOutput = "The -i or --input argument must be an .mzXML file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        testInputFile(testFile.name)

@pytest.fixture
def libraryFile():
    return LibraryFile()

def test__csodiaq__library_file_parsing_class__rejects_non_existing_file_values(libraryFile):
    testFile = 'this/path/does/not/exist'
    errorOutput = "The -l or --library argument must be an existing file (and not a directory)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        libraryFile(testFile)

def test__csodiaq__library_file_parsing_class__rejects_files_of_wrong_type(libraryFile):
    testFile = NamedTemporaryFile(prefix='csodiaq_input_test_file_', suffix=".txt")
    errorOutput = "The -l or --library argument must be a .tsv, .csv or .mgf file."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        libraryFile(testFile.name)

@pytest.fixture
def restrictedInt():
    return RestrictedInt("test", minValue=-2, maxValue=3)

def test__csodiaq__restricted_int_parsing_class__fails_when_value_cannot_be_assigned_to_int(restrictedInt):
    inputValue = "string"
    errorOutput = "The test argument must be an integer."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedInt(inputValue)

def test__csodiaq__restricted_int_parsing_class__fails_when_value_is_below_min_value(restrictedInt):
    inputValue = "-3"
    errorOutput = "The test argument must be an integer greater than or equal to -2."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedInt(inputValue)

def test__csodiaq__restricted_int_parsing_class__fails_when_value_is_above_max_value(restrictedInt):
    inputValue = "4"
    errorOutput = "The test argument must be an integer less than or equal to 3."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedInt(inputValue)

def test__csodiaq__restricted_int_parsing_class__succeeds_when_value_is_equal_to_min_value(restrictedInt):
    inputValue = "-2"
    output = restrictedInt(inputValue)
    assert isinstance(output, int)
    assert output == -2

def test__csodiaq__restricted_int_parsing_class__succeeds_when_value_is_greater_than_min_value_less_than_max_value(restrictedInt):
    inputValue = "0"
    output = restrictedInt(inputValue)
    assert isinstance(output, int)
    assert output == 0

def test__csodiaq__restricted_int_parsing_class__succeeds_when_value_is_equal_to_max_value(restrictedInt):
    inputValue = "3"
    output = restrictedInt(inputValue)
    assert isinstance(output, int)
    assert output == 3

@pytest.fixture
def restrictedFloat():
    return RestrictedFloat("test", minValue=-2, maxValue=3)

def test__csodiaq__restricted_float_parsing_class__fails_when_value_cannot_be_assigned_to_float(restrictedFloat):
    inputValue = "string"
    errorOutput = "The test argument must be a float."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedFloat(inputValue)

def test__csodiaq__restricted_float_parsing_class__fails_when_value_is_below_min_value(restrictedFloat):
    inputValue = "-3"
    errorOutput = "The test argument must be a float greater than or equal to -2."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedFloat(inputValue)

def test__csodiaq__restricted_float_parsing_class__fails_when_value_is_above_max_value(restrictedFloat):
    inputValue = "4"
    errorOutput = "The test argument must be a float less than or equal to 3."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        restrictedFloat(inputValue)

def test__csodiaq__restricted_float_parsing_class__succeeds_when_value_is_equal_to_min_value(restrictedFloat):
    inputValue = "-2"
    output = restrictedFloat(inputValue)
    assert isinstance(output, float)
    assert output == -2

def test__csodiaq__restricted_float_parsing_class__succeeds_when_value_is_greater_than_min_value_less_than_max_value(restrictedFloat):
    inputValue = "0"
    output = restrictedFloat(inputValue)
    assert isinstance(output, float)
    assert output == 0

def test__csodiaq__restricted_float_parsing_class__succeeds_when_value_is_equal_to_max_value(restrictedFloat):
    inputValue = "3"
    output = restrictedFloat(inputValue)
    assert isinstance(output, float)
    assert output == 3

@pytest.fixture
def idFiles():
    class fileObj:
        def __init__(self):
            self.parentDir = TemporaryDirectory(prefix="csodiaq_input_test_directory_")
            self.outputDir = TemporaryDirectory(
                prefix="csodiaq_output_test_directory_id_", dir=self.parentDir.name
            )
            self.inputFile = NamedTemporaryFile(prefix="csodiaq_input_test_file_id_", suffix=".mzXML")
            self.tramlCsvLibraryFile = NamedTemporaryFile(prefix="csodiaq_traml_library_file_csv_id_", suffix=".csv")
            self.tramlTsvLibraryFile = NamedTemporaryFile(prefix="csodiaq_traml_library_file_tsv_id_", suffix=".tsv")
            self.mgfLibraryFile = NamedTemporaryFile(prefix="csodiaq_mgf_library_file_id_", suffix=".mgf")
    return fileObj()


@pytest.fixture
def idArgs(idFiles):
    return [
        "id",
        "-o",
        idFiles.outputDir.name,
        "-i",
        idFiles.inputFile.name,
        "-l",
        idFiles.tramlCsvLibraryFile.name,
    ]

def test__csodiaq__set_command_line_settings__initialize_identification(parser, idArgs):
    parsedIdArgs = vars(parser.parse_args(idArgs))
    assert parsedIdArgs["command"] == "id"
    assert parsedIdArgs["output"]
    assert os.path.isdir(parsedIdArgs["output"])
    assert parsedIdArgs["input"]
    assert parsedIdArgs["matchTolerance"] == 30
    assert isinstance(parsedIdArgs["matchTolerance"], float)
    assert not parsedIdArgs["noCorrection"]
    assert parsedIdArgs["correctionDegree"] == 0
    assert not parsedIdArgs["histogram"]
    assert parsedIdArgs["maxQueryPeaks"] == np.inf

def test__csodiaq__set_command_line_settings__id_fails_when_no_output_dir_provided(parser, idArgs):
    idArgsWithoutOutput = idArgs[:1] + idArgs[3:]
    errorOutput = 'the following arguments are required: -o/--output'

    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgsWithoutOutput))

def test__csodiaq__set_command_line_settings__id_fails_when_no_input_file_provided(parser, idArgs):
    idArgsWithoutInput = idArgs[:3] + idArgs[5:]
    errorOutput = 'the following arguments are required: -i/--input'
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgsWithoutInput))

def test__csodiaq__set_command_line_settings__id_succeeds_with_multiple_input_files(parser, idArgs):
    idArgsWithMultipleInputs = idArgs + idArgs[3:5]
    args = vars(parser.parse_args(idArgsWithMultipleInputs))
    assert args["input"]
    assert len(args["input"]) == 2
    assert args["input"][0] == args["input"][1]

def test__csodiaq__set_command_line_settings__id_fails_when_no_library_file_provided(parser, idArgs):
    idArgsWithoutLibrary = idArgs[:5] + idArgs[7:]
    errorOutput = 'the following arguments are required: -l/--library'
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgsWithoutLibrary))

def test__csodiaq__set_command_line_settings__id_succeeds_with_custom_match_tolerance(parser, idArgs):
    idArgs += ['-t', '10']
    args = vars(parser.parse_args(idArgs))
    assert args['matchTolerance'] == 10

def test__csodiaq__set_command_line_settings__id_fails_when_match_tolerance_less_than_1(parser, idArgs):
    idArgs += ['-t', '0']
    errorOutput = "The matchTolerance argument must be a float greater than or equal to 1."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__id_fails_when_match_tolerance_greater_than_100(parser, idArgs):
    idArgs += ['-t', '61']
    errorOutput = "The matchTolerance argument must be a float less than or equal to 60."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__id_succeeds_with_custom_match_tolerance(parser, idArgs):
    idArgs += ['-t', '10']
    args = vars(parser.parse_args(idArgs))
    assert args['matchTolerance'] == 10

def test__csodiaq__set_command_line_settings__id_succeeds_with_no_correction_flag(parser, idArgs):
    idArgs += ['-nc']
    args = vars(parser.parse_args(idArgs))
    assert args['noCorrection']

def test__csodiaq__set_command_line_settings__id_fails_with_no_correction_arguments_added(parser, idArgs):
    idArgs += ['-nc', 'randomBadValue']
    errorOutput = "unrecognized arguments: randomBadValue"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__id_fails_when_correction_degree_less_than_p5(parser, idArgs):
    idArgs += ['-c', '0.4']
    errorOutput = "The correctionDegree argument must be a float greater than or equal to 0.5."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__id_fails_when_correction_degree_greater_than_2(parser, idArgs):
    idArgs += ['-c', '2.1']
    errorOutput = "The correctionDegree argument must be a float less than or equal to 2."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

def test__csodiaq__set_command_line_settings__id_succeeds_with_custom_correction_degree(parser, idArgs):
    idArgs += ['-c', '1']
    args = vars(parser.parse_args(idArgs))
    assert args['correctionDegree'] == 1

def test__csodiaq__set_command_line_settings__id_succeeds_with_histogram_flag(parser, idArgs):
    idArgs += ['-hist']
    args = vars(parser.parse_args(idArgs))
    assert args['histogram']

def test__csodiaq__set_command_line_settings__id_fails_with_histogram_argument_added(parser, idArgs):
    idArgs += ['-hist', 'randomBadValue']
    errorOutput = "unrecognized arguments: randomBadValue"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(idArgs))

@pytest.fixture
def scoreArgs():
    return ["score"]

@pytest.fixture
def parsedScoreArgs(parser, scoreArgs):
    return vars(parser.parse_args(scoreArgs))

def test__csodiaq__set_command_line_settings__initialize_scoring(parsedScoreArgs):
    assert parsedScoreArgs["command"] == "score"

@pytest.fixture
def reanalysisArgs():
    return ["targetedReanalysis"]

@pytest.fixture
def parsedReanalysisArgs(parser, reanalysisArgs):
    return vars(parser.parse_args(reanalysisArgs))

def test__csodiaq__set_command_line_settings__initialize_targeted_reanalysis(parsedReanalysisArgs):
    assert parsedReanalysisArgs["command"] == "targetedReanalysis"

#TODO: make test that checks if the noCorrection tag, histogram tag and correctionDegree values are in conflict
#TODO: evaluate if the -peaks tag should be included, or if we should move past that (currently just used to count the number peaks that contribute to the ionCount for the output)