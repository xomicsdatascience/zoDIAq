import pytest
import argparse
import re
import os
from csodiaq import set_command_line_settings
from csodiaq.csodiaqParser import OutputDirectory, InputQueryFile, LibraryFile, RestrictedInt, RestrictedFloat, IdentificationOutputDirectory, ScoringOutputDirectory
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
def identificationOutputDirectory():
    return IdentificationOutputDirectory()

def test__csodiaq__identification_output_directory_parsing_class__fails_when_not_a_directory(identificationOutputDirectory):
    testFile = NamedTemporaryFile(prefix="csodiaq_test_file_", suffix=".txt")
    errorOutput = "The -i or --input argument must be a directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        identificationOutputDirectory(testFile.name)

def test__csodiaq__identification_output_directory_parsing_class__fails_when_directory_has_no_csodiaq_identification_outputs(identificationOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="fullOutput.csv", dir=testDir.name, delete=False)
    testFile3 = NamedTemporaryFile(prefix="CsoDIAq-file_3", suffix=".csv", dir=testDir.name, delete=False)
    testFile4 = NamedTemporaryFile(prefix="stuffCsoDIAq-file_4", suffix="fullOutput.csv", dir=testDir.name, delete=False)
    errorOutput = "The -i or --input argument directory must contain .csv files that are outputs from the identification workflow in CsoDIAq."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        identificationOutputDirectory(testDir.name)

@pytest.fixture
def scoringOutputDirectory():
    return ScoringOutputDirectory()

def test__csodiaq__scoring_output_directory_parsing_class__fails_when_not_a_directory(scoringOutputDirectory):
    testFile = NamedTemporaryFile(prefix="csodiaq_test_file_", suffix=".txt")
    errorOutput = "The -i or --input argument must be a directory."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        scoringOutputDirectory(testFile.name)

def test__csodiaq__scoring_output_directory_parsing_class__succeeds_when_directory_has_a_single_peptide_fdr_file(scoringOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="peptideFDR.csv", dir=testDir.name, delete=False)
    outputDict = scoringOutputDirectory(testDir.name)
    assert 'peptide' in outputDict
    assert len(outputDict['peptide']) == 1
    assert outputDict['peptide'][0] == testFile2.name.split('/')[-1]

def test__csodiaq__scoring_output_directory_parsing_class__succeeds_when_directory_has_multiple_peptide_fdr_files(scoringOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="peptideFDR.csv", dir=testDir.name, delete=False)
    testFile3 = NamedTemporaryFile(prefix="csodiaq_test_file3_", suffix="peptideFDR.csv", dir=testDir.name, delete=False)
    outputDict = scoringOutputDirectory(testDir.name)
    assert 'peptide' in outputDict
    assert len(outputDict['peptide']) == 2
    assert set(outputDict['peptide']) == set([testFile2.name.split('/')[-1], testFile3.name.split('/')[-1]])

def test__csodiaq__scoring_output_directory_parsing_class__succeeds_when_directory_has_a_single_protein_fdr_file(scoringOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="proteinFDR.csv", dir=testDir.name, delete=False)
    outputDict = scoringOutputDirectory(testDir.name)
    assert 'protein' in outputDict
    assert len(outputDict['protein']) == 1
    assert outputDict['protein'][0] == testFile2.name.split('/')[-1]

def test__csodiaq__scoring_output_directory_parsing_class__succeeds_when_directory_has_multiple_protein_fdr_files(scoringOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="proteinFDR.csv", dir=testDir.name, delete=False)
    testFile3 = NamedTemporaryFile(prefix="csodiaq_test_file3_", suffix="proteinFDR.csv", dir=testDir.name, delete=False)
    outputDict = scoringOutputDirectory(testDir.name)
    assert 'protein' in outputDict
    assert len(outputDict['protein']) == 2
    assert set(outputDict['protein']) == set([testFile2.name.split('/')[-1], testFile3.name.split('/')[-1]])

def test__csodiaq__scoring_output_directory_parsing_class__fails_when_directory_has_no_csodiaq_identification_outputs(scoringOutputDirectory):
    testDir = TemporaryDirectory(prefix="csodiaq_test_directory_")
    testFile1 = NamedTemporaryFile(prefix="csodiaq_test_file1_", suffix=".txt", dir=testDir.name, delete=False)
    testFile2 = NamedTemporaryFile(prefix="csodiaq_test_file2_", suffix="FDR.csv", dir=testDir.name, delete=False)
    errorOutput = "The -i or --input argument directory must contain .csv files that are outputs from the scoring workflow in CsoDIAq (peptide or protein score outputs required)."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        scoringOutputDirectory(testDir.name)

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
def scoreFiles():
    class fileObj:
        def __init__(self):
            self.idOutputDir = TemporaryDirectory(prefix="test_csodiaq_id_output_dir")
            self.idOutputFile1 = NamedTemporaryFile(prefix="CsoDIAq-file", suffix="fullOutput.csv", dir=self.idOutputDir.name, delete=False)
            self.idOutputFile2 = NamedTemporaryFile(prefix="CsoDIAq-file", suffix="fullOutput.csv", dir=self.idOutputDir.name, delete=False)
    return fileObj()

@pytest.fixture
def scoreArgs(scoreFiles):
    return [
        "score",
        "-i",
        scoreFiles.idOutputDir.name,
    ]

@pytest.fixture
def parsedScoreArgs(parser, scoreArgs):
    return vars(parser.parse_args(scoreArgs))

def test__csodiaq__set_command_line_settings__initialize_scoring(scoreFiles, parsedScoreArgs):
    assert parsedScoreArgs["command"] == "score"
    assert isinstance(parsedScoreArgs["input"], dict)
    assert parsedScoreArgs["input"]["csodiaqDirectory"] == scoreFiles.idOutputDir.name
    assert len(parsedScoreArgs["input"]["idFiles"]) == 2
    assert parsedScoreArgs["score"] == "macc"

@pytest.mark.skip('we may want to make an fdr scoring metric based purely on cosine. Keeping this here in case that becomes used.')
def test__csodiaq__set_command_line_settings__score_succeeds_with_cosine_input(parser, scoreArgs):
    scoreArgs += ['-s', 'cosine']
    args = vars(parser.parse_args(scoreArgs))
    assert args["score"] == "cosine"

def test__csodiaq__set_command_line_settings__score_succeeds_with_macc_input(parser, scoreArgs):
    scoreArgs += ['-s', 'macc']
    args = vars(parser.parse_args(scoreArgs))
    assert args["score"] == "macc"

def test__csodiaq__set_command_line_settings__score_fails_with_other_input(parser, scoreArgs):
    scoreArgs += ['-s', 'shouldFail']
    errorOutput = "argument -s/--score: invalid choice: 'shouldFail' (choose from 'macc')"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(scoreArgs))

@pytest.fixture
def reanalysisFiles():
    class fileObj:
        def __init__(self):
            self.scoreOutputDir = TemporaryDirectory(prefix="test_csodiaq_score_output_dir")
            self.peptideFdrFile1 = NamedTemporaryFile(suffix="peptideFDR.csv", dir=self.scoreOutputDir.name, delete=False)
            self.peptideFdrFile2 = NamedTemporaryFile(suffix="peptideFDR.csv", dir=self.scoreOutputDir.name, delete=False)
            self.proteinFdrFile1 = NamedTemporaryFile(suffix="proteinFDR.csv", dir=self.scoreOutputDir.name, delete=False)
            self.proteinFdrFile2 = NamedTemporaryFile(suffix="proteinFDR.csv", dir=self.scoreOutputDir.name, delete=False)
    return fileObj()


@pytest.fixture
def reanalysisArgs(reanalysisFiles):
    return [
        "targetedReanalysis",
        "-i",
        reanalysisFiles.scoreOutputDir.name,
    ]

@pytest.fixture
def parsedReanalysisArgs(parser, reanalysisArgs):
    return vars(parser.parse_args(reanalysisArgs))

def test__csodiaq__set_command_line_settings__initialize_targeted_reanalysis(reanalysisFiles, parsedReanalysisArgs):
    assert parsedReanalysisArgs["command"] == "targetedReanalysis"
    assert parsedReanalysisArgs["input"]["csodiaqDirectory"] == reanalysisFiles.scoreOutputDir.name
    assert isinstance(parsedReanalysisArgs["input"], dict)
    assert set(parsedReanalysisArgs["input"]["peptide"]) == set([reanalysisFiles.peptideFdrFile1.name.split('/')[-1], reanalysisFiles.peptideFdrFile2.name.split('/')[-1]])
    assert set(parsedReanalysisArgs["input"]["protein"]) == set([reanalysisFiles.proteinFdrFile1.name.split('/')[-1], reanalysisFiles.proteinFdrFile2.name.split('/')[-1]])
    assert parsedReanalysisArgs["protein"] == 0
    assert not parsedReanalysisArgs["heavyIsotope"]

def test__csodiaq__set_command_line_settings__targeted_reanalysis_fails_when_below_1(parser, reanalysisArgs):
    reanalysisArgs += ['-p', '0']
    errorOutput = "argument -p/--protein: The protein argument must be an integer greater than or equal to 1."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(reanalysisArgs))

def test__csodiaq__set_command_line_settings__targeted_reanalysis_fails_when_float_provided(parser, reanalysisArgs):
    reanalysisArgs += ['-p', '0.2']
    errorOutput = "argument -p/--protein: The protein argument must be an integer."
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(reanalysisArgs))

def test__csodiaq__set_command_line_settings__targeted_reanalysis_succeeds_with_custom_protein_input(parser, reanalysisArgs):
    reanalysisArgs += ['-p', '1']
    args = vars(parser.parse_args(reanalysisArgs))
    assert args['protein'] == 1

def test__csodiaq__set_command_line_settings__targeted_reanalysis_succeeds_with_heavy_isotope_flag(parser, reanalysisArgs):
    reanalysisArgs += ['-heavy']
    args = vars(parser.parse_args(reanalysisArgs))
    assert args['heavyIsotope']

def test__csodiaq__set_command_line_settings__targeted_reanalysis_fails_with_heavy_isotope_argument_added(parser, reanalysisArgs):
    reanalysisArgs += ['-heavy', 'randomBadValue']
    errorOutput = "unrecognized arguments: randomBadValue"
    with pytest.raises(argparse.ArgumentTypeError, match=re.escape(errorOutput)):
        args = vars(parser.parse_args(reanalysisArgs))

#TODO: make test that checks if the noCorrection tag, histogram tag and correctionDegree values are in conflict
#TODO: make test that checks if a protein number is provided to targetedReanalysis then proteinFDRs should be present
#TODO: evaluate if the -peaks tag should be included, or if we should move past that (currently just used to count the number peaks that contribute to the ionCount for the output)