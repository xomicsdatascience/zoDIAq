import argparse
import os
import time
import numpy as np
from abc import ABC, abstractmethod
import re
import warnings


def set_args_from_command_line_input():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="zoDIAq Functions")
    guiParser = commandParser.add_parser(
        "gui", help="Launches the (optional) GUI application for using zoDIAq."
    )
    add_id_parser(commandParser)
    add_score_parser(commandParser)
    add_reanalysis_parser(commandParser)
    return parser


def add_id_parser(commandParser):
    idParser = commandParser.add_parser(
        "id",
        help="Identify peptides from a designated peptide library in mass spectrometry run (query) files.",
    )
    idParser.add_argument(
        "-o",
        "--output",
        type=_OutputDirectory("id"),
        required=True,
        help="Output directory to write output files to. A new directory will be created in this path.\nRequired.",
    )
    idParser.add_argument(
        "-i",
        "--input",
        type=_InputQueryFile(),
        required=True,
        action="append",
        help="mzXML input files from processed mass spectrometry .RAW files.\nRequired.",
    )
    idParser.add_argument(
        "-l",
        "--library",
        type=_LibraryFile(),
        required=True,
        help="File containing spectra of peptides to identify in query files.\nRequired.\nAccepts TraML (.csv or .tsv) and MGF (.mgf) formats.",
    )
    idParser.add_argument(
        "-t",
        "--matchTolerance",
        type=_RestrictedFloat("matchTolerance", minValue=1, maxValue=60),
        default=30.0,
        help="Tolerance between library and query peak m/z values to be considered a match before correction (in PPM).\nOptional.\nDefault value is 30.\nValue must be greater than 0.\nValue must be less than 60 as a generic cutoff to prevent needlessly taxing computational power.",
    )
    idParser.add_argument(
        "-nc",
        "--noCorrection",
        default=False,
        action="store_true",
        help="Disables correction of ppm tolerance when matching library and query peaks. \nOptional. This option is NOT recommended.",
    )
    idParser.add_argument(
        "-c",
        "--correctionDegree",
        type=_RestrictedFloat("correctionDegree", minValue=0.5, maxValue=2),
        default=0,
        help="Uses mean and standard deviation to correct ppm tolerance.\nOptional. A customized 'binning' correction method described in the paper is used by default.",
    )
    idParser.add_argument(
        "-hist",
        "--histogram",
        default=False,
        action="store_true",
        help="This flag indicates a histogram of the uncorrected PPM values (with lines for the chosen offset/tolerance) should be generated.\nOptional.",
    )
    idParser.add_argument(
        "-w",
        "--cancelWarnings",
        default=False,
        action="store_true",
        help="This flag indicates that warning errors should be oppressed.\nOptional.",
    )


def add_score_parser(commandParser):
    scoringParser = commandParser.add_parser(
        "score",
        help="Scores confidence in identified peptides and removes those above a False Discovery Rate (FDR) of 0.01.\nNOTE: The identification output must include decoys for accurate FDR scoring.",
    )
    scoringParser.add_argument(
        "-i",
        "--input",
        type=_IdentificationOutputDirectory(),
        required=True,
        help="Directory that contains outputs of the identification step of zoDIAq.\nRequired.",
    )
    scoringParser.add_argument(
        "-s",
        "--score",
        choices=["macc"],
        default="macc",
        help="Type of scoring to apply for calculating the False Discovery Rate (FDR). Default is MaCC score (see paper for explanation). \nOptional, default is 'macc' method. 'macc' is the only current choice.",
    )
    scoringParser.add_argument(
        "-p",
        "--proteinQuantMethod",
        choices=["maxlfq", "sum"],
        default="maxlfq",
        help="Method by which protein quantification metric is calculated (based on peptide quantities). \nOptional, default is 'maxlfq' method. Choices are 'maxlfq' or 'sum'.",
    )
    scoringParser.add_argument(
        "-min",
        "--minNumDifferences",
        type=_RestrictedInt("minNumDifferences", minValue=1, maxValue=2),
        default=2,
        help="Specific to the maxLFQ protein quantification method. Requires at minimum the given number of matches before a sample to sample ratio or difference is accepted.\nOptional, default is 2. Only 1 or 2 is accepted. \nThis flag will throw a warning error when paired with a protein quantification method other than 'maxlfq'.",
    )


def add_reanalysis_parser(commandParser):
    reanalysisParser = commandParser.add_parser(
        "targetedReanalysis",
        help="Creates files readable by a mass spectrometer for targeted reanalysis of identified peptides.",
    )
    reanalysisParser.add_argument(
        "-i",
        "--input",
        type=_ScoringOutputDirectory(),
        required=True,
        help="Directory that contains outputs of the scoring step of zoDIAq.\nRequired.",
    )
    reanalysisParser.add_argument(
        "-p",
        "--protein",
        type=_RestrictedInt("protein", minValue=1),
        default=0,
        help="Determines the maximum number of peptides per identified protein to include.\nOptional. Not setting this variable will result in evaluating peptides only with no reference to proteins.",
    )
    reanalysisParser.add_argument(
        "-heavy",
        "--heavyIsotope",
        default=False,
        action="store_true",
        help="This flag indicates that files for targeted re-analysis should include heavy fragment isotopes for SILAC quantification.\nOptional.",
    )
    reanalysisParser.add_argument(
        "-b",
        "--binValueProximity",
        type=_RestrictedFloat("binValueProximity", minValue=0.01),
        default=0.75,
        help="When setting bin values, this option indicates how close an m/z value must be to the bin value. Default is 0.75.\nOptional.\nNOTE: Multiple targeted m/z values may fall within a range that a mass spectrometer can identify in one scan. Thus, m/z values are binned to prevent redundant reanalysis.\nExample: let's say we have the m/z values of 199.5 and 200.5, a binValueProximity value of 0.75, and a bin value of 200.0.\nBoth of these m/z values would be in the same bin, as they are both with 0.75 of 200.0.",
    )


def check_for_conflicting_args(args):
    if args["command"] == "id" and args["histogram"] and args["noCorrection"]:
        raise argparse.ArgumentTypeError(
            "The histogram flag is invalidated by the noCorrection flag. Please inspect your input and remove one of the tags."
        )
    if args["command"] == "id" and args["correctionDegree"] and args["noCorrection"]:
        raise argparse.ArgumentTypeError(
            "The correctionDegree parameter is invalidated by the noCorrection flag. Please inspect your input and remove one of them."
        )
    if (
        args["command"] == "score"
        and args["proteinQuantMethod"] != "maxlfq"
        and args["minNumDifferences"] != 2
    ):
        warnings.warn(
            f"The minNumDifferences flag will only have an effect when paired with the 'maxlfq' proteinQuantMethod flag. You used it with the '{args['proteinQuantMethod']}' method, so this flag will be ignored.",
            UserWarning,
        )
    if (
        args["command"] == "targetedReanalysis"
        and args["protein"]
        and len(args["input"]["protein"]) == 0
    ):
        raise argparse.ArgumentTypeError(
            "The protein argument requires the presence of protein FDR files to function. Please run the protein scoring workflow or remove the protein argument from your commands."
        )


def get_output_name(commandName):
    return f"zodiaq-{commandName}-{time.strftime('%Y%m%d-%H%M%S')}"


def get_new_output_directory_path(newDirectoryLocation, commandName):
    newDirectoryName = get_output_name(commandName)
    if os.path.isdir(newDirectoryLocation):
        newDirectoryPath = os.path.join(newDirectoryLocation, newDirectoryName)
    else:
        newDirectoryHeader = os.path.basename(newDirectoryLocation)
        newDirectoryParentDirectory = os.path.dirname(newDirectoryLocation)
        newDirectoryPath = os.path.join(
            newDirectoryParentDirectory, f"{newDirectoryHeader}-{newDirectoryName}"
        )
    return newDirectoryPath


class _OutputDirectory:
    def __init__(self, commandName):
        self.commandName = commandName

    def __call__(self, newDirectoryLocation):
        if os.path.isfile(newDirectoryLocation):
            raise argparse.ArgumentTypeError(
                "The -o or --output argument must be a directory or to-be-created directory header, not an existing file."
            )
        if not os.path.isdir(os.path.dirname(newDirectoryLocation)):
            raise argparse.ArgumentTypeError(
                "The -o or --output argument directory requires an existing parent directory."
            )
        newDirectoryPath = get_new_output_directory_path(
            newDirectoryLocation, self.commandName
        )
        return newDirectoryPath


class _InputQueryFile:
    def __init__(self):
        self.allowedFileTypes = [".mzxml"]

    def __call__(self, inputQueryFile):
        if not os.path.isfile(inputQueryFile):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be an existing file (and not a directory)."
            )
        if os.path.splitext(inputQueryFile)[1].lower() not in self.allowedFileTypes:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be an .mzXML file."
            )
        return inputQueryFile


class _LibraryFile:
    def __init__(self):
        self.allowedFileTypes = ["csv", "tsv", "mgf"]

    def __call__(self, libraryFile):
        if not os.path.isfile(libraryFile):
            raise argparse.ArgumentTypeError(
                "The -l or --library argument must be an existing file (and not a directory)."
            )
        if libraryFile.split(".")[-1] not in self.allowedFileTypes:
            raise argparse.ArgumentTypeError(
                "The -l or --library argument must be a .tsv, .csv or .mgf file."
            )
        return libraryFile


class _RestrictedNumber(ABC):
    def __init__(self, type, minValue=-np.inf, maxValue=np.inf):
        self.type = type
        self.minValue = minValue
        self.maxValue = maxValue
        assert self.minValue <= self.maxValue

    @abstractmethod
    def return_expected_type(self):
        pass

    @abstractmethod
    def coerce_into_expected_type(self):
        pass

    def __call__(self, value):
        try:
            value = self.coerce_into_expected_type(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be {self.return_expected_type()}."
            )
        if value < self.minValue:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be {self.return_expected_type()} greater than or equal to {self.minValue}."
            )
        if value > self.maxValue:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be {self.return_expected_type()} less than or equal to {self.maxValue}."
            )
        return value


class _RestrictedInt(_RestrictedNumber):
    def return_expected_type(self):
        return "an integer"

    def coerce_into_expected_type(self, value):
        return int(value)


class _RestrictedFloat(_RestrictedNumber):
    def return_expected_type(self):
        return "a float"

    def coerce_into_expected_type(self, value):
        return float(value)

    def __call__(self, value):
        value = super().__call__(value)
        if self.type == "binValueProximity" and not round(value, 2) == value:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument cannot have values beyond 2 decimal places (mass spectrometers are typically not sensitive enough for that specificity)."
            )
        return value


class _ZodiaqOutputDirectory(ABC):
    def __call__(self, idDir):
        if not os.path.isdir(idDir):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be a directory."
            )
        outputDict = self.add_necessary_directory_contents(idDir)
        outputDict["zodiaqDirectory"] = idDir
        return outputDict

    def find_files_with_necessary_format(self, dir, regexPattern):
        filePattern = re.compile(regexPattern)
        dirContentList = os.listdir(dir)
        files = [file for file in dirContentList if filePattern.search(file)]
        return files

    @abstractmethod
    def add_necessary_directory_contents(self, directory):
        pass


class _IdentificationOutputDirectory(_ZodiaqOutputDirectory):
    def add_necessary_directory_contents(self, idDir):
        idFiles = self.find_files_with_necessary_format(
            idDir, r"^zoDIAq-file.*fullOutput\.csv"
        )
        if len(idFiles) == 0:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument directory must contain .csv files that are outputs from the identification workflow in zoDIAq."
            )
        return {"idFiles": idFiles}


class _ScoringOutputDirectory(_ZodiaqOutputDirectory):
    def add_necessary_directory_contents(self, scoreDir):
        peptideFdrFiles = self.find_files_with_necessary_format(
            scoreDir, r"peptideFDR\.csv$"
        )
        proteinFdrFiles = self.find_files_with_necessary_format(
            scoreDir, r"proteinFDR\.csv$"
        )
        if len(peptideFdrFiles) == 0 and len(proteinFdrFiles) == 0:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument directory must contain .csv files that are outputs from the scoring workflow in zoDIAq (peptide or protein score outputs required)."
            )

        return {
            "peptide": peptideFdrFiles,
            "protein": proteinFdrFiles,
        }
