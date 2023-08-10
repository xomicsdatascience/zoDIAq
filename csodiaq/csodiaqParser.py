import argparse
import os
import time
import numpy as np
from abc import ABC, abstractmethod
import re

def set_command_line_settings():
    parser = argparse.ArgumentParser(description="")
    commandParser = parser.add_subparsers(dest="command", help="CsoDIAq Functions")
    guiParser = commandParser.add_parser(
        'gui', help='Launches the (optional) GUI application for using CsoDIAq.')
    add_id_parser(commandParser)
    add_score_parser(commandParser)
    add_reanalysis_parser(commandParser)
    return parser

def add_id_parser(commandParser):
    idParser = commandParser.add_parser('id', help='Identify peptides from a designated peptide library in mass spectrometry run (query) files.')
    idParser.add_argument(
        "-o",
        "--output",
        type=OutputDirectory("id"),
        required=True,
        help="Output directory to write output files to. A new directory will be created in this path.\nRequired.",
    )
    idParser.add_argument(
        "-i",
        "--input",
        type=InputQueryFile(),
        required=True,
        action='append',
        help="mzXML input files from processed mass spectrometry .RAW files.\nRequired."
    )
    idParser.add_argument(
        "-l",
        "--library",
        type=LibraryFile(),
        required=True,
        help="File containing spectra of peptides to identify in query files.\nRequired.\nAccepts TraML (.csv or .tsv) and MGF (.mgf) formats."
    )
    idParser.add_argument(
        "-t",
        "--matchTolerance",
        type=RestrictedFloat("matchTolerance", minValue=1, maxValue=60),
        default=30.0,
        help="Tolerance between library and query peak m/z values to be considered a match before correction (in PPM).\nOptional.\nDefault value is 30.\nValue must be greater than 0.\nValue must be less than 60 as a generic cutoff to prevent needlessly taxing computational power."
    )
    idParser.add_argument(
        "-nc",
        "--noCorrection",
        default=False,
        action='store_true',
        help="Disables correction of ppm tolerance when matching library and query peaks. \nOptional. This option is NOT recommended."
    )
    idParser.add_argument(
        "-c",
        "--correctionDegree",
        type=RestrictedFloat("correctionDegree", minValue=0.5, maxValue=2),
        default=0,
        help="Uses mean and standard deviation to correct ppm tolerance.\nOptional. A customized 'binning' correction method described in the paper is used by default."
    )
    idParser.add_argument(
        "-hist",
        "--histogram",
        default=False,
        action='store_true',
        help="This flag indicates a histogram of the uncorrected PPM values (with lines for the chosen offset/tolerance) should be generated.\nOptional."
    )

def add_score_parser(commandParser):
    scoringParser = commandParser.add_parser('score', help='Scores confidence in identified peptides and removes those above a False Discovery Rate (FDR) of 0.01.\nNOTE: The identification output must include decoys for accurate FDR scoring.')
    scoringParser.add_argument(
        "-i",
        "--input",
        type=IdentificationOutputDirectory(),
        required=True,
        help="Directory that contains outputs of the identification step of CsoDIAq.\nRequired."
    )
    scoringParser.add_argument(
        "-s",
        "--score",
        choices=["macc"],
        default="macc",
        help="Type of scoring to apply for calculating the False Discovery Rate (FDR). Default is MaCC score (see paper for explanation). \nOptional. Choices are cosine and MaCC."
    )




def add_reanalysis_parser(commandParser):
    reanalysisParser = commandParser.add_parser('targetedReanalysis', help='Creates files readable by a mass spectrometer for targeted reanalysis of identified peptides.')
    reanalysisParser.add_argument(
        "-i",
        "--input",
        type=ScoringOutputDirectory(),
        required=True,
        help="Directory that contains outputs of the scoring step of CsoDIAq.\nRequired."
    )
    reanalysisParser.add_argument(
        "-p",
        "--protein",
        type=RestrictedInt("protein", minValue=1),
        default=0,
        help="Determines the maximum number of peptides per identified protein to include.\nOptional. Not setting this variable will result in "
    )
    reanalysisParser.add_argument(
        "-heavy",
        "--heavyIsotope",
        default=False,
        action='store_true',
        help="This flag indicates that files for targeted re-analysis should include heavy fragment isotopes for SILAC quantification.\nOptional."
    )



def create_output_name(commandName):
    return f"csodiaq-{commandName}-{time.strftime('%Y%m%d-%H%M%S')}"

def create_new_output_directory_path(newDirectoryLocation, commandName):
    if newDirectoryLocation[-1] == "/":
        newDirectoryLocation = newDirectoryLocation[:-1]
    newDirectoryName = create_output_name(commandName)
    if os.path.isdir(newDirectoryLocation):
        newDirectoryPath = os.path.join(newDirectoryLocation, newDirectoryName)
    else:
        newDirectoryHeader = newDirectoryLocation.split('/')[-1]
        newDirectoryParentDirectory = '/'.join(newDirectoryLocation.split('/')[:-1])
        newDirectoryPath = os.path.join(newDirectoryParentDirectory, f'{newDirectoryHeader}-{newDirectoryName}')
    return newDirectoryPath

class OutputDirectory:
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
        newDirectoryPath = create_new_output_directory_path(newDirectoryLocation, self.commandName)
        os.mkdir(newDirectoryPath)
        return newDirectoryPath

class InputQueryFile:
    def __init__(self):
        self.allowedFileTypes = ["mzXML"]

    def __call__(self, inputQueryFile):
        if not os.path.isfile(inputQueryFile):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be an existing file (and not a directory)."
            )
        if inputQueryFile.split(".")[-1] not in self.allowedFileTypes:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be an .mzXML file."
            )
        return inputQueryFile

class LibraryFile:
    def __init__(self):
        self.allowedFileTypes = ["csv","tsv","mgf"]

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

class RestrictedInt:
    def __init__(self, type, minValue=-np.inf, maxValue=np.inf):
        self.type = type
        self.minValue = minValue
        self.maxValue = maxValue
        assert self.minValue <= self.maxValue

    def __call__(self, intValue):
        try:
            intValue = int(intValue)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be an integer."
            )
        if intValue < self.minValue:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be an integer greater than or equal to {self.minValue}."
            )
        if intValue > self.maxValue:
            raise argparse.ArgumentTypeError(
                f"The {self.type} argument must be an integer less than or equal to {self.maxValue}."
            )
        return intValue

class RestrictedNumber(ABC):
    def __init__(self, type, minValue=-np.inf, maxValue=np.inf):
        self.type = type
        self.minValue = minValue
        self.maxValue = maxValue
        assert self.minValue <= self.maxValue

    @abstractmethod
    def return_expected_type(self): pass

    @abstractmethod
    def coerce_into_expected_type(self): pass

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

class RestrictedInt(RestrictedNumber):
    def return_expected_type(self):
        return 'an integer'

    def coerce_into_expected_type(self, value):
        return int(value)

class RestrictedFloat(RestrictedNumber):
    def return_expected_type(self):
        return 'a float'

    def coerce_into_expected_type(self, value):
        return float(value)

class CsodiaqOutputDirectory(ABC):
    def __call__(self, idDir):
        if not os.path.isdir(idDir):
            raise argparse.ArgumentTypeError(
                "The -i or --input argument must be a directory."
            )
        outputDict = self.add_necessary_directory_contents(idDir)
        outputDict['csodiaqDirectory'] = idDir
        return outputDict

    def find_files_with_necessary_format(self, dir, regexPattern):
        filePattern = re.compile(regexPattern)
        dirContentList = os.listdir(dir)
        files = [file for file in dirContentList if filePattern.search(file)]
        return files

    @abstractmethod
    def add_necessary_directory_contents(self, directory): pass

class IdentificationOutputDirectory(CsodiaqOutputDirectory):
    def add_necessary_directory_contents(self, idDir):
        idFiles = self.find_files_with_necessary_format(idDir, r'^CsoDIAq-file.*fullOutput\.csv')
        if len(idFiles) == 0:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument directory must contain .csv files that are outputs from the identification workflow in CsoDIAq."
            )
        return {'idFiles': idFiles}

class ScoringOutputDirectory(CsodiaqOutputDirectory):
    def add_necessary_directory_contents(self, scoreDir):

        peptideFdrFiles = self.find_files_with_necessary_format(scoreDir, r'peptideFDR\.csv$')
        proteinFdrFiles = self.find_files_with_necessary_format(scoreDir, r'proteinFDR\.csv$')
        if len(peptideFdrFiles) == 0 and len(proteinFdrFiles) == 0:
            raise argparse.ArgumentTypeError(
                "The -i or --input argument directory must contain .csv files that are outputs from the scoring workflow in CsoDIAq (peptide or protein score outputs required)."
            )

        return {
            'peptide': peptideFdrFiles,
            'protein': proteinFdrFiles,
        }
