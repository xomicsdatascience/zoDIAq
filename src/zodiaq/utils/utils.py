import os.path
import pandas as pd
import numpy as np
import re


def create_outfile_header(outputDir, queryFile, correction):
    outputZodiaqTag = "zoDIAq-file_"
    queryFileName = ".".join(queryFile.split("/")[-1].split(".")[:-1])
    outputFile = outputZodiaqTag + queryFileName
    outFileHeader = os.path.join(outputDir, outputFile)
    if correction != -1:
        outFileHeader += "_corrected"
    return outFileHeader


def format_protein_string_to_list(proteinString):
    """
    Returns a protein group string as a list of proteins.

    Parameters
    ----------
    proteinString : string
        Example: a protein group string would have the following format.
        '3/protein1/protein2/protein3'
        Where the number of proteins starts the string, followed by the proteins
        separated by slashes.

    Returns
    -------
    proteinList : list
        A list of protein names as strings, such as the following:
        [
            'protein1',
            'protein2',
            'protein3',
        ]
    """
    return proteinString.split("/")[1:]


def format_protein_list_to_string(proteinList):
    """
    Returns a list of proteins as a protein group string.

    Parameters
    ----------
    proteinList : list
        A list of protein names as strings, such as the following:
        [
            'protein1',
            'protein2',
            'protein3',
        ]

    Returns
    -------
        proteinString : string
        Example: a protein group string would have the following format.
        '3/protein1/protein2/protein3'
        Where the number of proteins starts the string, followed by the proteins
        separated by slashes.
    """

    return f"{len(proteinList)}/{'/'.join(sorted(proteinList))}"


def confirm_proteins_in_list_are_in_appropriate_format(proteinList):
    proteinPattern = re.compile(r"^\d+(?:\/.*[^\/])+$")
    for protein in proteinList:
        if not proteinPattern.search(protein) or not confirm_protein_count_is_accurate(
            protein
        ):
            return False
    return True


def confirm_protein_count_is_accurate(proteinStr):
    listFormat = proteinStr.split("/")
    count = listFormat[0]
    proteins = listFormat[1:]
    return len(proteins) == int(count)
