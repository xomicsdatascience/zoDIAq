import os.path
import pandas as pd
import numpy as np

def create_outfile_header(outputDir, queryFile, correction):
    outputCsodiaqTag = "CsoDIAq-file_"
    queryFileName = ".".join(queryFile.split("/")[-1].split(".")[:-1])
    outputFile = outputCsodiaqTag + queryFileName
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

    return f"{len(proteinList)}/{'/'.join(proteinList)}"
