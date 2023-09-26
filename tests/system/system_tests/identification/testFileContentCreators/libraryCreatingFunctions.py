import pandas as pd
import numpy as np
import os
from pyteomics import mgf


def get_parent_dir():
    return os.path.dirname(os.path.abspath(__file__))


def create_template_library_dataframe(
    minMz, maxMz, precursorMzDiff, peakMzDiff, peakIntensityDiff
):
    data = []
    for precursorMz in np.arange(minMz, maxMz, precursorMzDiff):
        if precursorMz < (maxMz + minMz) / 2:
            protein = f"1/protein{precursorMz}"
            id = f"id{precursorMz}"
        else:
            protein = f"1/DECOY_protein{precursorMz}"
            id = f"DECOY_id{precursorMz}"
        variableDict = {
            "precursorMz": precursorMz,
            "peptide": f"peptide{precursorMz}",
            "precursorCharge": 1,
            "id": id,
            "protein": protein,
        }
        for i in range(1, 11):
            peakMz = precursorMz - (precursorMzDiff / 2) + (i - 1) * peakMzDiff
            peakIntensity = i * peakIntensityDiff
            data.append(
                [
                    variableDict["precursorMz"],
                    variableDict["peptide"],
                    variableDict["precursorCharge"],
                    variableDict["id"],
                    variableDict["protein"],
                    peakMz,
                    peakIntensity,
                ]
            )

    df = pd.DataFrame(
        data,
        columns=[
            "precursorMz",
            "peptide",
            "precursorCharge",
            "id",
            "protein",
            "peakMz",
            "peakIntensity",
        ],
    )
    return df


def make_mgf_library_from_template_library_dataframe(templateLibraryDf):
    mgfFormattedSpectra = []
    for key, df in templateLibraryDf.groupby(["peptide", "precursorCharge"]):
        spectrumDict = {}
        paramsDict = {}
        paramsDict["pepmass"] = (df["precursorMz"].iat[0], None)
        paramsDict["seq"] = key[0]
        paramsDict["charge"] = int(key[1])
        paramsDict["title"] = df["id"].iat[0]
        paramsDict["protein"] = df["protein"].iat[0]
        spectrumDict["params"] = paramsDict
        spectrumDict["m/z array"] = np.array(df["peakMz"])
        spectrumDict["intensity array"] = np.array(df["peakIntensity"])
        mgfFormattedSpectra.append(spectrumDict)
    return mgfFormattedSpectra


spectrastColumns = [
    "PrecursorMz",
    "FullUniModPeptideName",
    "PrecursorCharge",
    "transition_group_id",
    "ProteinName",
    "ProductMz",
    "LibraryIntensity",
]

fragpipeColumns = [
    "PrecursorMz",
    "ModifiedPeptideSequence",
    "PrecursorCharge",
    "PeptideSequence",
    "ProteinId",
    "ProductMz",
    "LibraryIntensity",
]

prositColumns = [
    "PrecursorMz",
    "ModifiedPeptide",
    "PrecursorCharge",
    "StrippedPeptide",
    "FragmentLossType",
    "FragmentMz",
    "RelativeIntensity",
]
