import pandas as pd
import numpy as np
import os


def get_parent_dir():
    return os.path.dirname(os.path.abspath(__file__))


def create_template_library_dataframe(
    minMz, maxMz, precursorMzDiff, peakMzDiff, peakIntensityDiff
):
    data = []
    for precursorMz in range(minMz, maxMz, precursorMzDiff):
        if precursorMz < (maxMz + minMz) / 2:
            protein = f"1/protein{precursorMz}"
        else:
            protein = f"1/DECOY_protein{precursorMz}"
        variableDict = {
            "precursorMz": precursorMz,
            "peptide": f"peptide{precursorMz}",
            "precursorCharge": 1,
            "id": f"id{precursorMz}",
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
