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


def make_mz_window_spectra_summary_for_library_spectra(targetDecoyLibrarySpectra):
    return [
        {
            "title": "highCosineTargets",
            "scan": 1,
            "cosine": 0.99,
            "libSpectra": targetDecoyLibrarySpectra["targets"][:200],
        },
        {
            "title": "highCosineDecoys",
            "scan": 2,
            "cosine": 0.98,
            "libSpectra": targetDecoyLibrarySpectra["decoys"][:2],
        },
        {
            "title": "midCosineTargets",
            "scan": 3,
            "cosine": 0.97,
            "libSpectra": targetDecoyLibrarySpectra["targets"][200:400],
        },
        {
            "title": "midCosineDecoys",
            "scan": 4,
            "cosine": 0.96,
            "libSpectra": targetDecoyLibrarySpectra["decoys"][2:5],
        },
        {
            "title": "lowCosineTargets",
            "scan": 5,
            "cosine": 0.95,
            "libSpectra": targetDecoyLibrarySpectra["targets"][400:],
        },
        {
            "title": "lowCosineDecoys",
            "scan": 6,
            "cosine": 0.94,
            "libSpectra": targetDecoyLibrarySpectra["decoys"][5:],
        },
    ]


def separate_target_and_decoy_library_spectra(libraryDf):
    targets = []
    decoys = []
    for vars, df in libraryDf.groupby(
        [
            "precursorMz",
            "peptide",
            "precursorCharge",
            "id",
            "protein",
        ]
    ):
        if "DECOY" in vars[-1]:
            peptideList = decoys
        else:
            peptideList = targets
        tempDict = {}
        tempDict["precursorMz"] = vars[0]
        tempDict["peptide"] = vars[1]
        tempDict["precursorCharge"] = vars[2]
        tempDict["id"] = vars[3]
        tempDict["protein"] = vars[4]
        tempDict["mzs"] = np.array(df["peakMz"])
        tempDict["intensities"] = np.array(df["peakIntensity"])
        peptideList.append(tempDict)
    return {"targets": targets, "decoys": decoys}


def add_query_file_components_to_mz_window_spectra_summary(mzWindowSpectraBreakdown):
    for mzWindowGroup in mzWindowSpectraBreakdown:
        (
            mzWindowGroup["queryScanMz"],
            mzWindowGroup["windowWidth"],
        ) = find_mz_window_of_query_scan_from_spectra(mzWindowGroup["libSpectra"])
        mzWindowGroup["querySpectra"] = [
            create_query_spectrum_from_library_spectrum_and_cosine_score(
                spectrum, mzWindowGroup["cosine"]
            )
            for spectrum in mzWindowGroup["libSpectra"]
        ]
    return mzWindowSpectraBreakdown


def find_mz_window_of_query_scan_from_spectra(spectra):
    allPrecursorMzs = [spectrum["precursorMz"] for spectrum in spectra]
    minMz = min(allPrecursorMzs)
    maxMz = max(allPrecursorMzs)
    centralMz = (minMz + maxMz) / 2
    return (centralMz, (centralMz - minMz) * 2 + 0.0001)


def find_number_of_total_query_peaks_in_scan(spectra):
    lengths = [len(spectrum["mzs"]) for spectrum in spectra]
    return sum(lengths)


def create_vector_that_can_be_used_to_create_cosine_score(vectorA, cosineScore):
    unitVector = vectorA / np.linalg.norm(vectorA)
    positiveVector = np.array(np.full((len(vectorA)), 1000.0))
    perpendicularVector = positiveVector - positiveVector.dot(unitVector) * unitVector
    perpendicularVector = perpendicularVector / np.linalg.norm(perpendicularVector)
    vectorB = (
        cosineScore * unitVector + np.sqrt(1 - cosineScore**2) * perpendicularVector
    )
    return vectorB


def create_query_spectrum_from_library_spectrum_and_cosine_score(spectrum, cosineScore):
    return {
        "mzs": spectrum["mzs"],
        "intensities": create_vector_that_can_be_used_to_create_cosine_score(
            spectrum["intensities"], cosineScore
        ),
    }


spectrastColumns = [
    "PrecursorMz",
    "FullUniModPeptideName",
    "PrecursorCharge",
    "transition_group_id",
    "ProteinName",
    "ProductMz",
    "LibraryIntensity",
]
