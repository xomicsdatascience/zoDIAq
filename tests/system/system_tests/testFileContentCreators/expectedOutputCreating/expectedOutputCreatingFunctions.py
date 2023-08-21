import numpy as np
import pandas as pd


def make_output_df_from_spectra_breakdown(spectraBreakdown):
    data = []
    for cosineGroup in spectraBreakdown:
        cosine = cosineGroup["cosine"]
        numQueryPeaks = find_number_of_total_query_peaks_in_scan(
            cosineGroup["querySpectra"]
        )
        for i in range(len(cosineGroup["libSpectra"])):
            libSpectra = cosineGroup["libSpectra"][i]
            querySpectra = cosineGroup["querySpectra"][i]
            queryMzLargerThanPrecursorIdx = np.argwhere(
                np.array(querySpectra["mzs"]) > cosineGroup["queryScanMz"]
            )
            queryMzLargerThanPrecursorIntensities = np.array(
                querySpectra["intensities"]
            )[queryMzLargerThanPrecursorIdx]
            ionCount = 0
            if len(queryMzLargerThanPrecursorIntensities):
                ionCount = sum(queryMzLargerThanPrecursorIntensities)[0]
            if "DECOY" in libSpectra["protein"]:
                decoy = 1
            else:
                decoy = 0
            data.append(
                [
                    cosineGroup["scan"],
                    cosineGroup["queryScanMz"],
                    libSpectra["peptide"],
                    libSpectra["protein"],
                    decoy,
                    libSpectra["precursorMz"],
                    libSpectra["precursorCharge"],
                    cosine,
                    libSpectra["id"],
                    numQueryPeaks,
                    len(libSpectra["mzs"]),
                    len(libSpectra["mzs"]),
                    ionCount,
                    -cosineGroup["scan"],
                    cosineGroup["windowWidth"],
                    len(libSpectra["mzs"]) - len(queryMzLargerThanPrecursorIntensities),
                ]
            )
    return pd.DataFrame(
        data,
        columns=[
            "scan",
            "MzEXP",
            "peptide",
            "protein",
            "isDecoy",
            "MzLIB",
            "zLIB",
            "cosine",
            "name",
            "Peak(Query)",
            "Peaks(Library)",
            "shared",
            "ionCount",
            "CompensationVoltage",
            "totalWindowWidth",
            "exclude_num",
        ],
    )


def find_number_of_total_query_peaks_in_scan(spectra):
    lengths = [len(spectrum["mzs"]) for spectrum in spectra]
    return sum(lengths)
