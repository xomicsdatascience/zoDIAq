from csodiaq.loaders.LibraryLoaderStrategy import (
    LibraryLoaderStrategy,
    create_peaks_from_mz_intensity_lists_and_csodiaq_key_id,
    remove_low_intensity_peaks_below_max_peak_num,
    finalVariableNames,
)
from csodiaq.loaders.mappings import newColumns, oldColumnsSpectrast, oldColumnsFragpipe, oldColumnsProsit
import pandas as pd
import os


class LibraryLoaderStrategyTable(LibraryLoaderStrategy):
    """
    Concrete strategy implementation of the LibraryLoaderStrategy strategy class specific for
        loading comma- or tab-delimited library files (tables) into a standardized dictionary
        for use in CsoDIAq.

    Extended Summary
    ----------------
    While table-based libraries can all be read into a pandas dataframe, the column names can differ.
        This class accounts for THREE variations of library column titles (as derived from
        SpectraST/PanHuman, FragPipe, and Prosit outputs).
        The 'old' column name list of any of these three formats is converted into a single standard
        of 'new' column name list to avoid needing to switch back and forth through the code.

    Attributes
    ----------
    rawLibDf : pd.DataFrame
        The original table as loaded from the .tsv or .csv library file. The relevant columns of this
        table are renamed ('old' -> 'new') and irrelevant columns are removed during the reformatting.
    oldToNewColumnDict : dict
        A dictionary for mapping old column names to standardized new column names. This dictionary
        will differ between different library types (titled 'spectrast', 'fragpipe' and 'prosit'),
        as the library types have different 'old' column names.
    """

    def _load_raw_library_object_from_file(self, libraryFilePath: os.PathLike) -> None:
        if libraryFilePath.endswith(".tsv"):
            separator = "\t"
        else:
            separator = ","
        self.rawLibDf = pd.read_csv(libraryFilePath, sep=separator)
        self.oldToNewColumnDict = _set_old_to_new_column_dict(libraryFilePath)
        _assert_there_are_no_missing_columns(
            self.oldToNewColumnDict.keys(), self.rawLibDf.columns
        )

    def _format_raw_library_object_into_csodiaq_library_dict(self) -> dict:
        """abstract class implementation - see 'LibraryLoaderStrategy.py' for details"""
        maxPeakNum = 10
        reformattedLibDf = _reformat_raw_library_object_columns(
            self.rawLibDf, self.oldToNewColumnDict
        )
        organizedDataDict = _organize_data_by_csodiaq_library_dict_keys(
            reformattedLibDf
        )
        csodiaqLibraryDict = {}
        for csodiaqKeyIdx in range(len(organizedDataDict["csodiaqKeys"])):
            csodiaqKey, csodiaqValue = _create_csodiaq_library_entry(
                organizedDataDict, maxPeakNum, csodiaqKeyIdx
            )
            csodiaqLibraryDict[csodiaqKey] = csodiaqValue
        return csodiaqLibraryDict


def _assert_there_are_no_missing_columns(
    requiredColumns: list, presentColumns: list
) -> None:
    missingColumnValues = set(requiredColumns) - set(presentColumns)
    if len(missingColumnValues):
        raise ValueError(
            f'table library file is missing expected column(s). Missing values: [{", ".join(sorted(missingColumnValues))}])'
        )


def _reformat_raw_library_object_columns(
    df: pd.DataFrame, oldToNewColumnDict: dict
) -> pd.DataFrame:
    reformattedDf = df[oldToNewColumnDict.keys()]
    reformattedDf = reformattedDf.rename(columns=oldToNewColumnDict)
    reformattedDf["csodiaqLibKey"] = list(
        zip(
            reformattedDf["precursorMz"].tolist(), reformattedDf["peptideName"].tolist()
        )
    )
    return reformattedDf


def _organize_data_by_csodiaq_library_dict_keys(df: pd.DataFrame) -> dict:
    keys = sorted(set(df["csodiaqLibKey"]))
    mz = df.groupby("csodiaqLibKey")["peakMz"].apply(list).to_dict()
    intensities = df.groupby("csodiaqLibKey")["peakIntensity"].apply(list).to_dict()
    df.drop_duplicates(subset="csodiaqLibKey", inplace=True)
    df.set_index("csodiaqLibKey", drop=True, inplace=True)
    df.drop(
        ["precursorMz", "peakMz", "peptideName", "peakIntensity"], axis=1, inplace=True
    )
    metadata = df.to_dict(orient="index")
    return {
        "csodiaqKeys": keys,
        "mz": mz,
        "intensities": intensities,
        "metadata": metadata,
    }


def _create_csodiaq_library_entry(
    organizedDataDict: dict, maxPeakNum: int, csodiaqKeyIdx: int
) -> dict:
    csodiaqKey = organizedDataDict["csodiaqKeys"][csodiaqKeyIdx]
    peaks = create_peaks_from_mz_intensity_lists_and_csodiaq_key_id(
        organizedDataDict["mz"][csodiaqKey],
        organizedDataDict["intensities"][csodiaqKey],
        csodiaqKeyIdx,
    )
    reducedPeaks = remove_low_intensity_peaks_below_max_peak_num(peaks, maxPeakNum)
    isDecoy = int(
        "decoy" in organizedDataDict["metadata"][csodiaqKey]["proteinName"].lower()
    )
    return csodiaqKey, {
        finalVariableNames["precursorCharge"]: organizedDataDict["metadata"][
            csodiaqKey
        ]["precursorCharge"],
        finalVariableNames["identifier"]: organizedDataDict["metadata"][csodiaqKey][
            "identifier"
        ],
        finalVariableNames["proteinName"]: organizedDataDict["metadata"][csodiaqKey][
            "proteinName"
        ],
        finalVariableNames["peaks"]: sorted(reducedPeaks),
        finalVariableNames["csodiaqKeyIdx"]: csodiaqKeyIdx,
        finalVariableNames["isDecoy"]: isDecoy,
    }


def _set_old_to_new_column_dict(filePath) -> None:
    librarySource = _determine_library_source_from_file(filePath)
    if librarySource == "spectrast":
        oldColumns = oldColumnsSpectrast
        return dict(zip(oldColumns, newColumns))
    elif librarySource == "fragpipe":
        oldColumns = oldColumnsFragpipe
        return dict(zip(oldColumns, newColumns))
    elif librarySource == "prosit":
        oldColumns = oldColumnsProsit
        return dict(zip(oldColumns, newColumns))


def _determine_library_source_from_file(filePath) -> str:
    with open(filePath) as f:
        columns = f.readline()
        if "transition_group_id" in columns:
            return "spectrast"
        elif "ProteinId" in columns:
            return "fragpipe"
        elif "RelativeIntensity" in columns:
            return "prosit"
        else:
            raise ValueError(
                "The library table file provided does not match spectrast, fragpipe, or prosit formats."
            )
