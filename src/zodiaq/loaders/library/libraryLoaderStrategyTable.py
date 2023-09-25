from zodiaq.loaders.library.libraryLoaderStrategy import (
    LibraryLoaderStrategy,
    create_peaks_from_mz_intensity_lists_and_zodiaq_key_id,
    remove_low_intensity_peaks_below_max_peak_num,
    finalVariableNames,
)
from zodiaq.loaders.library.mappings import (
    newColumns,
    oldColumnsSpectrast,
    oldColumnsFragpipe,
    oldColumnsProsit,
)
import pandas as pd
import os


class LibraryLoaderStrategyTable(LibraryLoaderStrategy):
    """
    Concrete strategy implementation of the LibraryLoaderStrategy strategy class specific for
        loading comma- or tab-delimited library files (tables) into a standardized dictionary
        for use in zoDIAq.

    Extended Summary
    ----------------
    While table-based library can all be read into a pandas dataframe, the column names can differ.
        This class accounts for THREE variations of library column titles (as derived from
        SpectraST/PanHuman, FragPipe, and Prosit outputs).
        The 'old' column name list of any of these three formats is converted into a single standard
        of 'new' column name list to avoid needing to switch back and forth through the code.

    Attributes
    ----------
    rawUploadedLibraryObject : pd.DataFrame
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
        self.rawUploadedLibraryObject = pd.read_csv(libraryFilePath, sep=separator)
        self.oldToNewColumnDict = _set_old_to_new_column_dict(libraryFilePath)
        _assert_there_are_no_missing_columns(
            self.oldToNewColumnDict.keys(), self.rawUploadedLibraryObject.columns
        )

    def _format_raw_library_object_into_zodiaq_library_dict(self) -> dict:
        """abstract class implementation - see 'libraryLoaderStrategy.py' for details"""
        maxPeakNum = 10
        reformattedLibDf = _reformat_raw_library_object_columns(
            self.rawUploadedLibraryObject, self.oldToNewColumnDict
        )
        organizedDataDict = _organize_data_by_zodiaq_library_dict_keys(
            reformattedLibDf
        )
        zodiaqLibraryDict = {}
        for zodiaqKeyIdx in range(len(organizedDataDict["zodiaqKeys"])):
            zodiaqKey, zodiaqValue = _create_zodiaq_library_entry(
                organizedDataDict, maxPeakNum, zodiaqKeyIdx
            )
            zodiaqLibraryDict[zodiaqKey] = zodiaqValue
        return zodiaqLibraryDict


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
    reformattedDf["zodiaqLibKey"] = list(
        zip(
            reformattedDf["precursorMz"].tolist(), reformattedDf["peptideName"].tolist()
        )
    )
    return reformattedDf


def _organize_data_by_zodiaq_library_dict_keys(df: pd.DataFrame) -> dict:
    keys = sorted(set(df["zodiaqLibKey"]))
    mz = df.groupby("zodiaqLibKey")["peakMz"].apply(list).to_dict()
    intensities = df.groupby("zodiaqLibKey")["peakIntensity"].apply(list).to_dict()
    df.drop_duplicates(subset="zodiaqLibKey", inplace=True)
    df.set_index("zodiaqLibKey", drop=True, inplace=True)
    df.drop(
        ["precursorMz", "peakMz", "peptideName", "peakIntensity"], axis=1, inplace=True
    )
    metadata = df.to_dict(orient="index")
    return {
        "zodiaqKeys": keys,
        "mz": mz,
        "intensities": intensities,
        "metadata": metadata,
    }


def _create_zodiaq_library_entry(
    organizedDataDict: dict, maxPeakNum: int, zodiaqKeyIdx: int
) -> dict:
    zodiaqKey = organizedDataDict["zodiaqKeys"][zodiaqKeyIdx]
    peaks = create_peaks_from_mz_intensity_lists_and_zodiaq_key_id(
        organizedDataDict["mz"][zodiaqKey],
        organizedDataDict["intensities"][zodiaqKey],
        zodiaqKeyIdx,
    )
    reducedPeaks = remove_low_intensity_peaks_below_max_peak_num(peaks, maxPeakNum)
    isDecoy = int(
        "decoy" in organizedDataDict["metadata"][zodiaqKey]["proteinName"].lower()
    )
    return zodiaqKey, {
        finalVariableNames["precursorCharge"]: organizedDataDict["metadata"][
            zodiaqKey
        ]["precursorCharge"],
        finalVariableNames["identification"]: organizedDataDict["metadata"][zodiaqKey][
            "identification"
        ],
        finalVariableNames["proteinName"]: organizedDataDict["metadata"][zodiaqKey][
            "proteinName"
        ],
        finalVariableNames["peaks"]: sorted(reducedPeaks),
        finalVariableNames["zodiaqKeyIdx"]: zodiaqKeyIdx,
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
