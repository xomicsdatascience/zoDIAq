import pandas as pd
import os
from pyteomics import mgf
import re
from csodiaq.loaders.library_format_columns import naming_list, data_extraction
import numpy as np

"""
This file consists of the loaders for the different types of library files.
"""


def load_library(library_file: os.PathLike,
                 max_peaks: int = 10) -> dict:
    """
    Loads a library file and returns the contents as a dict. Attempts to automatically identify the type of file
    from the file extension and/or headers in the data.

    Parameters
    ----------
    library_file : os.PathLike
        Pathlike to the library file.
    max_peaks : int
        Number of peaks to retain.
    Returns
    -------
    dict
        Dictionary containing the data entries.
    """

    # Try to identify filetype by file extension
    library_file_str = str(library_file)
    if(library_file_str.endswith('.mgz')):
        return mgf_library_upload(library_file_str)
    elif(library_file_str.endswith(('.csv', '.tsv'))):
        # Table; identify standard, convert to CsoDIAq-compatible format
        if(library_file_str.endswith('.csv')):
            data = pd.read_csv(library_file)
        elif(library_file_str.endswith('.tsv')):
            data = pd.read_csv(library_file, sep='\t')
        data['LibraryIntensity'] = data['LibraryIntensity'].apply(np.sqrt, axis=1)
        _remap_table_columns(data)  # Renames columns to CsoDIAq standard

    # data table now has expected columns
    # Need to add new columns
    data['ID'] = list(zip(data['PrecursorMz'].tolist(),
                            data['FullUniModPeptideName'].tolist()))
    mz = data.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensities = data.groupby("ID")['LibraryIntensity'].apply(list).to_dict()
    data.drop_duplicates(subset="ID", inplace=True)
    data.set_index('ID', drop=True, inplace=True)

    data = data.to_dict(orient='index')
    idx = 0
    for k in data.keys():
        idx += 1
        mz_pairs = zip(intensities[k], mz[k])
        pair_sorted = sorted(mz_pairs, reverse=True)  # greatest intensities first; use mz as secondary sort index
        # NOTE: This favours peaks that have a higher m/z; if two peaks are equal intensities, the one with
        # higher m/z will be returned first.
        peak_intensities, peak_mzs = zip(*pair_sorted)  # separate the sorted values back into separate lists
        peak_tuple = [(peak_mz, peak_intensity, idx) for peak_mz, peak_intensity in zip(peak_mzs[:max_peaks], peak_intensities[:max_peaks])]
        data[k]['Peaks'] = peak_tuple
        data[k]['ID'] = idx

        data[k]['Decoy'] = int('decoy' in data[k]['ProteinName'].lower())  # set 'Decoy' to 0/1 if protein is decoy
    return data

def mgf_library_upload(file_name: str,
                       max_peaks: int = 10) -> dict:
    """
    Loads Mascot Generic Format (MGF) files and returns the contents as a formatted dict compatible with CsoDIAq.
    Parameters
    ----------
    file_name : str
        Filepath of the library file to load.
    max_peaks : int
        Maximum number of peaks to retain for a spectrum. Default: 10.
    Returns
    -------

    """
    libMGF = mgf.read(file_name)
    # smf.print_milestone('Enter library dictionary upload: ')
    lib = {}
    id = 0
    for spec in libMGF:
        id += 1
        key = (spec['params']['pepmass'][0], spec['params']['seq'])
        charge = int(re.sub('[+-]', '', str(spec['params']['charge'][0])))
        name = spec['params']['title']
        if 'protein' in spec['params']:
            protein = spec['params']['protein']
        else:
            protein = ''
        if 'DECOY' in name:
            decoy = 1
        else:
            decoy = 0
        mz = spec['m/z array']
        intensity = spec['intensity array']
        intensity = [x**0.5 for x in intensity]
        keyList = [id for x in mz]
        peaks = list(tuple(zip(mz, intensity, keyList)))
        peaks.sort(key=lambda x: x[1], reverse=True)
        if len(peaks) > max_peaks:
            peaks = peaks[:max_peaks]
        peaks.sort(key=lambda x: x[0])
        tempDict = {
            'PrecursorCharge': charge,
            'transition_group_id': name,
            'ProteinName': protein,
            'Peaks': peaks,
            'ID': id,
            'Decoy': decoy,
        }
        lib[key] = tempDict
    return lib


def _remap_table_columns(data: pd.DataFrame) -> None:
    """
    Modifies the input table columns in-place to have the CsoDIAq naming scheme.

    Parameters
    ----------
    data : pd.DataFrame
        Table with columns to rename.

    Returns
    -------
    None
    """
    col_mapper, conversion_func = _get_renaming(data.columns)
    conversion_func(data, col_mapper)
    return


def _get_renaming(cols: list):
    """
    Identifies the type of table from the column names and returns the renaming dict such that
    rename_dict['column_in_data'] -> 'csodiaq_col_name'.

    Parameters
    ----------
    cols : list
        List of column names to check.

    Returns
    -------
    rename_dict : dict
        Renaming dict with entries rename_dict['column_in_data'] -> 'csodiaq_col_name'

    Raises
    ------
    ValueError
        If the list of columns doesn't match any known library format.
    """
    # Convert to set for easy comparison
    cols_set = set(cols)
    idx = 0
    for reference_mapping, conversion_func in zip(naming_list, data_extraction):
        ref_set = set(reference_mapping.keys())
        if(ref_set.issubset(cols_set)):  # Are all required columns in the data's columns?
            return reference_mapping, conversion_func
        idx += 1
    else:
        raise ValueError(f'The input library does not match a known template. Please check the data.')


def _traml_column_headings(columns):
    if 'FullUniModPeptideName' in columns:  # SpectraST
        return {
            'type': 'SpectraST',
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'FullUniModPeptideName',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'transition_group_id',
            'ProteinName': 'ProteinName',
        }
    else:  # pan human
        return {
            'type': 'PanHuman',
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'ModifiedPeptideSequence',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'TransitionGroupId',
            'ProteinName': 'ProteinId',
        }