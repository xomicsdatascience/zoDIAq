import pandas as pd
import os
from pyteomics import mgf
import re
from csodiaq.loaders.library_format_columns import naming_list, data_extraction
from typing import Callable

"""
This file consists of the loaders for the different types of library files.
"""


def load_library(library_file: os.PathLike) -> dict:
    """
    Loads a library file and returns the contents as a dict. Attempts to automatically identify the type of file
    from the file extension and/or headers in the data.

    Parameters
    ----------
    library_file : os.PathLike
        Pathlike to the library file.

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
        _remap_table_columns(data)  # Renames columns to CsoDIAq standard
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


def traml_library_upload(file_name):
    if file_name.endswith('.tsv'):
        lib_df = pd.read_csv(file_name, sep='\t')
    else:
        lib_df = pd.read_csv(file_name)
    # smf.print_milestone('Enter library dictionary upload: ')

    # Pan human and spectraST libraries have different column names. This normalizes the columns.
    headings = _traml_column_headings(lib_df.columns)
    lib_df = lib_df.loc[:, lib_df.columns.intersection([headings['PrecursorMz'], headings['FullUniModPeptideName'], headings['PrecursorCharge'],
                                                        headings['ProductMz'], headings['LibraryIntensity'], headings['transition_group_id'], headings['ProteinName']])]

    lib_df = lib_df[[headings['PrecursorMz'], headings['FullUniModPeptideName'], headings['PrecursorCharge'],
                     headings['ProductMz'], headings['LibraryIntensity'], headings['transition_group_id'], headings['ProteinName']]]
    lib_df.columns = ['PrecursorMz', 'FullUniModPeptideName', 'PrecursorCharge',
                      'ProductMz', 'LibraryIntensity', 'transition_group_id', 'ProteinName']

    lib_df['LibraryIntensity'] = [
        x**0.5 for x in list(lib_df['LibraryIntensity'])]
    lib_df['ID'] = list(zip(lib_df['PrecursorMz'].tolist(),
                            lib_df['FullUniModPeptideName'].tolist()))

    mz_dict = lib_df.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensity_dict = lib_df.groupby(
        "ID")['LibraryIntensity'].apply(list).to_dict()
    lib_df.drop_duplicates(subset="ID", inplace=True)
    lib_df = lib_df.loc[:, lib_df.columns.intersection(
        ['ID', 'PrecursorCharge', 'transition_group_id', 'ProteinName'])]
    lib_df.set_index("ID", drop=True, inplace=True)
    lib = lib_df.to_dict(orient="index")

    # pan human library formats are different, including how the peptides are matched to proteins (esp. decoys). This section of code adjusts for this discrepancy.
    if headings['type'] == 'PanHuman':
        for key, value in lib.items():
            proteins = lib[key]['ProteinName'].split('/')
            num = proteins.pop(0)
            newProteins = [x for x in proteins if 'DECOY' not in x]
            proteinStr = str(len(newProteins))
            for x in newProteins:
                if 'DECOY' in num:
                    proteinStr += ('/DECOY_'+x)
                else:
                    proteinStr += ('/'+x)
            lib[key]['ProteinName'] = proteinStr

    id = 0
    for key in lib:
        id += 1
        mz, intensity = (list(t) for t in zip(
            *sorted(zip(mz_dict[key], intensity_dict[key]))))
        keyList = [id for i in range(len(mz))]
        peaks = list(tuple(zip(mz, intensity, keyList)))

        peaks.sort(key=lambda x: x[1], reverse=True)
        if len(peaks) > 10:
            peaks = peaks[:10]

        peaks.sort(key=lambda x: x[0])
        lib[key]['Peaks'] = peaks
        lib[key]['ID'] = id
        if 'DECOY' in lib[key]['ProteinName']:
            lib[key]['Decoy'] = 1
        else:
            lib[key]['Decoy'] = 0
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

    for reference_mapping, conversion_func in zip(naming_list, data_extraction):
        ref_set = set(reference_mapping.keys())
        if(ref_set.issubset(cols_set)):  # Are all required columns in the data's columns?
            return reference_mapping, conversion_func
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
