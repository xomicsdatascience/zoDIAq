import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from collections import defaultdict as dd
from pandas.api.types import is_numeric_dtype
import pyteomics.mzxml
from csodiaq.spectrum import Spectrum
import re
import itertools

# library_peptide_column_name = 'PeptideSequence'  # for later expansion to other library formats
library_peptide_column_name = 'FullUniModPeptideName'
library_precursormz_column_name = 'PrecursorMz'  # m/z of precursor
library_productmz_column_name = 'ProductMz'      # m/z of fragment
library_intensity_column_name = 'LibraryIntensity'  # intensity of fragment
csodiaq_query_scan_name = 'scan'  # CsoDIAq column name for scan id
csodiaq_mz_lib_name = 'MzLIB'  # CsoDIAq column name for library precursor m/z
csodiaq_protein_column_name = 'protein'
experimental_query_mz_name = 'm/z array'  # mz name for the input to CsoDIAq
experimental_query_intensity_name = 'intensity array'  # intensity name for the input to CsoDIAq
output_peptide_format = 'peptide_quantity_{x}'  # format for output column

def get_peptide_quantities(file_list: list,
                           library_file: str,
                           csodiaq_output_dir: str,
                           num_library_fragments: int = 10,
                           save_file: str = None):
    """
    Calculates the quantities of peptides that returned by CsoDIAq.
    Parameters
    ----------
    file_list : list
        List of MZXML files containing the scan data.
    library_file : str
        Path to the library file to use as reference. Expected to be .tsv
    csodiaq_output_dir : str
        Path to the output directory of CsoDIAq, containing the _proteinFDR.csv files.
    num_library_fragments : int
        Number of library fragments to use for matching. The fragments are ordered by library fragment intensity;
        10 fragments indicate that only the 10 most intense fragments will be used for quantification.
    save_file : str
        Optional. Path where to store DataFrame. If None, output is returned instead of saved.

    Returns
    -------
    None
    """

    # Get peptides predicted by CsoDIAq
    common_dataframe = extract_common_entries(csodiaq_output_dir,
                                              input_file_list=file_list,
                                              common_col=['peptide'],
                                              load_columns=[csodiaq_query_scan_name,
                                                            csodiaq_mz_lib_name, 'peptide', 'ionCount'],
                                              normalize=False)
    if common_dataframe.shape[0] == 0:
        return
    # Load library
    library_dataframe = pd.read_csv(library_file, sep='\t')

    # Get experimental spectra
    exp_data = []
    for f in file_list:
        exp_data.append(pyteomics.mzxml.read(f))

    # initialize new dataframe
    file_names_for_dataframe = format_filenames(file_list)
    # peptide_quant_df = pd.DataFrame(columns=[output_peptide_format.format(x=i) for i in range(len(exp_data))])
    peptide_quant_df = pd.DataFrame(columns=file_names_for_dataframe)

    # Go through each peptide common across files; quantify
    len_frame = common_dataframe.shape[0]
    for peptide_idx, peptide_dat in common_dataframe.iterrows():
        peptide = peptide_idx[0]
        precursor_mz = peptide_idx[1]

        # Go through each (pre-loaded) input, fetch relevant scan, quantify peptide
        peptide_intensity_data = []  # use to collect intensity to later append to peptide_quant_df

        peptide_library_spectrum = get_library_spectrum_for_peptide(library_dataframe,
                                                                    peptide=peptide,
                                                                    precursor_mz=precursor_mz,
                                                                    num_fragments=num_library_fragments)

        for exp_idx, exp in enumerate(exp_data):
            # Get scan id
            scan_id = peptide_dat[f'scan_{exp_idx}']
            scan_data = exp.get_by_id(str(int(scan_id)))
            scan_spectrum = Spectrum(mz=scan_data[experimental_query_mz_name],
                                     intensity=scan_data[experimental_query_intensity_name])
            # Get intensity at matched m/z
            scan_idx, library_idx = scan_spectrum.get_matching_mz_indices(spectrum_to_match=peptide_library_spectrum,
                                                                          match_tolerance_ppm=30)
            peptide_intensity_data.append(sum(scan_spectrum.intensity[scan_idx]))
        peptide_quant_df.loc[peptide] = peptide_intensity_data

    # Get mean + std across files
    peptide_quant_df['mean'] = peptide_quant_df.apply(np.mean, axis=1)
    peptide_quant_df['std'] = peptide_quant_df.apply(np.std, axis=1)
    # peptide_quant_df = peptide_quant_df[peptide_quant_df['mean'] > 0]
    # peptide_quant_df.drop(peptide_quant_df)
    if save_file is not None:
        peptide_quant_df.to_csv(save_file, index=True, header=True)
    else:
        return peptide_quant_df
    return


def get_protein_quantities(file_list: list,
                           library_file: str,
                           csodiaq_output_dir: str,
                           num_library_fragments: int = 10,
                           save_file: str = None) -> None:
    """
    Calculates the quantities of peptides that returned by CsoDIAq.
    Parameters
    ----------
    file_list : list
        List of MZXML files containing the scan data.
    library_file : str
        Path to the library file to use as reference. Expected to be .tsv
    csodiaq_output_dir : str
        Path to the output directory of CsoDIAq, containing the _proteinFDR.csv files.
    num_library_fragments : int
        Number of library fragments to use for matching. The fragments are ordered by library fragment intensity;
        10 fragments indicate that only the 10 most intense fragments will be used for quantification.
    save_file : str
        Optional. Path where to store DataFrame. If None, output is returned instead of saved.

    Returns
    -------
    None
    """

    common_dataframe = extract_common_entries(csodiaq_output_dir,
                                              input_file_list=file_list,
                                              common_col=['peptide', 'MzLIB', 'protein'],
                                              load_columns=[csodiaq_query_scan_name,
                                                            csodiaq_mz_lib_name, 'protein', 'ionCount', 'peptide'],
                                              normalize=False)
    common_dataframe.reset_index(inplace=True)
    common_dataframe.set_index(['peptide','MzLIB'], inplace=True)
    if common_dataframe.shape[0] == 0:
        return  # no common peptides
    # Get library spectra
    library_dataframe = pd.read_csv(library_file, sep='\t')

    # Get experimental spectra
    exp_data = []
    for f in file_list:
        exp_data.append(pyteomics.mzxml.read(f))  # These are the input files to CsoDIAq

    # initialize new dataframe
    file_names_for_dataframe = format_filenames(file_list)  # this reduces the column names

    protein_quant_df = pd.DataFrame(columns=file_names_for_dataframe)

    # iterate through each peptide that is common across the data; quantify each one
    # Get list of unique proteins
    protein_set = set(common_dataframe['protein'])
    for protein in protein_set:
        # Get all rows that correspond to this protein
        protein_df = common_dataframe[common_dataframe['protein'] == protein]
        # Go through each row; quantify all peptides associated with the protein
        protein_intensity_data = np.zeros(len(exp_data))
        for peptide_idx, peptide_dat in protein_df.iterrows():
            peptide = peptide_idx[0]
            precursor_mz = peptide_idx[1]
            peptide_library_spectrum = get_library_spectrum_for_peptide(library_dataframe,
                                                                        peptide=peptide,
                                                                        precursor_mz=precursor_mz,
                                                                        num_fragments=num_library_fragments)
            # Go through each experiment
            for exp_idx, exp in enumerate(exp_data):
                # Get scan
                scan_id = peptide_dat[f"scan_{exp_idx}"]
                scan_data = exp.get_by_id(str(int(scan_id)))
                scan_spectrum = Spectrum(mz=scan_data[experimental_query_mz_name],
                                         intensity=scan_data[experimental_query_intensity_name])
                # Get intensity at matched m/z
                scan_idx, library_idx = scan_spectrum.get_matching_mz_indices(
                    spectrum_to_match=peptide_library_spectrum,
                    match_tolerance_ppm=30)
                protein_intensity_data[exp_idx] += sum(scan_spectrum.intensity[scan_idx])
        protein_quant_df.loc[protein] = protein_intensity_data
        # Get mean + std across files
    protein_quant_df['mean'] = protein_quant_df.apply(np.mean, axis=1)
    protein_quant_df['std'] = protein_quant_df.apply(np.std, axis=1)

    if save_file is not None:
        protein_quant_df.to_csv(save_file, index=True, header=True)
    return protein_quant_df


def get_possible_protein_keys(protein_str: str) -> list:
    """
    Returns the n! factorial possible reorderings of the protein synonyms
    Parameters
    ----------
    protein_str : str
        String formatted for protein names (e.g. 3/sp|Q01813|PFKAP_HUMAN/sp|P17858|PFKAL_HUMAN/sp|P08237|PFKAM_HUMAN).

    Returns
    -------
    list
        List of possible re-orderings of the protein synonyms
    """

    # First value indicates how many synonyms there are
    synonyms = protein_str.split('/')
    num_synonyms = int(synonyms[0])
    if num_synonyms == 1:
        return [protein_str]  # only one name; done
    possible_orders = itertools.permutations(synonyms[1:], num_synonyms)
    protein_name_list = []
    # Construct reordered protein name
    for order in possible_orders:
        protein_name_list.append(f"{num_synonyms}/{'/'.join(order)}")
    return protein_name_list


def format_filenames(file_list: list) -> list:
    """
    Shortens the filenames and ensures that they're unique
    Parameters
    ----------
    file_list : list
        List of strings to format

    Returns
    -------
    list
        Lits of formatted strings
    """
    common_prefix = os.path.commonprefix(file_list)
    len_pre = len(common_prefix)
    return [f[len_pre:] for f in file_list]


def get_library_spectrum_for_peptide(library_dataframe: pd.DataFrame,
                                     peptide: str,
                                     precursor_mz: float,
                                     num_fragments: int = 10) -> Spectrum:
    """
    Returns the Spectrum object corresponding to specified peptide and precursor M/z.
    Parameters
    ----------
    library_dataframe : pd.DataFrame
        DataFrame containing the reference library.
    peptide : str
        Peptide sequence to match (corresponds to however the 'PeptideSequence' column is formatted in the library).
    precursor_mz : float
        Precursor M/z for the peptide.

    Returns
    -------
    Spectrum
        Spectrum of the matched peptide.
    """
    # Figure out which rows are relevant
    # First get precursor m/z (numerical comparison is faster)
    df_condition = library_dataframe[library_precursormz_column_name] == precursor_mz
    peptide_df = library_dataframe[df_condition]

    # Get matching peptide from smaller dataframe
    df_condition = peptide_df[library_peptide_column_name] == peptide
    peptide_df = peptide_df[df_condition]

    # We know that there should be a peptide here because CsoDIAq used the library to get the peptide in the first place
    # There is some bug where (UniMod:x) tags result in a small numerical error in the CsoDIAq library Mz value.
    # We'll check the library again, but with small bounds on the precursor_mz value
    if peptide_df.shape[0] == 0:  # If there are no matches, get
        print(peptide)
        df_condition = library_dataframe[library_precursormz_column_name] >= precursor_mz - 1e-5
        peptide_df = library_dataframe[df_condition]

        df_condition = peptide_df[library_precursormz_column_name] <= precursor_mz + 1e-5
        peptide_df = peptide_df[df_condition]

        df_condition = peptide_df[library_peptide_column_name] == peptide
        peptide_df = peptide_df[df_condition]

    # Extract values
    library_intensity = peptide_df[library_intensity_column_name].values
    library_mz = peptide_df[library_productmz_column_name].values
    return Spectrum(mz=library_mz, intensity=library_intensity, num_fragments=num_fragments)


def extract_common_entries(csodiaq_output_dir: str,
                           input_file_list: list,
                           common_col: list = ('peptide', 'MzLIB'),
                           load_columns: list = None,
                           normalize: bool = True) -> pd.DataFrame:
    """
    Looks at the specified column across multiple files and reduces to the values that are common
    across all files.
    Parameters
    ----------
    csodiaq_output_dir : str
        Directory containing the output from CsoDIAq.
    input_file_list : list
        List of files input to CsoDIAq
    common_col : str
        Column containing the value to reduce across the input files. Default: 'peptide'
    load_columns : list
        List of columns to load.
    normalize : bool
        Whether to normalize across matched columns. Default: True
    Returns
    -------
    None
    """
    # Get the output files
    fdr_files = gather_matching_proteinfdr_files(csodiaq_output_dir=csodiaq_output_dir,
                                                 file_list=input_file_list)
    if len(fdr_files) == 0:
        return pd.DataFrame()
    # Load and combine DataFrames
    common_dataframe = load_multi_dataframe(fdr_files, index_col=common_col, load_columns=load_columns)

    # Normalize quantities across files; columns differentiated by only the trailing "_#" are assumed to match
    if normalize:
        common_dataframe = normalize_columns(common_dataframe, inplace=True)
    return common_dataframe


def gather_matching_proteinfdr_files(csodiaq_output_dir: str,
                                     file_list: list) -> list:
    """
    Gathers the _proteinFDR files matching the input list. The order in which the files are returned matches the order
    in file_list.
    Parameters
    ----------
    csodiaq_output_dir : str
        Path to CsoDIAq's output directory.
    file_list : list
        List of files used as input to CsoDIAq.

    Returns
    -------
    list
        List of matched _proteinFDR.csv files
    """
    matched_file_list = [None for _ in file_list]
    # start_idx = 0
    for output in os.listdir(csodiaq_output_dir):
        for file_idx, file in enumerate(file_list):
            if _is_proteinfdr_match(file, output):
                matched_file_list[file_idx] = os.path.join(csodiaq_output_dir, output)
                # start_idx = file_idx
                break
    return matched_file_list

def _is_proteinfdr_match(input_file: str, protein_fdr_file: str):
    """Checks whether the two files match"""
    if '_' not in protein_fdr_file:
        return False
    fdr_basename = os.path.basename(protein_fdr_file)
    fdr_split = fdr_basename.split('_')[1]
    input_basename = os.path.basename(input_file)
    last_extension_sep_idx = input_basename.rfind(os.path.extsep)
    if last_extension_sep_idx == -1:
        last_extension_sep_idx = len(input_basename)  # file may not have extension; rfind returns -1; fix it
    input_ref_name = input_basename[:last_extension_sep_idx]  # remove .mzxml (or whatever format)
    # Check that input filename is in potential match
    return (input_ref_name == fdr_split) and (fdr_basename.endswith('_proteinFDR.csv'))


def gather_proteinfdr_files(output_dir: str) -> list:
    """
    Goes through the output_directory and returns the files ending in "FDR"
    Parameters
    ----------
    output_dir : str
        Path to the CsoDIAq output directory
    Returns
    -------
    list
        List of FDR files
    """
    fdr_files = []
    for f in os.listdir(output_dir):
        if f.endswith('proteinFDR.csv'):
            fdr_files.append(os.path.join(output_dir, f))
    return fdr_files


def load_multi_dataframe(csv_list: list,
                         index_col: str,
                         load_columns: list = None,
                         merge_how='inner') -> pd.DataFrame:
    """
    Loads each DataFrame in the list, merges them, and returns the result.
    Parameters
    ----------
    csv_list : list
        List of csv files to load into DataFrames
    index_col : str
        Name of the column to use as an index; this determines how the DataFrames are merged.
    load_columns : list
        Columns to load. Specifying the columns of interest will speed up loading and reduce the memory footprint.
    merge_how : str
        How to merge the DataFrames. 'inner' is recommended; this will reduce the returned DataFrame to contain
        only the elements common across all DataFrames.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the loaded and merged files.
    """

    merged_df = None
    print(f"index_col: {index_col}")
    for file_idx, file in enumerate(csv_list):
        # load_columns == None results in default behaviour (load all)
        # loaded_df = pd.read_csv(file, usecols=load_columns, index_col=index_col[0])  # PATCH
        print(file)
        loaded_df = pd.read_csv(file, usecols=load_columns, index_col=index_col)
        # loaded_df.reset_index(inplace=True)
        # loaded_df.set_index(list(index_col), inplace=True)
        # Rename columns so they'll be unique after merge (since we're looping, merge alterations will get messy
        rename_cols = {}
        for col in loaded_df.columns:
            rename_cols[col] = f"{col}_{file_idx}"
        loaded_df.rename(columns=rename_cols, inplace=True)

        # rename peptides with (Unimod: x)


        # Merge current file with previous one
        if merged_df is None:
            merged_df = loaded_df
        else:
            # merged_df = merged_df.merge(loaded_df, on=index_col, how=merge_how).drop_duplicates()
            merged_df = merged_df.merge(loaded_df, left_index=True, right_index=True, how=merge_how).drop_duplicates()
        # if 'peptide' in load_columns or 'peptide' in index_col:
        #     pattern='\(UniMod:[0-9]*\)'
        #     merged_df.reset_index(inplace=True)
        #     merged_df['peptide'] = merged_df['peptide'].apply(lambda x: re.sub(pattern, '', x))
        #     merged_df.set_index(index_col, inplace=True)
    return merged_df

def normalize_columns(df: pd.DataFrame,
                      inplace: bool = True) -> pd.DataFrame:
    """
    Divides by the maximum value across corresponding columns; columns that differ only by the trailing "_#"
    (e.g. "_2" and "_3") are assumed to correspond to one another.
    Parameters
    ----------
    df : DataFrame
        DataFrame with columns to normalize.
    inplace : bool
        Whether to perform the normalizing in place. If True, returned DataFrame will be a pointer to the input
         (modified) DataFrame. Default: True
    Returns
    -------
    pd.DataFrame
        DataFrame with columns normalized across corresponding columns.
    """

    # Match columns
    matched_columns = _match_columns(list(df.columns))
    if inplace:
        df_to_use = df
    else:
        df_to_use = df.copy()
    # Iterate through each set, find max, then divide
    for col_list in matched_columns:
        # Skip if not numeric
        if not is_numeric_dtype(df_to_use[col_list]):
            continue
        max_across_cols = df[col_list].apply(np.max)
        df_to_use[col_list] /= max_across_cols
    return df_to_use

def _match_columns(columns: list) -> list:
    """
    Return lists of columns that correspond to one another (differ only in the trailing underscore).
    Parameters
    ----------
    columns : list
        List of columns to match

    Returns
    -------
    list
        List of lists, each containing columns matching one another
    """

    match_dict = dd(list)
    matched_columns = []
    for col in columns:
        # First get common name
        split_col = col.split('_')[:-1]  # last entry is the differentiating info
        # Recover previous underscores
        common_col = '_'.join(split_col)
        match_dict[common_col].append(col)  # All entries with common col will be listed here

    # Extract lists from dict
    for _, col_list in match_dict.items():
        matched_columns.append(col_list)
    return matched_columns