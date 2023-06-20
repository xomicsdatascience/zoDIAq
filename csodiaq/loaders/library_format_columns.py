"""
This file specifies how to map the columns of different standards onto the internal one used by CsoDIAq.
"""
import pandas as pd


spectrast_cols = {
    'PrecursorMz': 'PrecursorMz',
    'FullUniModPeptideName': 'FullUniModPeptideName',
    'PrecursorCharge': 'PrecursorCharge',
    'ProductMz': 'ProductMz',
    'LibraryIntensity': 'LibraryIntensity',
    'transition_group_id': 'transition_group_id',
    'ProteinName': 'ProteinName'
}

panhuman_cols = {
    'PrecursorMz': 'PrecursorMz',
    'ModifiedPeptideSequence': 'FullUniModPeptideName',
    'PrecursorCharge': 'PrecursorCharge',
    'ProductMz': 'ProductMz',
    'LibraryIntensity': 'LibraryIntensity',
    'TransitionGroupId': 'transition_group_id',
    'ProteinId': 'ProteinName'
}

fragpipe_cols = {
    'PrecursorMz': 'PrecursorMz',
    'ModifiedPeptideSequence': 'FullUniModPeptideName',
    'PrecursorCharge': 'PrecursorCharge',
    'ProductMz': 'ProductMz',
    'LibraryIntensity': 'LibraryIntensity',
    'ProteinId': 'ProteinName'
}


def rename(dataframe: pd.DataFrame,
           col_rename_dict: dict) -> None:
    '''
    Renames the columns in `dataframe` according to renaming dict. Modifications are made in-place.
    Parameters
    ----------
    dataframe : pd.DataFrame
        Dataframe with columns to rename
    col_rename_dict : dict
        Mapping dictionary for renaming columns
    Returns
    -------
    None
    '''
    dataframe.rename(columns=col_rename_dict, inplace=True)
    return


def fragpipe_conv(dataframe: pd.DataFrame,
                  col_rename_dict: dict) -> None:
    '''
    Examines the dataframe and extracts the relevant information to create the expected format for the column.
    Modifications are made in-place.
    Parameters
    ----------
    dataframe : pd.DataFrame
        DataFrame to with the data to extract/rename.
    col_rename_dict : dict
        Mapping dictionary for renaming columns

    Returns
    -------
    None
    '''
    # Only need to create transition_group_id:
    # [idx]_[b/y][position]_[fragment_charge]_[peptide_sequence]_[precursor_charge]
    fstring = '0_{ion_type}{frag_pos}_{frag_charge}_{peptide_sequence}_{precursor_charge}'
    row_extract = lambda row: fstring.format(ion_type=row.Annotation[0],
                                             frag_pos=row.Annotation[1:row.Annotation.index('^')],
                                             frag_charge=row.Annotation[row.Annotation.index('^')+1:],
                                             peptide_sequence=row.PeptideSequence,
                                             precursor_charge=row.PrecursorCharge)
    dataframe['transition_group_id'] = dataframe.apply(row_extract, axis=1)
    rename(dataframe, col_rename_dict)
    return

naming_list = [spectrast_cols, panhuman_cols, fragpipe_cols]
data_extraction = [rename, rename, fragpipe_conv]
