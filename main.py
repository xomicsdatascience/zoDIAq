import pandas as pd
import cosine_similarity_functions as csf
pd.set_option("display.max_rows", None, "display.max_columns", None)

'''
#  For slicing the library spectra into smaller files
## full list
temp_df = pd.read_csv('Part 2/consensus_transitions_lowres_decoys.tsv', sep='\t')

## for sequences
lib_df = temp_df.loc[temp_df['FullUniModPeptideName']=='APIIAVTR']

## for ranges
#lib_df = temp_df.query( 'FullUniModPeptideName == 418 and PrecursorMz < 423' )

## write to file
lib_df.to_csv('Part 2/condensed_APIIAVTR.csv')
'''

lib_file = 'condensed_APIIAVTR.csv'
exp_file = '20190411_DI2A_1to16_n1b.mzXML'

lib = csf.tramlFileConversionCSV(lib_file)
final_df = csf.expSpectraAnalysis( exp_file, lib )
final_df.to_csv( 'Data/cosine_output1.csv' )
