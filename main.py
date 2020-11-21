import pandas as pd
import cosine_similarity_functions as csf
pd.set_option("display.max_rows", None, "display.max_columns", None)
from timeit import default_timer as timer
from datetime import timedelta

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
lib_file2 = 'consensus_transitions_lowres_decoys.tsv'
exp_file = '20190411_DI2A_1to16_n1b.mzXML'
out_file = 'Data/output2_APIIAVTR_fullSpect_directWrite-removeDf.csv'


t0 = timer()
print('Enter lib upload/conversion:')
print(timedelta(seconds=t0))
lib = csf.tramlFileConversionCSV(lib_file)
t1 = timer()
print('enter spectra comparison:')
print(timedelta(seconds=t1))
final_df = csf.expSpectraAnalysis( exp_file, out_file, lib )
t2 = timer()
print('done')
print(timedelta(seconds=t2))
#final_df.to_csv( 'Data/cosine_output2.csv' )
