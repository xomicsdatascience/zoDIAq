import pandas as pd
import cosine_similarity_functions as csf
pd.set_option("display.max_rows", None, "display.max_columns", None)
from timeit import default_timer as timer
from datetime import timedelta
import random as rand

'''
#  For slicing the library spectra into smaller files
## full list
temp_df = pd.read_csv('consensus_transitions_lowres_decoys.tsv', sep='\t')
## for sequences
#lib_df = temp_df.loc[(temp_df['FullUniModPeptideName']=='FDPSSVSEWEPHWR')] #| (temp_df['FullUniModPeptideName']=='PLAEGTPR')]

## for ranges
lib_df = temp_df.query( 'PrecursorMz > 422 and PrecursorMz < 423' )

## write to file
lib_df.to_csv('condensed_422to423.csv')

'''

lib_file = 'condensed_APIIAVTR_PLAEGTPR.csv'
lib_file2 = 'iproph-speclib_con_decoys31.tsv'
exp_file = '20190411_DI2A_1to16_n1b.mzXML'
out_file = 'Data/fullOutput_lib31_match3_ppm10.csv'

t0 = timer()
print('#Enter lib upload/conversion:')
print('#'+str(timedelta(seconds=t0)))
lib = csf.traml_library_upload_csv(lib_file2)
t1 = timer()
print('#enter spectra comparison:')
print('#'+str(timedelta(seconds=t1)))
final_df = csf.query_spectra_analysis( exp_file, out_file, lib, 3, 10 )
t2 = timer()
print('#done:')
print('#'+str(timedelta(seconds=t2)))
#final_df.to_csv( 'Data/cosine_output2.csv' )
