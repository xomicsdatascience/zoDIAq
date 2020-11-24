import pandas as pd
import cosine_similarity_functions as csf
pd.set_option("display.max_rows", None, "display.max_columns", None)
from timeit import default_timer as timer
from datetime import timedelta

'''
#  For slicing the library spectra into smaller files
## full list
temp_df = pd.read_csv('consensus_transitions_lowres_decoys.tsv', sep='\t')

## for sequences
lib_df = temp_df.loc[(temp_df['FullUniModPeptideName']=='APIIAVTR') | (temp_df['FullUniModPeptideName']=='PLAEGTPR')]

## for ranges
#lib_df = temp_df.query( 'FullUniModPeptideName == 418 and PrecursorMz < 423' )

## write to file
lib_df.to_csv('condensed_APIIAVTR_PLAEGTPR.csv')
'''

lib_file = 'condensed_APIIAVTR_PLAEGTPR.csv'
lib_file2 = 'consensus_transitions_lowres_decoys.tsv'
exp_file = '20190411_DI2A_1to16_n1b.mzXML'
out_file = 'Data/fullOutput.csv'

t0 = timer()
print('Enter lib upload/conversion:')
print(timedelta(seconds=t0))
lib = csf.tramlFileConversionCSV(lib_file2)
t1 = timer()
print('enter spectra comparison:')
print(timedelta(seconds=t1))
final_df = csf.expSpectraAnalysis( exp_file, out_file, lib )
t2 = timer()
print('done')
print(timedelta(seconds=t2))
#final_df.to_csv( 'Data/cosine_output2.csv' )

'''

t1 = [1, 4, 3, 5]
t2 = [2, 6, 7, 8]
t3 = [9, 2, 6, 1]
t4 = list(tuple(zip(t1,t2,t3)))
#t4 = list(sorted(zip(t1, t2, t3)))))
print(t4)
'''
'''
o1 = pd.read_csv('Data/output2_fullLib_fullSpect_directWrite-removeDf.csv')
o2 = pd.read_csv('Data/output_fullLib_fullSpect_combineLib.csv')

top = min(len(o1.index),len(o2.index))

for i in range(top):
    v1 = o1.loc[i]['cosine']
    v2 = o2.loc[i]['cosine']
    if v1 == v2:
        print()
        print('original: '+str(v1))
        print('new: '+str(v2))
'''
