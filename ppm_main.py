import pandas as pd
import ppm_optimization_step as pos
pd.set_option("display.max_rows", None, "display.max_columns", None)
from timeit import default_timer as timer
from datetime import timedelta
import random as rand
import csv

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
folder = 'Data/Output/'

# subset
#tag = 'subset_APIIAVTR-PLAEGTPR_lib31_ppmOptimization_sqrt_'
#libFile = 'Data/Input/condensed31_APIIAVTR-PLAEGTPR.csv'

# full
tag = 'full_lib31_ppmOptimization_sqrt_'
libFile = 'Data/Input/iproph-speclib_con_decoys31.tsv'

# everything else
expFile = 'Data/Input/20190411_DI2A_1to16_n1b.mzXML'
outFile = 'Output.csv'
ppmFullFile = 'Ppm.csv'
histogramFile = 'Histogram.png'

t0 = timer()
print('#Enter lib upload/conversion:')
print('#'+str(timedelta(seconds=t0)))
lib = pos.traml_library_upload_csv(libFile)

t1 = timer()
print('#enter spectra comparison:')
print('#'+str(timedelta(seconds=t1)))
pos.query_spectra_analysis( expFile,
                            folder + tag + outFile,
                            lib,
                            3,
                            10,
                            folder + tag + ppmFullFile )

t2 = timer()
print('#enter FDR and offset calculation:')
print('#'+str(timedelta(seconds=t2)))
df = pd.read_csv(folder + tag + outFile)
bestMatchNum, bestFDR = pos.find_best_matchNum_fdr(df, 0.01)
offset, sd = pos.find_query_offset_sd(folder + tag + outFile,
                                        folder + tag + ppmFullFile,
                                        bestMatchNum,
                                        bestFDR,
                                        folder + tag + histogramFile)
print('#offset: '+str(offset))
print('#standard deviation: '+str(sd))

t3 = timer()
print('#enter corrected spectra comparison:')
print('#'+str(timedelta(seconds=t3)))
pos.query_spectra_analysis( expFile,
                            folder +'corrected_'+ tag + outFile,
                            lib,
                            3,
                            2*sd,
                            folder +'corrected_'+ tag + ppmFullFile,
                            ppmYOffset = -offset )

t2 = timer()
print('#enter corrected FDR and offset calculation:')
print('#'+str(timedelta(seconds=t2)))
df = pd.read_csv(folder +'corrected_'+ tag + outFile)
bestMatchNum, bestFDR = pos.find_best_matchNum_fdr(df, 0.01)
offset, sd = pos.find_query_offset_sd(folder +'corrected_'+ tag + outFile,
                                    folder +'corrected_'+ tag + ppmFullFile,
                                    bestMatchNum,
                                    bestFDR,
                                    folder + 'corrected_' + tag + histogramFile)
print('#offset: '+str(offset))
print('#standard deviation: '+str(sd))
