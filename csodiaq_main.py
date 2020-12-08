import pandas as pd
import csodiaq_menu_functions as menu
import csodiaq_figure_functions as figure
import csodiaq_base_functions as cbf
import os
from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np

# experimental file
expFile = 'Data/Input/20190411_DI2A_1to16_n1b.mzXML'

# subset
libFile = 'Data/Input/condensed31_APIIAVTR-PLAEGTPR.csv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-31peaks-APIIAVTR-PLAEGTPR_exp-n1b_sqrt.csv'

# full - 31
#libFile = 'Data/Input/iproph-speclib_con_decoys31.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-31peaks_exp-n1b_sqrt.csv'

# full - 20
#libFile = 'Data/Input/iproph-speclib_con_decoys20.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks_exp-n1b_sqrt.csv'

# full - 10
#libFile = 'Data/Input/iproph-speclib_con_decoys10.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-10peaks_exp-n1b.csv'

# full - 20 - highres
#libFile = 'Data/Input/iproph-speclib_con_decoys20_highres.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks-highres_exp-n1b.csv'

# full - 20 - no loss, pt2mz
#libFile = 'Data/Input/iproph-speclib_con_decoys20_noloss_pt2mz.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks-noloss-pt2mz_exp-n1b.csv'

# full - 20 - no loss, pt2mz, 300to2000
#libFile = 'Data/Input/human_20peaks_noloss_300to2000_pt2mz.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-20peaks-noloss-pt2mz-300to2000_exp-n1b.csv'

# full - 20 - no loss, pt2mz, 400to2000
#libFile = 'Data/Input/human_20peaks_noloss_400to2000_pt2mz.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-20peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

# full - 10 - no loss
#libFile = 'Data/Input/human_10peaks_noloss.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss_exp-n1b.csv'

# full - 10 - no loss, 400to2000
#libFile = 'Data/Input/human_10peaks_noloss_400to2000.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss-400to2000_exp-n1b.csv'

# full - 10 - no loss, pt2mz, 400to2000
#libFile = 'Data/Input/human_10peaks_noloss_400to2000_pt2mz.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

# full - 5 - no loss, pt2mz, 400to2000
#libFile = 'Data/Input/human_5peaks_noloss_400to2000_pt2mz.tsv'
#outFile = 'Data/Output/csodiaq_lib-human-5peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

# full - MGF
#libFile = 'Data/Input/human.faims.fixed.decoy.mgf'
#outFile = 'Data/Output/csodiaq_lib10-mgf-human-faims_exp-n1b.csv'

# subset - MGF
libFile2 = 'Data/Input/condensed_APIIAVTR.mgf'
outFile = 'Data/Output/csodiaq_lib10-mgf-human-faims_exp-n1b.csv'


'''
menu.write_csodiaq_output(libFile, expFile, outFile)
menu.filter_optimal_match_csodiaq_output(outFile)
menu.write_ppm_spread(outFile)
menu.write_ppm_offset_tolerance(outFile, hist=True)

menu.write_csodiaq_output(libFile, expFile, outFile, corrected=True)
menu.filter_optimal_match_csodiaq_output(outFile, corrected=True)
menu.write_ppm_spread(outFile, corrected=True)
menu.write_ppm_offset_tolerance(outFile, corrected=True, hist=True)
'''

'''
lib1 = cbf.traml_library_upload_csv(libFile)
key1 = list(lib1.keys())[0]
lib2 = cbf.mgf_library_upload(libFile2)
key2 = list(lib2.keys())[0]
def plot_spec(SPECTRA, COLOR):
    plt.vlines(SPECTRA.columns, np.repeat(0, len(SPECTRA.columns)), SPECTRA, colors=COLOR)
print(key1)
print(key2)

peak1 = lib1[key1]['Peaks']
peak2 = lib2[key2]['Peaks']

df1 = pd.DataFrame(columns=[x[0] for x in peak1])
df2 = pd.DataFrame(columns=[x[0] for x in peak2])
df1.loc[0] = [x[1] for x in peak1]
df2.loc[0] = [x[1] for x in peak2]

plot_spec(df1, 'red')
plot_spec(-df2, 'green')
plt.show()
'''
