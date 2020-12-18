import pandas as pd
import csodiaq_menu_functions as menu
import csodiaq_figure_functions as figure
import csodiaq_base_functions as cbf
import os
from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np
import re

# experimental file
expFile = 'Data/Input/20190411_DI2A_1to16_n1b.mzXML'
#expFile = 'Data/Input/20190405_MCF7_FAIMS_18_2.mzXML'

# subset
#libFile = 'Data/Input/condensed31_APIIAVTR-PLAEGTPR.csv'
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


# Library Variants

## full - 20 - highres
##libFile = 'Data/Input/iproph-speclib_con_decoys20_highres.tsv'
##outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks-highres_exp-n1b.csv'

## full - 20 - no loss, pt2mz
##libFile = 'Data/Input/iproph-speclib_con_decoys20_noloss_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks-noloss-pt2mz_exp-n1b.csv'

## full - 20 - no loss, pt2mz, 300to2000
##libFile = 'Data/Input/human_20peaks_noloss_300to2000_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-20peaks-noloss-pt2mz-300to2000_exp-n1b.csv'

## full - 20 - no loss, pt2mz, 400to2000
##libFile = 'Data/Input/human_20peaks_noloss_400to2000_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-20peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

## full - 10 - no loss
##libFile = 'Data/Input/human_10peaks_noloss.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss_exp-n1b.csv'

## full - 10 - no loss, 400to2000
##libFile = 'Data/Input/human_10peaks_noloss_400to2000.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss-400to2000_exp-n1b.csv'

## full - 10 - no loss, pt2mz, 400to2000
##libFile = 'Data/Input/human_10peaks_noloss_400to2000_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-10peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

## full - 5 - no loss, pt2mz, 400to2000
libFile = 'Data/Input/human_5peaks_noloss_400to2000_pt2mz.tsv'
outFile = 'Data/Output/csodiaq_lib-human-5peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

# full - MGF
##libFile = 'Data/Input/human.faims.fixed.decoy.mgf'
##outFile = 'Data/Output/csodiaq_lib31-mgf-human-faims_exp-n1b.csv'

## subset - MGF
##libFile = 'Data/Input/condensed_APIIAVTR.mgf'
##outFile = 'Data/Output/csodiaq_lib31-mgf-human-faims-APIIAVTR_exp-n1b.csv'

# full - 31 - no loss, pt2mz, 400to2000
##libFile = 'Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-31peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

## full - 31 to 20 test - no loss, pt2mz, 400to2000
##libFile = 'Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv'
##outFile = 'Data/Output/csodiaq_lib-human-31to20peaks-noloss-pt2mz-400to2000_exp-n1b.csv'

'''
#menu.write_csodiaq_output(libFile, expFile, outFile)
menu.filter_optimal_match_csodiaq_output(outFile)
menu.write_ppm_spread(outFile)
menu.write_ppm_offset_tolerance(outFile, hist=True)

menu.write_csodiaq_output(libFile, expFile, outFile, corrected=True)
menu.filter_optimal_match_csodiaq_output(outFile, corrected=True)
menu.write_ppm_spread(outFile, corrected=True)
menu.write_ppm_offset_tolerance(outFile, corrected=True, hist=True)
'''
#menu.write_csodiaq_fdr_outputs(outFile)
menu.write_csodiaq_fdr_outputs(outFile, corrected=True)

#figure.create_venn_diagrams()
