import pandas as pd
import csodiaq_menu_functions as menu
import csodiaq_figure_functions as figure
import os


# experimental file
expFile = 'Data/Input/20190411_DI2A_1to16_n1b.mzXML'

# subset
#libFile = 'Data/Input/condensed31_APIIAVTR-PLAEGTPR.csv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-31peaks-APIIAVTR-PLAEGTPR_exp-n1b_sqrt.csv'

# full - 31
libFile = 'Data/Input/iproph-speclib_con_decoys20.tsv'
outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks_exp-n1b_sqrt.csv'

# full - 20
#libFile = 'Data/Input/iproph-speclib_con_decoys20.tsv'
#outFile = 'Data/Output/csodiaq_lib-iproph1-20peaks_exp-n1b_sqrt.csv'

#menu.draw_histogram('Data/Output/csodiaq_lib-iproph1-31peaks_exp-n1b_sqrt_corrected1.csv')


menu.write_csodiaq_output(libFile, expFile, outFile, 3)
menu.filter_optimal_match_csodiaq_output(outFile)
menu.write_ppm_spread(outFile)
menu.write_ppm_offset_sd(outFile)
figure.draw_histogram(outFile)


c=1
menu.write_csodiaq_output(libFile, expFile, outFile, 3, correctedNumSD=c)
menu.filter_optimal_match_csodiaq_output(outFile, correctedNumSD=c)
menu.write_ppm_spread(outFile, correctedNumSD=c)
menu.write_ppm_offset_sd(outFile, correctedNumSD=c)
figure.draw_histogram(outFile, correctedNumSD=c)
