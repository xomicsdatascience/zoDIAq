import pandas as pd
import csodiaq_menu_functions as menu
import csodiaq_figure_functions as figure
import csodiaq_base_functions as cbf
import os
from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np
import re
from os import listdir
from os.path import isfile, join
from pyteomics import mzxml, mgf
import statistics

print("Index,TimeElapsed,NumPeaks,ExpPrecursorMz,NumLibrarySpectra,outputFile")
'''
mypath = 'Data/Input/100mzxml/'
onlyfiles = [mypath+f for f in listdir(mypath) if isfile(join(mypath, f))]
expDict = {}
for x in onlyfiles:
    key = re.sub('Data/Input/100mzxml/20200719_MAGIC_MCF7_1128repro_(\d+)\.mzXML', r'exp-100reps-rep\1', x)
    expDict[key] = x
'''
expDict = {
'exp-n1b':'Data/Input/20190411_DI2A_1to16_n1b.mzXML',
#'exp-MCF7':'Data/Input/20190405_MCF7_FAIMS_18_2.mzXML'
}


libDict = {
'lib-human-noloss-400to2000-pt2mz-31peaks_':'Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv'
#'lib-faims-mgf_':'Data/Input/human.faims.fixed.decoy.mgf'
#'lib-human-noloss-400to2000-pt2mz-allTop_':'Data/Input/human_31peaks_noloss_400to2000_pt2mz_allTop.csv'
#'lib-human-noloss-400to2000-pt2mz-goodTop_':'Data/Input/human_31peaks_noloss_400to2000_pt2mz_goodTop.csv'
}

for libTag in libDict:
    libFile = libDict[libTag]
#    lib = cbf.library_file_to_dict(libFile)

    for exp in sorted(expDict):
        expFile = expDict[exp]
        outFile = 'Data/Output/csodiaq_'+libTag+exp+'.csv'
#        menu.write_csodiaq_output(lib, expFile, outFile)
#        menu.write_ppm_offset_tolerance(outFile, hist=True)

#        menu.write_csodiaq_output(lib, expFile, outFile, corrected=True)
#        menu.write_csodiaq_fdr_outputs(outFile, corrected=True)
#        menu.write_DISPA_targeted_reanalysis_files(outFile)

        menu.heavy_light_quantification(outFile, libFile, 'Data/Input/jesse/')
