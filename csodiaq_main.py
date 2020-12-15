import pandas as pd
import csodiaq_menu_functions as menu
import csodiaq_figure_functions as figure
import csodiaq_base_functions as cbf
import os
from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2
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

#menu.write_csodiaq_fdr_outputs(outFile, corrected=True)

caleb_peptide_df = pd.read_csv('Data/Output/csodiaq_lib-human-5peaks-noloss-pt2mz-400to2000_exp-n1b_corrected_2SD_peptideFDR.csv')
caleb_protein_df = pd.read_csv('Data/Output/csodiaq_lib-human-5peaks-noloss-pt2mz-400to2000_exp-n1b_corrected_2SD_proteinFDR.csv')
jesse_peptide_df = pd.read_csv('Data/Input/peptide_matches_Jesse.csv')
jesse_protein_df = pd.read_csv('Data/Input/protein_matches_Jesse.csv')

caleb_peptides = sorted(list(caleb_peptide_df['peptide']))
caleb_proteins = sorted(list(set(caleb_protein_df['protein'])),reverse=True)
jesse_peptides = sorted(list(jesse_peptide_df['Peptide']))
jesse_proteins = sorted(list(jesse_protein_df['Protein']))

print(caleb_peptides[:10])
print(jesse_peptides[:10])
print(caleb_proteins[:10])
print(jesse_proteins[:10])

#UniMod:1 = +42.01057
#UniMod:4 = +57.0215
#UniMod:5 = +43.0058
#UniMod:35 = +15.9949


unimodDict = {
    '(UniMod:4)':'+57.0215',
    '(UniMod:5)':'+42.01057',
    '(UniMod:35)':'+15.9949'
}

for i in range(len(caleb_peptides)):
    caleb_peptides[i] = re.sub('\(UniMod:\d+\)','',caleb_peptides[i])

for i in range(len(jesse_peptides)):
    jesse_peptides[i] = re.sub('\+\d+\.\d+','',jesse_peptides[i])


print(len(caleb_proteins))
print(len(set(caleb_proteins)))

caleb_proteins = [re.sub('(.*)\|(.*)\|(.*_HUMAN)(.*)', r'\2_\3', x) for x in caleb_proteins]
jesse_proteins = [re.sub('(.*)\|(.*)\|(.*_HUMAN)(.*)', r'\2_\3', x) for x in jesse_proteins]

print(caleb_proteins[4])

protDict = {}
for x in caleb_proteins:
    if x not in protDict:
        protDict[x] = 1
    else:
        protDict[x] += 1
for x in protDict:
    if protDict[x] > 1: print(x)

print(len(caleb_proteins))
print(len(set(caleb_proteins)))

venn2([set(caleb_peptides),set(jesse_peptides)], set_labels = ["Caleb's Output", "Jesse's Output"])
plt.title('Comparing Peptide Identification Outputs\n')
plt.savefig('peptide_comparison_venn.png')

plt.clf()

venn2([set(caleb_proteins),set(jesse_proteins)], set_labels = ["Caleb's Output", "Jesse's Output"])
plt.title('Comparing Protein Identification Outputs\n')
plt.savefig('protein_comparison_venn.png')
