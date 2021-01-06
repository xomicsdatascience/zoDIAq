import pandas as pd
import re
import csodiaq_base_functions as cbf
import csv
from statistics import pstdev, median
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib_venn import venn2

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_meta_analysis_files(inFile, corrected=False):
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    overallDf = pd.read_csv(inFile).sort_values('cosine', ascending=False)
    overallDf = overallDf.reset_index(drop=True)

    outFileSpectrum = re.sub('Data/Output/(.*).csv', r'Data/Figures/FDRGraphs/\1_FDRGraph_spectral.csv', inFile)
    outFilePeptide = re.sub('Data/Output/(.*).csv', r'Data/Figures/FDRGraphs/\1_FDRGraph_peptide.csv', inFile)
    outFileProtein = re.sub('Data/Output/(.*).csv', r'Data/Figures/FDRGraphs/\1_FDRGraph_protein.csv', inFile)
    match_fdrSat_graph(overallDf, 0.01, outFileSpectrum)
    match_fdrSat_graph(overallDf, 0.01, outFilePeptide)
    match_fdrSat_graph(overallDf, 0.01, outFileProtein)
'''
    peptideDf = overallDf.copy()
    peptideDf = peptideDf.drop_duplicates(subset='peptide', keep='first')
    peptideDf = peptideDf.reset_index(drop=True)

    proteinDf = overallDf.copy()
    print(len(proteinDf))
    proteinDf.drop_duplicates(subset='protein', keep='first')
    print(len(proteinDf))
    proteinDf = proteinDf.reset_index(drop=True)
    print('-------------')
'''

'''
Function:
Purpose:
Parameters:
Returns:
'''
def match_fdrSat_graph(df, fdrMax, outFile, DEBUG=False):
    matches = sorted(list(set(df['shared'])))

    if DEBUG:
        print('****************************************')
        print(outFile)
        print('****************************************')
    fdrSat = []
    numPeaks = []
    cosine = []
    subsetDecoys = []
    totalDecoys = []
    for m in matches:

        if 'peptide' in outFile: tempDf = cbf.add_fdr_to_csodiaq_output(df, filterType='peptide', bestMatchNum=m)
        elif 'protein' in outFile: tempDf = cbf.generate_protein_csodiaq_fdr_output(df, bestMatchNum=m); tempDf = tempDf.drop_duplicates(subset='leadingProtein', keep='first'); tempDf = tempDf.reset_index(drop=True)
        else: tempDf = cbf.add_fdr_to_csodiaq_output(df, bestMatchNum=m)

        dec = [x for x in tempDf['protein'] if 'DECOY' in x]
        fdrSat.append(len(tempDf))
        subsetDecoys.append(len(dec))
        if len(tempDf) != 0: cos = tempDf['cosine'].loc[len(tempDf)-1]
        else: cos = 0
        cosine.append(cos)

    graphDf = pd.DataFrame({'matches':matches, 'FDRCutoff': fdrSat, 'cosine': cosine, 'FDRDecoys': subsetDecoys})
    if DEBUG:
        print('----------------------------')
        print(graphDf)
    graphDf.to_csv(outFile, index=False)

'''
Function:
Purpose:
Parameters:
Returns:
'''
def draw_histogram(inFile, correctedNumSD=0):
    ppmSpreadFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    with open(ppmSpreadFile, newline='') as f: ppmList = list(csv.reader(f))[0]
    ppmList = [float(x) for x in ppmList]
    data = [ppmList]
    labels = ['uncorrected']
    mean = (sum(ppmList)/len(ppmList))
    sd = pstdev(ppmList)
    if correctedNumSD:
        inFile = re.sub('(.*).csv', r'\1_corrected'+str(correctedNumSD)+'.csv', inFile)
        ppmSpreadFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
        with open(ppmSpreadFile, newline='') as f: ppmList = list(csv.reader(f))[0]
        ppmList = [float(x) for x in ppmList]
        data.append(ppmList)
        labels.append('corrected')
        mean = (sum(ppmList)/len(ppmList))
        sd = pstdev(ppmList)

    plt.figure(figsize=(8,6))
    plt.hist(data, bins=200, alpha=0.8, density=True, label=labels)
    plt.xlabel("PPM Difference", size=14)
    plt.ylabel("Count", size=14)
    plt.title("Histogram of PPM spread")
    plt.legend(loc='upper right')
    plt.axvline(x=mean, color='black', linestyle = 'dashed')
    plt.axvline(x=sd+mean, color='green', linestyle = 'dashed')
    plt.axvline(x=(-sd+mean), color='green', linestyle = 'dashed')
    plt.axvline(x=2*sd+mean, color='red', linestyle = 'dashed')
    plt.axvline(x=(-2*sd+mean), color='red', linestyle = 'dashed')

    outFile = re.sub('Data/Output/(.*).csv', r'Data/Figures/\1_histogram.png', inFile)
    plt.savefig(outFile)

'''
Function:
Purpose:
Parameters:
Returns:
'''
def draw_histogram_decoy(inFile):
    ppmSpreadFile = re.sub('(.*).csv', r'\1_ppmSpreadDecoy.csv', inFile)
    data = []
    with open(ppmSpreadFile, newline='') as f:
        read = csv.reader(f)
        for r in read: data.append(list([float(x) for x in r]))

    print(len(data[0]))
    print(len(data[1]))
    plt.figure(figsize=(8,6))
#    plt.hist(data, bins=200, stacked=True, density=True)
    plt.hist(data[0], bins=200)
    plt.hist(data[1], bins=200)
    plt.xlabel("PPM Difference", size=14)
    plt.ylabel("Count", size=14)
    plt.title("Histogram of PPM spread")
    outFile = re.sub('Data/Output/(.*).csv', r'Data/Figures/\1_histogram_decoys.png', inFile)
    plt.show()
    #plt.savefig(outFile)

def create_venn_diagrams():
    caleb_peptide_df = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-6peaks_exp-n1b_corrected_peptideFDR.csv')
    caleb_protein_df = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-6peaks_exp-n1b_corrected_proteinFDR.csv')
    jesse_peptide_df = pd.read_csv('Data/Input/peptide_matches_Jesse.csv')
    jesse_protein_df = pd.read_csv('Data/Input/protein_matches_Jesse.csv')

    caleb_peptides = sorted(list(caleb_peptide_df['peptide']))
    caleb_proteins = sorted(list(set(caleb_protein_df['leadingProtein'])),reverse=True)
    jesse_peptides = sorted(list(jesse_peptide_df['Peptide']))
    jesse_proteins = sorted(list(jesse_protein_df['Protein']))

    #print(caleb_peptides[:10])
    #print(jesse_peptides[:10])
    #print(caleb_proteins[:10])
    #print(jesse_proteins[:10])

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

    caleb_peptides = list(set(caleb_peptides))
    for i in range(len(jesse_peptides)):
        jesse_peptides[i] = re.sub('\+\d+\.\d+','',jesse_peptides[i])
    jesse_peptides = list(set(jesse_peptides))


    #print(len(caleb_proteins))
    #print(len(set(caleb_proteins)))
    caleb_proteins1 = [re.sub('(\d+/)(DECOY_0_)?(sp\|\w{6}\|)', r'\2\3', x) for x in caleb_proteins if (x[0]+x[1])=='1/']
    caleb_proteins = caleb_proteins1 + [x for x in caleb_proteins if (x[0]+x[1])!='1/']
    jesse_proteins = [re.sub('(.*)(DECOY_0_)?(sp\|\w{6}\|)(.*)', r'\2\3', x) for x in jesse_proteins]
    #for x in sorted(caleb_proteins): print(x)

    #print(len(caleb_proteins))
    #print(len(set(caleb_proteins)))
    #for x in sorted(caleb_peptides): print(x)
    #for i in range(10): print('-----------------')
    #for x in sorted(jesse_peptides): print(x)
    venn2([set(caleb_peptides),set(jesse_peptides)], set_labels = ["csoDIAq", "MSPLIT-DIA"])
    plt.title('Comparing Peptide Identification Outputs (stripped sequences, 6peaks)\n')
    plt.savefig('Data/Figures/peptide_comparison_venn.png')

    plt.clf()

    venn2([set(caleb_proteins),set(jesse_proteins)], set_labels = ["csoDIAq", "MSPLIT-DIA"])
    plt.title('Comparing Protein Identification Outputs (6peaks)\n')
    plt.savefig('Data/Figures/protein_comparison_venn.png')
