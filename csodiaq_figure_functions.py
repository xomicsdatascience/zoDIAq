import pandas as pd
import re
import csodiaq_base_functions as cbf
import csv
from statistics import pstdev, median
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_meta_analysis_files(inFile):
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
        tempDf = df[df['shared'] >= m].reset_index(drop=True)

        if 'peptide' in outFile: tempDf = tempDf.drop_duplicates(subset='peptide', keep='first'); tempDf = tempDf.reset_index(drop=True)
        elif 'protein' in outFile: tempDf = tempDf.drop_duplicates(subset='protein', keep='first'); tempDf = tempDf.reset_index(drop=True)

        fdrRows, dec = cbf.fdr_calculation(tempDf, 0.01)

        numPeaks.append(len(tempDf))
        fdrSat.append(fdrRows)
        subsetDecoys.append(dec)
        totalDecoys.append(len( [ x for x in list(tempDf['protein']) if 'DECOY' in x ] ))
        if fdrRows != 0: cos = tempDf['cosine'].loc[fdrRows-1]
        else: cos = 0
        cosine.append(cos)

    graphDf = pd.DataFrame({'matches':matches, 'FDRCutoff': fdrSat, 'total': numPeaks, 'cosine': cosine, 'FDRDecoys': subsetDecoys, 'totalDecoys': totalDecoys})
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
