import csodiaq_base_functions as cbf
from timeit import default_timer as timer
from datetime import timedelta
import re
import pandas as pd
from os import listdir
from os.path import isfile, join
from statistics import pstdev, median
import csv
import matplotlib.pyplot as plt
import numpy as np
import statistics
from scipy.stats import linregress


'''
Function: write_csodiaq_output()
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_output(lib, expFile, outFile, corrected=False ):
    if corrected:
        offsetFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', outFile)
        df = pd.read_csv(offsetFile)
        ppmTol = df['tolerance'].loc[0]
        ppmOffset=(-df['offset'].loc[0])
        outFile = re.sub('(.*).csv', r'\1_corrected.csv', outFile)
    else:
        ppmTol=10,
        ppmOffset=0
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
#    lib = cbf.library_file_to_dict(libFile, numLibPeaks)
    print('#enter spectra comparison:')
    print('#'+str(timedelta(seconds=timer())))
    cbf.pooled_all_query_spectra_analysis( expFile,
#    cbf.query_spectra_analysis( expFile,
                                outFile,
                                ppmFile,
                                lib,
                                ppmTol,
                                ppmYOffset=ppmOffset)

    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_ppm_offset_tolerance(inFile, corrected=0, hist=False):

    print('#enter ppm offset and tolerance calculations:')
    print('#'+str(timedelta(seconds=timer())))

    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    df = pd.read_csv(inFile).sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
    hits, decoys = cbf.fdr_calculation(df)
    df = df.loc[:hits]

    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    ppmList = cbf.return_ppm_spread(df, ppmFile)

    outFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', inFile)
    if hist: histFile = re.sub('Data/Output/(.*).csv', r'Data/Figures/\1_histogram.png', inFile)
    else: histFile = 0

    offset, tolerance = cbf.find_offset_tol(ppmList, histFile)
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(['offset','tolerance'])
        writer.writerow([offset, tolerance])
    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))


'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_fdr_outputs(inFile, corrected=False):
    print('#enter csodiaq FDR Calculation:')
    print('#'+str(timedelta(seconds=timer())))
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    spectralFile = re.sub('(.*).csv', r'\1_spectralFDR.csv', inFile)
    peptideFile = re.sub('(.*).csv', r'\1_peptideFDR.csv', inFile)
    proteinFile = re.sub('(.*).csv', r'\1_proteinFDR.csv', inFile)

    cbf.write_csodiaq_fdr_outputs(inFile, spectralFile, peptideFile, proteinFile)

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_DISPA_targeted_reanalysis_files(inFile, proteins=1, trypsin=True):
    print('#enter DISPA targeted reanalysis file writing:')
    print('#'+str(timedelta(seconds=timer())))
    inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    header = re.sub('(.*).csv', r'\1_mostIntenseTargs_', inFile)
    if proteins: inFile = re.sub('(.*).csv', r'\1_proteinFDR.csv', inFile)
    else: inFile = re.sub('(.*).csv', r'\1_peptideFDR.csv', inFile)
    cbf.return_DISPA_targeted_reanalysis_dfs(header, inFile, proteins, trypsin)


'''
Function:
Purpose:
Parameters:
Returns:
'''
def quantify(inFile, libFile, expFolder):
    inFile = re.sub('(.*).csv', r'\1_corrected_proteinFDR.csv', inFile)
    #files = [expFolder+f for f in listdir(expFolder) if isfile(join(expFolder, f))]
    files = [
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to8_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to4_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to2_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to1_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_2to1_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_4to1_01.mzXML',
        'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_8to1_01.mzXML'
    ]
    ratios = [1, 2, 4, 8]
    lc = calc_all_variables(pd.read_csv('Data/MQpeptides_Quant.csv'), 'LC', ratios)
    X = lc['median']

    print('#LC average standard deviation: '+str(np.mean(lc['stdev'])))
    print('#LC total peptides quantified: '+str(sum(lc['numPeps'])))
    data = []
    '''
    var = [1, 'median', 1]
    print(var)
    dfDIA = cbf.quantify(inFile, libFile, files, var)
    dia = calc_all_variables(dfDIA, 'DIA', ratios)
    Y = dia['median']
    data.append([3]+var+list(linregress(X, Y))+[np.mean(dia['stdev']),sum(dia['numPeps'])])
    print(data[-1])
    finalDf = pd.DataFrame(data, columns=['libPeaks', 'minMatch', 'mode', 'corrStDev', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr', 'avrStDev', 'numPeps'])
    finalDf.to_csv('Data/QuantifyCompare/variables/compare.csv', index=False)
    '''

    #minMatch = [1, 2]
    #minMatch = [1, 2, 3, 4]
    #minMatch = [1, 2, 3, 5, 7, 9]
    #m = ['mean', 'median', 'intensity', 'weighted']
    #stdev = [1, 2]
    minMatch = [1]
    m = ['median']
    stdev = [1]
    for x in minMatch:
        for y in m:
            for z in stdev:
                var = [x, y, z]
                #print(var)
                dfDIA = cbf.quantify(inFile, libFile, files, var)
                temp = dfDIA.drop('peptide', axis=1)
                temp.to_csv('Data/CalebOutput.csv')
                plt.clf()
                temp.boxplot(column=files)
                plt.axhline(y=-3, color='blue', linestyle='dashed')
                plt.axhline(y=-2, color='blue', linestyle='dashed')
                plt.axhline(y=-1, color='blue', linestyle='dashed')
                plt.axhline(y=0, color='blue', linestyle='dashed')
                plt.axhline(y=1, color='blue', linestyle='dashed')
                plt.axhline(y=2, color='blue', linestyle='dashed')
                plt.axhline(y=3, color='blue', linestyle='dashed')
                #plt.ylim(-10, 10)
                plt.show()




                dia = calc_all_variables(dfDIA, 'DIA', ratios)
                print(dia['numPeps'])
                dia.to_csv('Data/lib3-1-median-1.csv', index=False)
                Y = dia['median']
                data.append([3]+var+list(linregress(X, Y))+[np.mean(dia['stdev']),sum(dia['numPeps'])])
                print(data[-1])
    finalDf = pd.DataFrame(data, columns=['libPeaks', 'minMatch', 'mode', 'corrStDev', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr', 'avrStDev', 'numPeps'])
    finalDf.to_csv('Data/QuantifyCompare/variables/compare.csv', index=False)

def calc_key_variables(df, r, type, ori=1):
    if ori:
        if type=='LC': col = 'Ratio H/L 1to'+str(r)
        else: col = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to'+str(r)+'_01.mzXML'
        ratio = r
    else:
        if type=='LC': col = 'Ratio H/L '+str(r)+'to1'
        else: col = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_'+str(r)+'to1_01.mzXML'
        ratio = -r
    if r==1: ratio=0


    if type=='LC': df[col] = np.log2(df[col])
    app = df[~df[col].isnull()][col]
    if len(app) > 0: return [[np.median(app),statistics.pstdev(app),ratio, type, len(app)]]

    return [[0,0, ratio, type, 0]]

def calc_all_variables(df, type, ratios):
    data = []
    for ratio in ratios:
        data += calc_key_variables(df, ratio, type)
        if ratio != 1: data += calc_key_variables(df, ratio, type, ori=0)
    finalDf = pd.DataFrame(data, columns=['median','stdev','ratio','type', 'numPeps'])
    finalDf = finalDf.sort_values('ratio', ascending=False).reset_index(drop=True)
    return finalDf
