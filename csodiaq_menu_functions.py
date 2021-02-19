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
#from scipy.stats import linregress
import pickle


'''
Function: write_csodiaq_output()
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_output(lib, expFile, outFile, initialTol=10, corrected=False ):
    print('#enter spectra comparison:', flush=True)
    print('#'+str(timedelta(seconds=timer())), flush=True)
    if corrected:
        offsetFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', outFile)
        df = pd.read_csv(offsetFile)
        ppmTol = df['tolerance'].loc[0]
        ppmOffset=(df['offset'].loc[0])
        outFile = re.sub('(.*).csv', r'\1_corrected.csv', outFile)
    else:
        ppmTol=initialTol,
        ppmOffset=0
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
    cbf.pooled_all_query_spectra_analysis( expFile,
                                outFile,
                                ppmFile,
                                lib,
                                ppmTol,
                                ppmOffset)


'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_ppm_offset_tolerance(inFile, corrected=0, hist=False):
    print('#enter ppm offset and tolerance calculations:', flush=True)
    print('#'+str(timedelta(seconds=timer())), flush=True)

    #if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    df = pd.read_csv(inFile, sep=',').sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
    hits, decoys = cbf.fdr_calculation(df)
    df = df.loc[:hits]

    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    ppmList = cbf.return_ppm_spread(df, ppmFile)

    outFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', inFile)
    if hist: histFile = re.sub('(.*).csv', r'\1_histogram.png', inFile)
    else: histFile = 0

    offset, tolerance = cbf.find_offset_tol(ppmList, histFile, stdev=corrected)
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
    print('#enter csodiaq FDR Calculation:', flush=True)
    print('#'+str(timedelta(seconds=timer())), flush=True)
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
    print('#enter DISPA targeted reanalysis file writing:', flush=True)
    print('#'+str(timedelta(seconds=timer())), flush=True)
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
def heavy_light_quantification(inFile, libFile, expFolder):
    print('#enter SILAC Quantification:', flush=True)
    print('#'+str(timedelta(seconds=timer())), flush=True)
    inFile = re.sub('(.*).csv', r'\1_corrected_mostIntenseTargs_allCVs.csv', inFile)
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
    dfLC = pd.read_csv('Data/MQpeptides_Quant.csv')
    lc = calc_all_variables(dfLC, 'LC', ratios)
    X = lc['median']

    lcCols = [
        'Ratio H/L 1to8',
        'Ratio H/L 1to4',
        'Ratio H/L 1to2',
        'Ratio H/L 1to1',
        'Ratio H/L 2to1',
        'Ratio H/L 4to1',
        'Ratio H/L 8to1',
    ]

    ratioStr = [
        '1:8',
        '1:4',
        '1:2',
        '1:1',
        '2:1',
        '4:1',
        '8:1',
    ]

    boxLC = dfLC.loc[:, dfLC.columns.intersection(lcCols)]
    boxLC = boxLC[lcCols]
    boxLC.columns = ratioStr
    plt.clf()
    set_plot_settings('','')
    boxLC.boxplot(column=ratioStr)
    plt.axhline(y=-3, color='blue', linestyle='dashed')
    plt.axhline(y=-2, color='blue', linestyle='dashed')
    plt.axhline(y=-1, color='blue', linestyle='dashed')
    plt.axhline(y=0, color='blue', linestyle='dashed')
    plt.axhline(y=1, color='blue', linestyle='dashed')
    plt.axhline(y=2, color='blue', linestyle='dashed')
    plt.axhline(y=3, color='blue', linestyle='dashed')
    plt.ylim(-6.2, 6.2)
    plt.savefig('Data/QuantifyCompare/BoxPlots/LCOriginal.png')

    plt.clf()
    set_plot_settings('','')
    finalDf = pd.read_csv('Data/jesseOutput.csv')
    finalDf.columns = ['scan'] + ratioStr
    finalDf.boxplot(column=ratioStr)
    plt.axhline(y=-3, color='blue', linestyle='dashed')
    plt.axhline(y=-2, color='blue', linestyle='dashed')
    plt.axhline(y=-1, color='blue', linestyle='dashed')
    plt.axhline(y=0, color='blue', linestyle='dashed')
    plt.axhline(y=1, color='blue', linestyle='dashed')
    plt.axhline(y=2, color='blue', linestyle='dashed')
    plt.axhline(y=3, color='blue', linestyle='dashed')
    plt.ylim(-6.2, 6.2)

    plt.savefig('Data/QuantifyCompare/BoxPlots/jesseOriginal.png')


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
    graph = True
    #libPeaks = [3, 5, 10, 'all']
    libPeaks = ['all']
    #minMatches = [[1, 2],[0, 1, 3, 4],[0, 1, 2, 3, 5],[0, 1, 2, 3, 5]]
    minMatches = [[0]]
    #m = ['mean', 'median', 'intensity', 'weighted']
    m = ['median']
    #stdev = [-1, 0, 0.5, 1]
    stdev = [0]

    for w in range(len(libPeaks)):
        #fragDict, libDict = cbf.make_quant_dicts('Data/Input/TempHold/mostintense_quantmzlist.txt', 'Data/Input/human.faims.fixed.decoy.mgf', files, libPeaks[w])
        #fragDict, libDict = cbf.make_quant_dicts(inFile, libFile, files, libPeaks[w])
        #fragDict, libDict = cbf.make_quant_dicts('Data/Input/TempHold/mostintense_quantmzlist.txt', libFile, files, libPeaks[w])
        #pickle.dump(libDict, open( "Data/Input/TempHold/lib_"+str(libPeaks[w])+"_traml.p", "wb" ))
        #pickle.dump(libDict, open( "Data/Input/TempHold/lib_"+str(libPeaks[w])+"_mgf.p", "wb" ))
        #pickle.dump(fragDict, open( "Data/Input/TempHold/frag_mgf.p", "wb" ))
        #pickle.dump(fragDict, open( "Data/Input/TempHold/frag_traml.p", "wb" ))
        minMatch = minMatches[w]
        libDict = pickle.load(open( "Data/Input/TempHold/lib_"+str(libPeaks[w])+"_mgf.p", "rb" ))
        #libDict = pickle.load(open( "Data/Input/TempHold/lib_"+str(libPeaks[w])+"_traml.p", "rb" ))
        fragDict = pickle.load(open( "Data/Input/TempHold/frag_mgf.p", "rb" ))
        #fragDict = pickle.load(open( "Data/Input/TempHold/frag_traml.p", "rb" ))
        #print(len(fragDict))
        for x in minMatch:
            for y in m:
                for z in stdev:
                    var = [x, y, z]
                    strVar = '_'.join([str(libPeaks[w])]+[str(v) for v in var])

                    #print(var)
                    dfDIA = cbf.heavy_light_quantification(fragDict, libDict, files, var)

                    if graph:
                        temp = dfDIA.drop('peptide', axis=1)
                        temp.to_csv('Data/QuantifyCompare/RatioFiles/lib'+strVar+'.csv', index=False)
                        plt.clf()
                        set_plot_settings('','')

                        temp.columns = ['scan']+ratioStr
                        temp.boxplot(column=ratioStr)
                        plt.axhline(y=-3, color='blue', linestyle='dashed')
                        plt.axhline(y=-2, color='blue', linestyle='dashed')
                        plt.axhline(y=-1, color='blue', linestyle='dashed')
                        plt.axhline(y=0, color='blue', linestyle='dashed')
                        plt.axhline(y=1, color='blue', linestyle='dashed')
                        plt.axhline(y=2, color='blue', linestyle='dashed')
                        plt.axhline(y=3, color='blue', linestyle='dashed')
                        plt.ylim(-6.2, 6.2)
                        plt.savefig('Data/QuantifyCompare/BoxPlots/lib'+strVar+'.png')




                    dia = calc_all_variables(dfDIA, 'DIA', ratios)
                    print(dia['numPeps'])
                    dia.to_csv('Data/QuantifyCompare/RatioCompareFiles/lib'+strVar+'.csv', index=False)
                    Y = dia['median']
                    data.append([libPeaks[w]]+var+list(linregress(X, Y))+[np.mean(dia['stdev']),sum(dia['numPeps'])])
                    lineDf = pd.concat([lc, dia])
                    ratioDict = {
                            -8:-3,
                            -4:-2,
                            -2:-1,
                            0:0,
                            2:1,
                            4:2,
                            8:3
                        }
                    lineDf['ratio'] = [ratioDict[x] for x in lineDf['ratio']]


                    plt.clf()
                    fig, ax = plt.subplots()
                    for key, group in lineDf.groupby('type'):
                        group.plot('ratio', 'median', yerr='stdev',
                                   label=key, ax=ax)
                    plt.ylim(-4.5, 4.5)
                    plt.savefig('Data/QuantifyCompare/LineGraphs/lib'+strVar+'.png')
                    plt.close("all")
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

def set_plot_settings(xlabel, ylabel, wide=True):
    #if wide: plt.figure(figsize=(18,12))
    #else: plt.figure(figsize=(12,12))
    #plt.axhline(linewidth=4, y=1, color='black')
    #plt.axhline(linewidth=4, y=10, color='black')
    #plt.axvline(linewidth=4, x=1.5, color='black')
    #plt.axvline(linewidth=4, x=10, color='black')
    #plt.xlabel(xlabel, fontsize = 36, weight='bold')
    #plt.ylabel(ylabel, fontsize = 36, weight='bold')
    plt.tick_params(axis="x", labelsize=18)
    plt.tick_params(axis="y", labelsize=18)
