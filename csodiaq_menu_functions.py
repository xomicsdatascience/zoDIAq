import csodiaq_base_functions as cbf
from timeit import default_timer as timer
from datetime import timedelta
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import linecache
import pickle

'''
Function: write_csodiaq_output()
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_output(lib, expFile, outFile, initialTol=10, corrected=False, queryPooling=0 ):
    print('enter spectra comparison:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    if corrected:
        offsetFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', outFile)
        df = pd.read_csv(offsetFile)
        ppmTol = df['tolerance'].loc[0]
        ppmOffset=(df['offset'].loc[0])
        outFile = re.sub('(.*).csv', r'\1_corrected.csv', outFile)
    else:
        ppmTol=initialTol
        ppmOffset=0
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
    if not queryPooling: queryPooling = np.inf
    initialOutFile = re.sub('(.*).csv', r'\1_delete.csv', outFile)
    cbf.pooled_spectra_analysis(  expFile,
                                        initialOutFile,
                                        lib,
                                        ppmTol,
                                        ppmOffset,
                                        queryPooling)

    spectraKeys = cbf.generate_valid_FDR_spectra_keys(initialOutFile)
    print('spectraKeys: ' + str(len(spectraKeys)))
    pickle.dump(spectraKeys, open('C:/Users/ccranney/Desktop/Caleb_Files/data/output/spectraKeys.p','wb'))

    ppmList = cbf.pooled_spectra_analysis(  expFile,
                                        outFile,
                                        lib,
                                        ppmTol,
                                        ppmOffset,
                                        queryPooling,
                                        spectraKeys=spectraKeys)

    return ppmList

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_ppm_offset_tolerance(inFile, ppmList, corrected=0, hist=False):
    print('enter ppm offset and tolerance calculations:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    '''
    #if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    df = pd.read_csv(inFile, sep=',').sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
    hits, decoys = cbf.fdr_calculation(df)
    df = df.loc[:hits]
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    ppmList = cbf.return_ppm_spread(df, ppmFile)
    '''
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
def write_csodiaq_fdr_outputs(inFile, proteins, corrected=False):
    print('enter csodiaq FDR Calculation:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    spectralFile = re.sub('(.*).csv', r'\1_spectralFDR.csv', inFile)
    peptideFile = re.sub('(.*).csv', r'\1_peptideFDR.csv', inFile)
    if proteins: proteinFile = re.sub('(.*).csv', r'\1_proteinFDR.csv', inFile)
    else: proteinFile = ''

    cbf.write_csodiaq_fdr_outputs(inFile, spectralFile, peptideFile, proteinFile)

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_DISPA_targeted_reanalysis_files(inFile, proteins=1, trypsin=True, corrected=False, heavy=False):
    print('enter DISPA targeted reanalysis file writing:')
    print(str(timedelta(seconds=timer())), flush=True)
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    header = re.sub('(.*).csv', r'\1_mostIntenseTargs', inFile)
    if proteins: inFile = re.sub('(.*).csv', r'\1_proteinFDR.csv', inFile)
    else: inFile = re.sub('(.*).csv', r'\1_peptideFDR.csv', inFile)
    cbf.return_DISPA_targeted_reanalysis_dfs(header, inFile, proteins, trypsin, heavy)

    print('Complete:')
    print(str(timedelta(seconds=timer())), flush=True)
'''
Function:
Purpose:
Parameters:
Returns:
'''
def heavy_light_quantification(inFile, libFile, expFiles, outDir, libPeaks, massTol, minMatches, ratioType, correction, hist):
    print('enter SILAC Quantification:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    fragDict, libDict = cbf.make_quant_dicts(inFile, libFile, expFiles, libPeaks)
    dfDIA = cbf.heavy_light_quantification(fragDict, libDict, expFiles, outDir, massTol, minMatches, ratioType, correction, hist)
    dfDIA.to_csv(outDir + 'CsoDIAq_output_SILAC_Quantification.csv')
