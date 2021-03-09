import csodiaq_base_functions as cbf
from timeit import default_timer as timer
from datetime import timedelta
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


'''
Function: write_csodiaq_output()
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_output(lib, expFile, outFile, initialTol=10, corrected=False, queryPooling=True ):
    print('enter spectra comparison:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
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
    if queryPooling:
        cbf.pooled_all_query_spectra_analysis( expFile,
                                    outFile,
                                    ppmFile,
                                    lib,
                                    ppmTol,
                                    ppmOffset)
    else:
        cbf.pooled_library_only_spectra_analysis( expFile,
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
    print('enter ppm offset and tolerance calculations:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)

    #if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    fields = ['scan','peptide','protein','zLIB','MaCC_Score']
    df = pd.read_csv(inFile, sep=',', usecols=fields).sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
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
    print('enter csodiaq FDR Calculation:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
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
def write_DISPA_targeted_reanalysis_files(inFile, proteins=1, trypsin=True, corrected=False):
    print('enter DISPA targeted reanalysis file writing:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
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
def heavy_light_quantification(inFile, libFile, expFiles, outDir, libPeaks, massTol, minMatches, ratioType, correction, hist):
    print('enter SILAC Quantification:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    fragDict, libDict = cbf.make_quant_dicts(inFile, libFile, expFiles, libPeaks)
    dfDIA = cbf.heavy_light_quantification(fragDict, libDict, expFiles, outDir, massTol, minMatches, ratioType, correction, hist)
    dfDIA.to_csv(outDir + 'CsoDIAq_output_SILAC_Quantification.csv')
