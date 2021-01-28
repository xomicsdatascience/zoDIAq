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
    files = [expFolder+f for f in listdir(expFolder) if isfile(join(expFolder, f))]
    dfLC = pd.read_csv('Data/MQpeptides_Quant.csv')


    var = [5, 'median', 1]
    dfDIA = cbf.quantify(inFile, libFile, files, var)

    

    '''
    minMatch = [3]
    m = ['mean', 'median']
    stdev = [0.5, 1, 2]
    for x in minMatch:
        for y in m:
            for z in stdev:
                var = [x, y, x-1, z]
                cbf.quantify(inFile, libFile, files, var)
                if x > 2:
                    var = [x, y, x//2, z]
                    cbf.quantify(inFile, libFile, files, var)
    '''
