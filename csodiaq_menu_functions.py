import csodiaq_base_functions as cbf
from timeit import default_timer as timer
from datetime import timedelta
import re
import pandas as pd
import os
from statistics import pstdev, median
import csv
import matplotlib.pyplot as plt



'''
Function: write_csodiaq_output()
Purpose:
Parameters:
Returns:
'''
def write_csodiaq_output(lib, expFile, outFile, corrected=False, numLibPeaks=31 ):
    print('#Enter lib upload/conversion:')
    print('#'+str(timedelta(seconds=timer())))
    if corrected:
        offsetFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', outFile)
        df = pd.read_csv(offsetFile)
        ppmTol = df['tolerance'].loc[0]
        ppmOffset=(-df['offset'].loc[0])
        outFile = re.sub('(.*).csv', r'\1_corrected.csv', outFile)
    else:
        ppmTol=10,
        ppmOffset=0
    if numLibPeaks > 31: print("Number of library peaks cannot exceed 30, reducing to 30"); numLibPeaks = 30
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
#    lib = cbf.library_file_to_dict(libFile, numLibPeaks)
    print('#enter spectra comparison:')
    print('#'+str(timedelta(seconds=timer())))
    cbf.query_spectra_analysis( expFile,
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
    df = pd.read_csv(inFile).sort_values('csoDIAq_Score', ascending=False).reset_index(drop=True)
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

def write_ppm_spread(inFile, corrected=0):
    print('#enter csodiaq ppm spread file creation:')
    print('#'+str(timedelta(seconds=timer())))
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)


    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    outFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    tag = re.sub('Data/Output/(.*).csv', r'\1_filteredBestMatch', inFile)
    files = list(os.listdir('Data/Output'))
    for x in files:
        if tag in x: inFile = 'Data/Output/'+x
    cbf.write_ppm_spread(inFile, ppmFile, outFile)
    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))
'''
'''
Function:
Purpose:
Parameters:
Returns:

def write_ppm_offset_tolerance(inFile, corrected=0, hist=False):
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    ppmSpreadFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    with open(ppmSpreadFile, newline='') as f: ppmList = list(csv.reader(f))[0]
    ppmList = [float(x) for x in ppmList]
    outFile = re.sub('(.*).csv', r'\1_offset_tolerance.csv', inFile)
    if hist: histFile = re.sub('Data/Output/(.*).csv', r'Data/Figures/\1_histogram.png', inFile)
    else: histFile = 0

    offset, tolerance = cbf.find_offset_tol(ppmList, histFile)
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(['offset','tolerance'])
        writer.writerow([offset, tolerance])
'''
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
