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
def write_csodiaq_output(libFile, expFile, outFile, corrected=False ):
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
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
    lib = cbf.library_file_to_dict(libFile)
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
def filter_optimal_match_csodiaq_output(inFile, corrected=0):
    print('#enter csodiaq output filtering:')
    print('#'+str(timedelta(seconds=timer())))

    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    df = pd.read_csv(inFile).sort_values('cosine', ascending=False).reset_index(drop=True)
    bestMatchNum, bestFDR = cbf.find_best_matchNum_fdr(df, 0.01)

    outFile = re.sub('(.*).csv', r'\1_filteredBestMatch'+str(bestMatchNum)+'.csv', inFile)
    df.iloc[:bestFDR].to_csv(outFile, index=False)

    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))
    return outFile

'''
Function:
Purpose:
Parameters:
Returns:
'''
def write_ppm_spread(inFile, corrected=0):
    print('#enter csodiaq ppm spread file creation:')
    print('#'+str(timedelta(seconds=timer())))
    if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)

    tag = re.sub('Data/Output/(.*).csv', r'\1', inFile)
    r = re.compile(tag+'_filteredBestMatch')
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    outFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    cbf.write_ppm_spread(inFile, ppmFile, outFile)
    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))

'''
Function:
Purpose:
Parameters:
Returns:
'''
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

def write_ppm_spread_decoy(inFile):
    print('#enter csodiaq ppm spread file creation:')
    print('#'+str(timedelta(seconds=timer())))

    tag = re.sub('Data/Output/(.*).csv', r'\1', inFile)
    r = re.compile(tag+'_filteredBestMatch')
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    decoyFile = re.sub('(.*).csv', r'\1_ppmSpreadDecoy.csv', inFile)

    cbf.write_ppm_spread_decoy(inFile, ppmFile, decoyFile)
    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))
