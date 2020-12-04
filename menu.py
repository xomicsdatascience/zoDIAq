import ppm_optimization_step as pos
from timeit import default_timer as timer
from datetime import timedelta
import re
import pandas as pd
import os
from statistics import pstdev
import csv
import matplotlib.pyplot as plt


def write_csodiaq_output(libFile, expFile, outFile, matches, correctedNumSD=0 ):
    print('#Enter lib upload/conversion:')
    print('#'+str(timedelta(seconds=timer())))
    if correctedNumSD:
        offsetFile = re.sub('(.*).csv', r'\1_offset_SD.csv', outFile)
        df = pd.read_csv(offsetFile)
        ppmTol = correctedNumSD*df['standardDeviation'].loc[0]
        ppmOffset=(-df['offset'].loc[0])
        outFile = re.sub('(.*).csv', r'\1_corrected'+str(correctedNumSD)+'.csv', outFile)
    else:
        ppmTol=10,
        ppmOffset=0
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', outFile)
    lib = pos.traml_library_upload_csv(libFile)
    print('#enter spectra comparison:')
    print('#'+str(timedelta(seconds=timer())))
    pos.query_spectra_analysis( expFile,
                                outFile,
                                ppmFile,
                                lib,
                                matches,
                                ppmTol,
                                ppmYOffset=ppmOffset)

    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))

def filter_optimal_match_csodiaq_output(inFile, correctedNumSD=0):
    print('#enter csodiaq output filtering:')
    print('#'+str(timedelta(seconds=timer())))

    if correctedNumSD: inFile = re.sub('(.*).csv', r'\1_corrected'+str(correctedNumSD)+'.csv', inFile)
    df = pd.read_csv(inFile).sort_values('cosine', ascending=False).reset_index(drop=True)
    bestMatchNum, bestFDR = pos.find_best_matchNum_fdr(df, 0.01)

    outFile = re.sub('(.*).csv', r'\1_filteredBestMatch'+str(bestMatchNum)+'.csv', inFile)
    df.iloc[:bestFDR].to_csv(outFile, index=False)

    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))
    return outFile

def write_ppm_spread(inFile, correctedNumSD=0):
    print('#enter csodiaq ppm spread file creation:')
    print('#'+str(timedelta(seconds=timer())))
    if correctedNumSD: inFile = re.sub('(.*).csv', r'\1_corrected'+str(correctedNumSD)+'.csv', inFile)

    tag = re.sub('Data/Output/(.*).csv', r'\1', inFile)
    r = re.compile(tag+'_filteredBestMatch')
    filterFile = 'Data/Output/'+list(filter(r.match,os.listdir('Data/Output/')))[0]
    ppmFile = re.sub('(.*).csv', r'\1_unfilteredPpmPerRow.csv', inFile)
    outFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    pos.write_ppm_spread(inFile, ppmFile, outFile)
    print('#Complete')
    print('#'+str(timedelta(seconds=timer())))

def write_ppm_offset_sd(inFile, correctedNumSD=0):
    if correctedNumSD: inFile = re.sub('(.*).csv', r'\1_corrected'+str(correctedNumSD)+'.csv', inFile)
    ppmSpreadFile = re.sub('(.*).csv', r'\1_ppmSpread.csv', inFile)
    with open(ppmSpreadFile, newline='') as f: ppmList = list(csv.reader(f))[0]
    ppmList = [float(x) for x in ppmList]
    outFile = re.sub('(.*).csv', r'\1_offset_SD.csv', inFile)
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(['offset','standardDeviation'])
        writer.writerow([(sum(ppmList)/len(ppmList)), pstdev(ppmList)])

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

    outFile = re.sub('(.*).csv', r'\1_histogram.png', inFile)
    plt.savefig(outFile)
