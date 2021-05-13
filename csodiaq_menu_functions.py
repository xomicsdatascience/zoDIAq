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
Function:
Purpose:
Parameters:
Returns:
'''
def write_DISPA_targeted_reanalysis_files(inFile, proteins, heavy):
    print('enter DISPA targeted reanalysis file writing:')
    print(str(timedelta(seconds=timer())), flush=True)
    header = re.sub('(.*).csv', r'\1_mostIntenseTargs', inFile)
    if proteins: inFile = re.sub('(.*).csv', r'\1_proteinFDR.csv', inFile)
    else: inFile = re.sub('(.*).csv', r'\1_peptideFDR.csv', inFile)
    cbf.return_DISPA_targeted_reanalysis_dfs(header, inFile, proteins, heavy)

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
