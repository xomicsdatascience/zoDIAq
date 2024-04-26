import pandas as pd
import os
from collections import Counter
import re

parentDir = '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM_ms1_filtering_with_ppm_correction'
dirs = sorted([os.path.join(parentDir, x, 'fdrScores-macc-maxlfq') for x in os.listdir(parentDir) if x != '.DS_Store'])
dirs += ['/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM_HeLa-zodiaq-id-20240319-075759/fdrScores-macc-maxlfq']

filePattern = r'zoDIAq-file_400ng-uL_HeLa_SLIM-DIA_FT_2080ms_2.04s_MS2_30k_NCE35_([\d\.]+)min_([\d\.]+)_([\d\.]+)_([\d\.]+)_DIA_corrected_fullOutput_\w+.csv'
dirPattern = r'top(\d+)Peaks-zodiaq-id-\d+-\d+/fdrScores-macc-maxlfq'

data = []
for dir in dirs:
    allFiles = os.listdir(dir)
    spectralFiles = sorted([x for x in allFiles if 'spectralFDR' in x])
    peptideFiles = sorted([x for x in allFiles if 'peptideFDR' in x])
    for i in range(len(spectralFiles)):
        spectralDf = pd.read_csv(os.path.join(dir, spectralFiles[i]))
        peptideDf = pd.read_csv(os.path.join(dir, peptideFiles[i]))
        spectralDf = spectralDf[spectralDf['spectralFDR'] < 0.01]
        peptideDf = peptideDf[peptideDf['peptideFDR'] < 0.01]
        settings = [float(x) for x in re.findall(filePattern, spectralFiles[i])[0]]
        topNPeaks = re.findall(dirPattern, dir)
        if len(topNPeaks):
            topNPeaks = int(topNPeaks[0])
        else:
            topNPeaks = 'noTopPeakFilter'
        data.append([topNPeaks, len(spectralDf), len(peptideDf)] + settings + [spectralFiles[i], dir])

df = pd.DataFrame(data, columns=['topNPeaks', 'spectralCountBelowFDR', 'peptideCountBelowFDR', 'numMinutes', 'setting1', 'setting2', 'setting3', 'file', 'dir'])
#df.to_csv('/Users/cranneyc/Desktop/topNPeakAnalysis_ppmCorrection.csv', index=False)
#df.to_csv('/Users/cranneyc/Desktop/topNPeakAnalysis_ppmCorrection.csv', index=False)
print(df)




