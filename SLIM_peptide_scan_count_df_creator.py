import pandas as pd
import os
from collections import Counter
import re

dirs = [
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM_HeLa-zodiaq-id-20240319-075759/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top010ms1Peaks-zodiaq-id-20240419-130648/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top025ms1Peaks-zodiaq-id-20240419-130057/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top050ms1Peaks-zodiaq-id-20240419-074339/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top100ms1Peaks-zodiaq-id-20240419-072557/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top200ms1Peaks-zodiaq-id-20240419-074240/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top300ms1Peaks-zodiaq-id-20240419-075012/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top400ms1Peaks-zodiaq-id-20240419-075115/fdrScores-macc-maxlfq',
    '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/delete_SLIM_library_filtering_experiment_top600ms1Peaks-zodiaq-id-20240419-075430/fdrScores-macc-maxlfq',
]
filePattern = r'zoDIAq-file_400ng-uL_HeLa_SLIM-DIA_FT_2080ms_2.04s_MS2_30k_NCE35_([\d\.]+)min_([\d\.]+)_([\d\.]+)_([\d\.]+)_DIA_corrected_fullOutput_\w+.csv'
dirPattern = r'delete_SLIM_library_filtering_experiment_top(\d+)ms1Peaks-zodiaq-id-\d+-\d+/fdrScores-macc-maxlfq'

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
df.to_csv('/Users/cranneyc/Desktop/topNPeakAnalysis.csv', index=False)
print(df)




