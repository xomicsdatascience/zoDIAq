import pandas as pd
import os
from collections import Counter

#dir = '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM_summed_matchDf-zodiaq-id-20240119-075329/fdrScores-macc-maxlfq'
#dir = '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM_summedQueryPeaks-zodiaq-id-20240119-111811/fdrScores-macc-maxlfq'
dir = '/Users/cranneyc/Documents/Projects/csodiaqOutput/currentDevBranchOutput/SLIM-zodiaq-id-20240102-194549/fdrScores-macc-maxlfq'
allFiles = os.listdir(dir)
spectralFiles = sorted([x for x in allFiles if 'spectralFDR' in x])
allPeptides = set()
allScans = set()

for file in spectralFiles:
    df = pd.read_csv(os.path.join(dir, file))
    allPeptides = allPeptides.union(df['peptide'])
    allScans = allScans.union(df['scan'])

peptideDf = pd.DataFrame(index=sorted(allPeptides), columns=sorted(spectralFiles))
scanDf = pd.DataFrame(index=sorted(allScans), columns=sorted(spectralFiles))

for file in spectralFiles:
    df = pd.read_csv(os.path.join(dir, file))
    peptideCount = Counter(df['peptide'])
    scanCount = Counter(df['scan'])
    peptideDf[file] = [peptideCount.get(item, 0) for item in peptideDf.index]
    scanDf[file] = [scanCount.get(item, 0) for item in scanDf.index]
    #df['peptideCount'] = [peptideCount[peptide] for peptide in df['peptide']]
    #print(file)
    #print(df)
    #df.to_csv(os.path.join(dir, file), index=False)

peptideDf.to_csv(os.path.join(dir, '00peptide_count_in_each_file.csv'))
scanDf.to_csv(os.path.join(dir, '00num_peptides_per_scan_in_each_file.csv'))



