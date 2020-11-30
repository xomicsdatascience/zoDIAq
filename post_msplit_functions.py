import pandas as pd

overallDf = pd.read_csv('Data/fullOutput_lib31_match3_ppm10.csv').sort_values('cosine', ascending=False)
#peptideDf.to_csv('delete.csv')
#peptideDf = pd.read_csv('delete.csv')
overallDf = overallDf.reset_index(drop=True)

peptideDf = overallDf.drop_duplicates(subset='peptide', keep='first')
peptideDf = peptideDf.reset_index(drop=True)

proteinDf = overallDf.drop_duplicates(subset='protein', keep='first')
proteinDf = proteinDf.reset_index(drop=True)


def fdr_calculation(df, cutoff=0.01):
    fdr = []
    numDecoys = 0
    for i in range(len(df)):
        if 'DECOY' in df.loc[i]['protein']:
            numDecoys += 1
        new = numDecoys/(i+1)
        if new > cutoff:
            if len(fdr) < 1/cutoff:
                return []
            return fdr
        fdr.append(new)
    return fdr

def match_fdrSat_graph(df, fdrMax, outFile):
    matches = sorted(list(set(df['shared'])))

    fdrSat = []
    numPeaks = []
    cosine = []
    for m in matches:
        tempDf = df[df['shared'] >= m].reset_index(drop=True)
        fdr = fdr_calculation(tempDf)
        numPeaks.append(len(tempDf))
        fdrSat.append(len(fdr))
        if len(fdr) != 0: cos = df['cosine'].loc[len(fdr)-1]
        else: cos = 0
        cosine.append(cos)


    graphDf = pd.DataFrame({'matches':matches, 'FDRCutoff': fdrSat, 'total': numPeaks, 'cosine': cosine})
    graphDf.to_csv(outFile, index=False)

match_fdrSat_graph(overallDf, 0.01, 'Data/Figures/FDRGraph_lib31_overall.csv')
match_fdrSat_graph(peptideDf, 0.01, 'Data/Figures/FDRGraph_lib31_peptide.csv')
match_fdrSat_graph(proteinDf, 0.01, 'Data/Figures/FDRGraph_lib31_protein.csv')


#overallDf['FDR'] = fdr_calculation(overallDf)
#peptideDf['FDR'] = fdr_calculation(peptideDf)
#proteinDf['FDR'] = fdr_calculation(proteinDf)
