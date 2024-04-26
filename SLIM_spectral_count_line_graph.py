import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv('/Users/cranneyc/Desktop/topNPeakAnalysis.csv')
df = pd.read_csv('/Users/cranneyc/Desktop/topNPeakAnalysis_ppmCorrection.csv')


# Merge the four category columns into a single column
df['categories'] = df[['numMinutes', 'setting1', 'setting2', 'setting3']].astype(str).apply(lambda x: ','.join(x), axis=1)
tempNoPeakDf = df[df['topNPeaks']=='noTopPeakFilter']
tempPeakDf = df[df['topNPeaks']!='noTopPeakFilter']
df = pd.concat([tempPeakDf, tempNoPeakDf])


# Plotting
plt.figure(figsize=(10, 6))
for category, group in df.groupby('categories'):
    plt.plot(group['topNPeaks'], group['spectralCountBelowFDR'], marker='o', label=category)
    #plt.plot(group['topNPeaks'], group['peptideCountBelowFDR'], marker='o', label=category)

plt.xlabel('top N Peaks (ms1 filtering)')
plt.ylabel('Count')
plt.title('Spectral counts below 0.01 FDR for different settings')
#plt.title('Peptide counts below 0.01 FDR for different settings')
plt.legend()
plt.show()
