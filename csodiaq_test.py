import pandas as pd
import os
import csodiaq_base_functions as cbf
import re
from collections import defaultdict
import statistics
import csv
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import numpy as np

def histogram(df, fileHeader, col):
    target_scores = []
    decoy_scores = []

    temp = df.sort_values(col, ascending=False)
    temp = temp.drop_duplicates(subset='peptide', keep='first').reset_index(drop=True)

    for j in range(len(temp)):
        if 'DECOY' in temp.loc[j]['protein']: decoy_scores.append(float(temp.loc[j][col]))
        else: target_scores.append(float(temp.loc[j][col]))

    f, d = cbf.fdr_calculation(temp)
    line = temp.loc[f][col]

    plt.clf()
    plt.figure()
    plt.hist([target_scores, decoy_scores], 50, stacked=True, density=True)
    plt.axvline(x=line, color='black', linestyle = 'dashed')
    plt.title(col + ', ' + str(f))
    plt.savefig(fileHeader+col+'.png')

'''
delete = ['Data/Output/' + x for x in list(os.listdir('Data/Output/')) if ('FDR' not in x and '31' not in x)]
for x in delete: os.remove(x)
files = ['Data/Output/' + x for x in list(os.listdir('Data/Output/')) if ('peptideFDR' in x and 'corrected' in x)]
def defdict():
    return [0,0]
d = defaultdict(defdict)
good = set()
bad = set()
g = True
for x in files:
    tempDf = pd.read_csv(x)
    match = int(re.sub('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-6peaks_exp-100reps-rep(\d{2,3})_corrected_peptideFDR\.csv', r'\1', x))
    if len(tempDf) > 2000: g = True; good.add(match)
    else: g = False; bad.add(match)

    peptides = list(tempDf['peptide'])
    for p in peptides:
        if g: d[p][0]+= 1
        else: d[p][1]+= 1


for k in d: d[k] = (d[k][0],d[k][1])
print('Number of Peptides: '+str(len(d)))
print('Good: '+str(len(good)))
print('Bad: '+str(len(bad)))
inin = 0
inout = 0
outin = 0
outout = 0

bb = []
bg = []

all_top = []
good_top = []

for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True):
    print([k,v])
#    if v[0] > (len(good)/4)*3 and v[1] > (len(bad)/4)*3: inin += 1
#    elif v[0] > (len(good)/4)*3 and v[1] < len(bad)/4: inout += 1
#    elif v[0] < len(good)/4 and v[1] > (len(bad)/4)*3: outin += 1
    if v[0] > len(good)/2 and v[1] > len(bad)/2: inin += 1; all_top.append(k)
    elif v[0] > len(good)/2 and v[1] < len(bad)/2: inout += 1; good_top.append(k)
    elif v[0] < len(good)/2 and v[1] > len(bad)/2: outin += 1
    else:
        outout += 1
        bg.append(v[0])
        bb.append(v[1])


print(inin)
print(inout)
print(outin)
print(outout)

print(sum(bg)/len(bg))
print(statistics.pstdev(bg))
print(sum(bb)/len(bb))
print(statistics.pstdev(bb))


d2 = defaultdict(int)

def fdr_calculation(df, fdrList=False):
    # initializing the two return values at 0
    fdrValues = []
    numDecoys = 0
    df.fillna("nan",inplace=True)
    # for every row in the dataframe
    for i in range(len(df)):
        # current criteria for 'decoys' is to have 'decoy' in the protein name. This may change in the future.
        if 'DECOY' in df.loc[i]['Name']:
            numDecoys += 1
        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(i+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if len(fdrValues) < 1/0.01:
                if fdrList: return [], 0
                else: return 0, 0
            if fdrList: return fdrValues, numDecoys-1
            else: return len(fdrValues), numDecoys-1
        fdrValues.append(curFDR)
    if fdrList: return fdrValues, numDecoys-1
    else: return len(fdrValues), numDecoys-1

head = 'Data/Input/100reps_searchoutput/'
files = [head+x for x in list(os.listdir(head))]
for x in files:
    tempDf = pd.read_csv(x, sep='\t').sort_values('cosine', ascending=False).reset_index(drop=True)
    tempDf = tempDf.drop_duplicates(subset='Peptide', keep='first').reset_index(drop=True)
    run = re.sub('Data/Input/100reps_searchoutput/2da10ppm20200719_MAGIC_MCF7_1128repro_(\d{2,3})\.txt', r'\1', x)
    hits, decoys = fdr_calculation(tempDf)
    tempDf = tempDf[:hits]
    peptides = list(tempDf['Peptide'])
    for p in peptides:
        d2[p] += 1

msplitPeptides = []
for k in d2:
    if d2[k] > 75: msplitPeptides.append(k)

print(len(msplitPeptides))

data = [
    sorted(list(good)),
    sorted(list(bad)),
    all_top,
    good_top,
    msplitPeptides
]

file = open('Data/Output/peptides.csv', 'w+', newline ='')

# writing the data into the file
with file:
    write = csv.writer(file)
    write.writerows(data)




data = []
with open('Data/Output/peptides.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        data.append(row)

all_top = set(data[2])
good_top = set(data[3])
msplitPeptides = data[4]
'''
'''
peptideFDR = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-100reps-rep16_extraScores.csv')

h = 'Data/Figures/testScoreHistogram_'
histogram(peptideFDR, h, 'cosine')
histogram(peptideFDR, h, 'matchesXcosine')
histogram(peptideFDR, h, 'matches2rootXcosine')
histogram(peptideFDR, h, 'matches3rootXcosine')
histogram(peptideFDR, h, 'matches4rootXcosine')
'''

'''
lib_df = pd.read_csv('Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv', sep='\t')
all_top_df = lib_df[lib_df['PeptideSequence'].isin(all_top)]
good_top_df = lib_df[lib_df['PeptideSequence'].isin(good_top)]
all_top_df.to_csv('Data/Input/human_31peaks_noloss_400to2000_pt2mz_allTop.csv')
good_top_df.to_csv('Data/Input/human_31peaks_noloss_400to2000_pt2mz_goodTop.csv')

for i in range(len(all_top)):
    all_top[i] = re.sub('\(UniMod:\d+\)','',all_top[i])
all_top = list(set(all_top))

for i in range(len(good_top)):
    good_top[i] = re.sub('\(UniMod:\d+\)','',good_top[i])
good_top = list(set(good_top))

for i in range(len(msplitPeptides)):
    msplitPeptides[i] = re.sub('\+\d+\.\d+','',msplitPeptides[i])
msplitPeptides = list(set(msplitPeptides))


venn3([set(all_top),set(good_top),set(msplitPeptides)], set_labels = ["csoDIAq, both clusters","csoDIAq, good cluster only", "MSPLIT-DIA"])
plt.title('Comparing Peptide Identification Outputs (stripped sequences, 6peaks)\n')
plt.savefig('Data/Figures/peptide_comparison_clusters_msplit75.png')
'''
'''
files = ['Data/Output/'+x for x in list(os.listdir('Data/Output/')) if 'extraScores' not in x]
#for x in files: os.remove(x)
'''
'''
files = list(os.listdir('Data/Output/'))
delete = set()
for x in files:
    temp = re.search('csodiaq_lib-human-noloss-400to2000-pt2mz-\d{1,2}peaks_exp-100reps-rep\d{2,3}(_corrected)?\.csv',x)#r'csodiaq_lib-human-noloss-400to2000-pt2mz-\1peaks_exp-100reps-rep\2\3.csv',x)
    if not temp: delete.add('Data/Output/'+x)
#for x in sorted(delete): print(x)
for x in delete: os.remove(x)
'''

files = ['Data/Output/'+x for x in list(os.listdir('Data/Output/'))]
for file in sorted(files):
    print(file)
    df = pd.read_csv(file)
    sl = 11
    scores = [[] for i in range(sl)]
    for i in range(len(df)):
        match = int(df.loc[i]['shared'])
        cosine = float(df.loc[i]['cosine'])
        for j in range(sl):
            scores[j].append((match**(1/(j+1)))*cosine)

    for i in range(sl):
        title = 'match'+str(i+1)+'rootXcosine'
        df[title] = scores[i]


    newFile = re.sub('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-(\d{1,2})peaks_exp-100reps-rep(\d{2,3})(_corrected)?\.csv',r'Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-\1peaks_exp-100reps-rep\2\3_extraScores.csv',file)
    print(newFile+'\n')
    df.to_csv(newFile)

'''
df = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-100reps-rep16_corrected.csv')
scores = [
    [],
    [],
    [],
    []
]
for i in range(len(df)):
    match = int(df.loc[i]['shared'])
    cosine = float(df.loc[i]['cosine'])
    scores[0].append(match*cosine)
    scores[1].append((match**(1/2))*cosine)
    scores[2].append((match**(1/3))*cosine)
    scores[3].append((match**(1/4))*cosine)

df['matchesXcosine'] = scores[0]
df['matches2rootXcosine'] = scores[1]
df['matches3rootXcosine'] = scores[2]
df['matches4rootXcosine'] = scores[3]
df.to_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-100reps-rep16_corrected_extraScores.csv')
'''
