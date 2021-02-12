import pandas as pd
#import os
import csodiaq_base_functions as cbf
#import re
from collections import defaultdict
#import statistics
#import csv
#from matplotlib_venn import venn3
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure
import numpy as np
#from matplotlib.ticker import PercentFormatter
#import math
#from scipy.stats import linregress
import pickle
import umap
import seaborn as sns
from sklearn.preprocessing import StandardScaler

def set_plot_settings(xlabel, ylabel, wide=True):
    if wide: pyplot.figure(figsize=(18,12))
    else: pyplot.figure(figsize=(12,12))
    pyplot.xlim(1.5,10)
    pyplot.ylim(1,10)
    pyplot.axhline(linewidth=4, y=1, color='black')
    pyplot.axhline(linewidth=4, y=10, color='black')
    pyplot.axvline(linewidth=4, x=1.5, color='black')
    pyplot.axvline(linewidth=4, x=10, color='black')
    pyplot.xlabel(xlabel, fontsize = 36, weight='bold')
    pyplot.ylabel(ylabel, fontsize = 36, weight='bold')
    pyplot.tick_params(axis="x", labelsize=45)
    pyplot.tick_params(axis="y", labelsize=45)

def histogram(df, fileHeader, col, title, l=False):
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
    set_plot_settings(title, col, 'Frequency')

    plt.hist([target_scores, decoy_scores], 50, stacked=True, density=True)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    if l: plt.axvline(x=line, color='black', linestyle = 'dashed')

    plt.savefig(fileHeader+col+'.png')

def corr(X, Y):
    """Computes the Pearson correlation coefficient and a 95% confidence
    interval based on the data in X and Y."""

    r = np.corrcoef(X, Y)[0,1]
    f = 0.5*np.log((1+r)/(1-r))
    se = 1/np.sqrt(len(X)-3)
    ucl = f + 2*se
    lcl = f - 2*se

    lcl = (np.exp(2*lcl) - 1) / (np.exp(2*lcl) + 1)
    ucl = (np.exp(2*ucl) - 1) / (np.exp(2*ucl) + 1)

    return r,lcl,ucl

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
peptideFDR = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-100reps-rep16_corrected_2SD.csv')

h = 'Data/Figures/testScoreHistogram_newScore_'
histogram(peptideFDR, h, 'cosine', 'Distribution of Cosine Scores')
histogram(peptideFDR, h, 'shared', 'Distribution of Matched Fragments')
histogram(peptideFDR, h, 'csoDIAq_Score', 'Distribution of csoDIAq Scores', l=True)
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
protFiles = set(['Data/Output/'+x for x in list(os.listdir('Data/Output/')) if 'proteinFDR' in x])
files = ['Data/Output/'+x for x in list(os.listdir('Data/Output/'))]
delete = set()
for x in files:
    temp = re.search('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-\d{1,2}peaks_exp-100reps-rep\d{2,3}(?:_corrected_2SD)?\.csv',x)#r'csodiaq_lib-human-noloss-400to2000-pt2mz-\1peaks_exp-100reps-rep\2\3.csv',x)

    if not temp and x not in protFiles: delete.add(x)
#keep = sorted([x for x in files if x not in delete])
#for x in keep: print(x)
for x in delete: os.remove(x)
'''
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
'''
df = pd.read_csv('Data/Input/human_31peaks_noloss_400to2000_pt2mz_subset.csv')
d = dict(tuple(df.groupby(['FullUniModPeptideName','PrecursorCharge'])))
for key in d:
    print(key)
    print(d[key])
    print('-'*20)
'''
'''
cols = [
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to8_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to4_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to2_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to1_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_2to1_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_4to1_01.mzXML',
    'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_8to1_01.mzXML'
]
finalResults = []

libs = [3, 10, 20]
rgx = re.compile(r'Data/oldOutput/lib\d+/test_(\d+)_(mean|median)_(\d+)_(P?\d+)\.csv')
rgx2 = re.compile(r'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_(\d)to(\d)_01\.mzXML')
for lib in libs:
    header = 'Data/oldOutput/lib'+str(lib)+'/'
    files = [header + x for x in list(os.listdir(header))]
    for file in files:
        df = pd.read_csv(file)
        reg = rgx.match(file)
        match = reg.group(1)
        m = reg.group(2)
        expMin = reg.group(3)
        stdev = reg.group(4)
        if stdev=='P5': stdev='0.5'
        main = [lib, match, m, expMin, stdev]
        outFile = 'Data/QuantifyCompare/lib'+str(lib)+'_match'+match+'_'+m+'_expressMin'+expMin+'_stdev'+stdev+'.csv'
        print(outFile)
        for col in cols:
            proportions = rgx2.match(col)
            prop = float(proportions.group(1))/float(proportions.group(2))
            values = []
            v1 = len(df[df[col]==-100])
            v2 = len(df[(df[col]<=-0.1) & (df[col]>-100)])
            v3 = len(df[(df[col]>=-0.1) & (df[col]<0.1)])
            v4 = len(df[(df[col]>=0.1) & (df[col]<100)])
            v5 = len(df[df[col]==100])
            v6 = len(df[df[col].isna()])
            values = [v1, v2, v3, v4, v5, v6]
            finalResults.append([prop] + main + values)
fDf = pd.DataFrame(finalResults, columns=['light:heavy', 'lib', 'matches', 'calculation', 'minimumExpression', 'stdev', 'onlyHeavy','mostlyHeavy','aboutSame','mostlyLight','onlyLight','neither'])
fDf.to_csv('Data/QuantifyCompare/full.csv', index=False)
'''

dfLC = pd.read_csv('Data/MQpeptides_Quant.csv')
dfDIA = pd.read_csv('Data/lib3-1-median-1.csv')



#col1 = 'Ratio H/L 1to8'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to8_01.mzXML'

#col1 = 'Ratio H/L 1to4'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to4_01.mzXML'

#col1 = 'Ratio H/L 1to2'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to2_01.mzXML'

#col1 = 'Ratio H/L 1to1'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to1_01.mzXML'

#col1 = 'Ratio H/L 2to1'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_2to1_01.mzXML'

#col1 = 'Ratio H/L 4to1'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_4to1_01.mzXML'

#col1 = 'Ratio H/L 8to1'
#col2 = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_8to1_01.mzXML'

'''
def calc_med_stdev(dfLC, dfDIA, ratio, ori=1, expr=False):
    if ori:
        colLC = 'Ratio H/L 1to'+str(RATIO)
        colDIA = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to'+str(RATIO)+'_01.mzXML'
        ratio = RATIO
    else:
        colLC = 'Ratio H/L '+str(RATIO)+'to1'
        colDIA = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_'+str(RATIO)+'to1_01.mzXML'
        ratio = -RATIO
    if RATIO==1: ratio=0

    print(ratio)
    ret = []
    dfLC[colLC] = np.log2(dfLC[colLC])
    appLC = dfLC[~dfLC[colLC].isnull()][colLC]
    appDIA = dfDIA[(~dfDIA[colDIA].isnull())][colDIA]
    ret.append([np.median(appLC),statistics.pstdev(appLC),ratio, 'LC'])

    e = np.median(dfDIA[~dfDIA[colDIA].isnull()][colDIA])
    if expr: ret.append([e,statistics.pstdev(appDIA),ratio, 'DIA'])
    else: ret.append([np.median(appDIA),statistics.pstdev(appDIA),ratio, 'DIA'])

    return ret
'''
'''
RATIOS = [1, 2, 4, 8]
data = []
e = False
for RATIO in RATIOS:

    data += calc_med_stdev(dfLC, dfDIA, RATIO, expr=e)
    if RATIO != 1: data += calc_med_stdev(dfLC, dfDIA, RATIO, ori=0, expr=e)
'''
'''
# For Making Line Graphs
def calc_key_variables(df, r, type, ori=1):
    if ori:
        if type=='LC': col = 'Ratio H/L 1to'+str(r)
        else: col = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to'+str(r)+'_01.mzXML'
        ratio = r
    else:
        if type=='LC': col = 'Ratio H/L '+str(r)+'to1'
        else: col = 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_'+str(r)+'to1_01.mzXML'
        ratio = -r
    if r==1: ratio=0


    if type=='LC': df[col] = np.log2(df[col])
    app = df[~df[col].isnull()][col]
    med = np.median(app)

    ratioDict = {
        -8:-3,
        -4:-2,
        -2:-1,
        0:0,
        2:1,
        4:2,
        8:3
    }
    #if type == 'DIA': med = med*3
    return [[med,statistics.pstdev(app),ratioDict[ratio], type]]

def calc_all_variables(df, type, ratios):
    data = []
    for ratio in ratios:
        data += calc_key_variables(df, ratio, type)
        if ratio != 1: data += calc_key_variables(df, ratio, type, ori=0)
    finalDf = pd.DataFrame(data, columns=['median','stdev','ratio','type'])
    finalDf = finalDf.sort_values('ratio', ascending=False).reset_index(drop=True)
    return finalDf

RATIOS = [1, 2, 4, 8]

ratioDict = {
        -8:-3,
        -4:-2,
        -2:-1,
        0:0,
        2:1,
        4:2,
        8:3
    }
lc = calc_all_variables(dfLC, 'LC', RATIOS)
dia = dfDIA.drop(labels=['numPeps'], axis=1)
diaRats = list(dia['ratio'])
for i in range(len(diaRats)): diaRats[i] = ratioDict[diaRats[i]]
dia['ratio'] = diaRats

diaMeans = list(dia['median'])
diaMeans = [1.265**x for x in diaMeans]
dia['median'] = diaMeans


print(lc)
print(dia)


X = lc['median']
Y = dia['median']
lr = ['test']+list(linregress(X, Y))
print(pd.DataFrame([lr], columns=['name', 'slope', 'intercept', 'rvalue', 'pvalue', 'stderr' ]))
print(lr)
pyplot.scatter(X,Y)
pyplot.ylim(-3.3, 3.3)
pyplot.show()
'''
'''
#For Making Line Graphs
finalDf = pd.concat([lc, dia])

finalDf.to_csv('Data/median.csv')

fig, ax = pyplot.subplots()

for key, group in finalDf.groupby('type'):
    group.plot('ratio', 'median', yerr='stdev',
               label=key, ax=ax)

pyplot.show()
'''
'''
dictLC = dict(zip(dfLC['Sequence'], dfLC[col1]))
dictDIA = dict(zip(dfDIA['peptide'], dfDIA[col2]))

intersect = set.intersection(set(dictLC.keys()),set(dictDIA.keys()))

X, Y = [], []
for key in intersect:
    LCcheck = ~np.isnan(dictLC[key])
    DIAcheck = (dictDIA[key]!=5 and dictDIA[key]!=-5 and ~np.isnan(dictDIA[key]))
    if LCcheck and DIAcheck:
        X.append(dictLC[key])
        Y.append(dictDIA[key])



pyplot.scatter(X,Y)
pyplot.xlabel('LC')
pyplot.ylabel('DIA')
pyplot.show()
'''
'''
# Full Compare Compilation
head = 'Data/QuantifyCompare/variables/'#lib3_compare.csv'
#files = [head+x for x in list(os.listdir(head))]

files = [
    'Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_mgf.csv',
    'Data/QuantifyCompare/variables/compare.csv'
]

dfList = []
for file in files:
    dfList.append(pd.read_csv(file))
finalDf = pd.concat(dfList)
finalDf.to_csv(head+'full_compare.csv', index=False)
'''
'''
libDict = pickle.load(open( "Data/Input/TempHold/lib_3_.p", "rb" ))
fragDict = pickle.load(open( "Data/Input/TempHold/frag.p", "rb" ))

jesse = pd.read_csv("Data/jesseQuant.csv")
jesseLights = {}
jesseHeavys = {}
jesseCount = defaultdict(int)
for i in range(len(jesse)):
    index = jesse.loc[i]['quantscan']
    jesseCount[index] += 1
    lights = jesse.loc[i]['ylight'].strip("[]").replace("'","").split(', ')
    heavys = jesse.loc[i]['yheavy'].strip("[]").replace("'","").split(', ')
    if len(lights[0]) != 0: jesseLights[int(index)] = set([round(float(x),1) for x in lights])
    else: jesseLights[int(index)] = set()
    if len(heavys[0]) != 0: jesseHeavys[int(index)] = set([round(float(x),1) for x in heavys])
    else: jesseHeavys[int(index)] = set()



calebLights = {}
calebHeavys = {}

for key in libDict:
    lights = set()
    heavys = set()
    for peak in libDict[key]:
        if peak[2][0]=='light': lights.add(round(peak[0],1))
        else: heavys.add(round(peak[0],1))
    calebLights[int(key)] = lights
    calebHeavys[int(key)] = heavys

Jkeys = set(list(jesseLights.keys()))
Ckeys = set(list(calebLights.keys()))
'''
'''
# comparing peptides between outputs
jesse = pd.read_csv('Data/Input/TempHold/mostintense_quantmzlist.txt', sep='\t')
caleb = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-n1b_corrected_mostIntenseTargs_allCVs.csv')

jesse_peptides = set(list(tuple(zip(jesse['Peptide'], jesse['CV']))))
caleb_peptides = set(list(tuple(zip(caleb['peptide'], caleb['CompensationVoltage']))))

both = set.intersection(jesse_peptides, caleb_peptides)
print(len(jesse_peptides))
print(len(caleb_peptides))
print(len(both))
'''

# Finding best variables
allVariables = pd.read_csv('Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_mgf.csv')
#allVariables = pd.read_csv('Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_traml.csv')
print('no filter:')
print(len(allVariables))
print('\n')

# filter by average standard deviation
allVariables = allVariables[allVariables['avrStDev'] < 1]
print('filter out average standard deviations above 1:')
print(len(allVariables))
print('\n')

# filter by number of peptides
allVariables = allVariables[allVariables['numPeps'] > (3864*.8)]
#allVariables = allVariables[allVariables['numPeps'] > (2926*.75)]
print('filter out peptides < 75%:')
print(len(allVariables))
print('\n')

# filter slopes not close to 1
allVariables = allVariables[(allVariables['slope'] > .95) & (allVariables['slope'] < 1.05)]
print('filter out slopes outside .5 from 1:')
print(len(allVariables))
print('\n')

print(allVariables)
topHits = allVariables.index


# UMAP attempt
allVariables = pd.read_csv('Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_mgf.csv')
#t1 = ['one stdev' if x=='1' else x for x in allVariables['corrStDev']]
#allVariables['corrStDev'] = t1
#print(t1)
#allVariables.to_csv('Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_mgf.csv', index=False)
allVariables = allVariables[(allVariables['mode']=='mean') | (allVariables['mode']=='median')]#.reset_index(drop=True)
#topHits = [101, 133, 149, 181, 213, 229]



reducer = umap.UMAP()
print(1)
data = allVariables[
    ['slope',
#    'intercept',
#    'rvalue',
#    'pvalue',
#    'stderr',
    'avrStDev',
    'numPeps']
].values
scaled_data = StandardScaler().fit_transform(data)
embedding = reducer.fit_transform(scaled_data)


'''
#interactive UMAP
from bokeh.plotting import figure, show, output_notebook
from bokeh.models import HoverTool, ColumnDataSource, CategoricalColorMapper
from bokeh.palettes import Spectral10

digits_df = pd.DataFrame(embedding, columns=('x', 'y'))
print([allVariables.loc[i] for i in range(len(allVariables))])
digits_df['digit'] = ['_'.join([str(x) for x in list(allVariables.loc[i])]) for i in range(len(allVariables))]


datasource = ColumnDataSource(digits_df)


plot_figure = figure(
    title='UMAP projection of the Digits dataset',
    plot_width=600,
    plot_height=600,
    tools=('pan, wheel_zoom, reset')
)

plot_figure.add_tools(HoverTool(tooltips="""
<div>
    <div>
        <img src='@image' style='float: left; margin: 5px 5px 5px 5px'/>
    </div>
    <div>
        <span style='font-size: 16px; color: #224499'>Digit:</span>
        <span style='font-size: 18px'>@digit</span>
    </div>
</div>
"""))

plot_figure.circle(
    'x',
    'y',
    source=datasource,
    color=dict(field='digit'),
    line_alpha=0.6,
    fill_alpha=0.6,
    size=4
)
show(plot_figure)
'''

#continue UMAP
x_values = embedding[:, 0]
y_values = embedding[:, 1]
#pickle.dump(x_values, open( "Data/x.p", "wb" ))
#pickle.dump(y_values, open( "Data/y.p", "wb" ))
x_values = pickle.load(open( "Data/x.p", "rb" ))
y_values = pickle.load(open( "Data/y.p", "rb" ))

ins = ['numPeps','avrStDev','slope']
outs = ['corrStDev','libPeaks','minMatch']
edgeColors = ['#DC143C' if x in topHits else 'white' for x in allVariables.index]
edgeWidth = [6 if x in topHits else 0.5 for x in allVariables.index]
palDict = {
    'avrStDev':'viridis',
    'numPeps':'viridis_r',
    'slope':'viridis_r'
    #'numPeps':'Purples',
    #'slope':'YlOrBr'
}


for i in ins:
    for o in outs:

        labels = allVariables[i]
        markers = allVariables[o]
        print(len(allVariables))
        #sns.set_palette("viridis")
        h = None
        if i=='slope': h = (None, 1)


        set_plot_settings('','')
        #style=markers,
        ax = sns.scatterplot(x=x_values, y=y_values, hue=labels, style=markers, palette=palDict[i], edgecolor=edgeColors, hue_norm=h, s=750, linewidth=edgeWidth)
        norm = plt.Normalize(allVariables[i].min(), allVariables[i].max())
        sm = plt.cm.ScalarMappable(cmap=palDict[i], norm=norm)
        sm.set_array([])

        ax.get_legend().remove()
        cbar = ax.figure.colorbar(sm,orientation="horizontal")
        cbar.ax.tick_params(labelsize=45)
        #plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='City Area')
        plt.savefig('Data/QuantifyCompare/UMAP/'+i+'_'+o+'.png')
        plt.close("all")

'''
#violin plots
# Chosen: all, custom, median, custom
ins = ['numPeps','avrStDev','slope']

colDict = {
    'avrStDev':'blueviolet',
    'numPeps':'skyblue',
    'slope':'red'
}
yDict = {
    'avrStDev':0.9849648653900464,
    'numPeps':3489,
    'slope':0.9636899678653778
}


allVariables = pd.read_csv('Data/QuantifyCompare/variables/full_compare_noSqrRt_HalvingBot10_mgf.csv')
for i in ins:
    df = allVariables.loc[:, allVariables.columns.intersection([i])]#'avrStDev','slope','numPeps'
    df = df.melt(var_name='groups', value_name='vals')
    sns.violinplot(x="groups", y="vals", data=df, color=colDict[i])
    plt.axhline(y=yDict[i], color='red', linestyle = 'dashed') #slope
    plt.savefig('Data/QuantifyCompare/ViolinPlots/'+i+'.png')
    plt.close("all")
'''
'''
# Trying to match peptides between LC and DIA (fail)
def calc_heavy_mz(seq, mz, z):
    hK=8.014199 ## mass of heavy lysine
    hR=10.00827 ## mass of heavy arg

    nK = seq.count('K')
    nR = seq.count('R')
    #nK = len(seq) - len(re.sub('K','',seq))
    #nR = len(seq) - len(re.sub('R','',seq))

    heavyMz = mz + (nK*hK)/z + (nR*hR)/z
    return heavyMz

# Matching Peptides
from pyteomics import mass

dfLC = pd.read_csv('Data/MQpeptides_Quant.csv')
dfLC = dfLC.loc[:, dfLC.columns.intersection(['Sequence','Charges'])]
dfDIA = pd.read_csv('Data/Input/TempHold/mostintense_quantmzlist.txt', sep='\t')


LCDict = defaultdict(list)

print('LC start')
for index, row in dfLC.iterrows():
    seq = row['Sequence']
    charges = [int(c) for c in row['Charges'].split(';')]
    for c in charges:
        light = mass.calculate_mass(sequence=seq, charge=c)
        heavy = calc_heavy_mz(seq, light, c)
        key = (round(light,2),round(heavy,2),c)
        LCDict[key].append(seq)

print('DIA start')
DIADict = defaultdict(list)
for index, row in dfDIA.iterrows():
    key = (round(row['prec_light_mz'],2),round(row['prec_heavy_mz'],2),int(row['z']))
    DIADict[key].append(row['Peptide'])

print('key start')
LCKeys = set(LCDict.keys())
DIAKeys = set(DIADict.keys())
both = set.intersection(LCKeys, DIAKeys)
for key in both:
    print(key)
    print(LCDict[key])
    print(DIADict[key])
    print('\n')

print(len(both))
'''
