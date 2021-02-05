import pandas as pd
import os
import csodiaq_base_functions as cbf
import re
from collections import defaultdict
import statistics
import csv
from matplotlib_venn import venn3
import matplotlib.pyplot as pyplot
from matplotlib.pyplot import figure
import numpy as np
from matplotlib.ticker import PercentFormatter
import math
from scipy.stats import linregress
import pickle


def set_plot_settings(title, xlabel, ylabel):
    plt.figure(figsize=(10,10))
    plt.title(title, fontsize = 32)
    plt.xlabel(xlabel, fontsize = 24)
    plt.ylabel(ylabel, fontsize = 24)
    plt.tick_params(axis="x", labelsize=18)
    plt.tick_params(axis="y", labelsize=18)

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
files = [head+x for x in list(os.listdir(head))]
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
jesse = pd.read_csv('Data/Input/TempHold/mostintense_quantmzlist.txt', sep='\t')
caleb = pd.read_csv('Data/Output/csodiaq_lib-human-noloss-400to2000-pt2mz-31peaks_exp-n1b_corrected_mostIntenseTargs_allCVs.csv')

jesse_peptides = set(list(tuple(zip(jesse['Peptide'], jesse['CV']))))
caleb_peptides = set(list(tuple(zip(caleb['peptide'], caleb['CompensationVoltage']))))

both = set.intersection(jesse_peptides, caleb_peptides)
print(len(jesse_peptides))
print(len(caleb_peptides))
print(len(both))
