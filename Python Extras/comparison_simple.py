from pyteomics import mgf
import random
import pickle
import os
import pandas as pd
import re

random.seed(0)
'''
#Reading MGF file, selecting 10 random (non-Decoy)
mg = {}
count = 0
with mgf.read('Data/Input/human.faims.fixed.decoy.mgf') as reader:
    for spectrum in reader:
        if count % 10000 == 0:
            print(count)
        count += 1
        if ('DECOY' not in spectrum['params']['title']) and (not bool(re.match("\+\d+\.\d+",spectrum['params']['seq']))):
            mg[spectrum['params']['title']] = spectrum
cond = random.sample(list(mg.keys()), 500)
mg2 = {}
for x in cond: mg2[x] = mg[x]

print(mg2.keys())
with open('Data/Input/mgf_condensed_500.pkl','wb') as f:
    pickle.dump(mg2,f)
for i in mg2:
    print(mg2[i]['params']['seq'])
'''

#Loading in condensed MGF
with open('Data/Input/mgf_condensed_500.pkl','rb') as f:
    mg = pickle.load(f)

'''
#Reading big TRAML file, writing small TRAML file
seqs = []
charges = []
for key in mg:
    seqs.append(mg[key]['params']['seq'])
    charges.append(int(mg[key]['params']['charge'][0]))
ts = pd.read_csv("Data/Input/iproph-speclib_con_decoys31.tsv", sep='\t')
ts2 = ts[ts['FullUniModPeptideName'].isin(seqs) ]
#ts2.to_csv('condensed_seq.csv')
#ts2 = pd.read_csv('condensed_seq.csv')
temp = [ts2[ts2['FullUniModPeptideName'] == x] for x in seqs]
temp = [ temp[i][ temp[i]['PrecursorCharge'] == charges[i] ] for i in range(len(temp))]

ts3 = pd.concat(temp)

ts3.to_csv('Data/Input/condensed_500.csv')

'''

#Reading small TRAML file
ts = pd.read_csv("Data/Input/condensed_500.csv")
count = 0
for x in mg:
    seq = mg[x]['params']['seq']
    print(seq)
    mgf_mzs = mg[x]['m/z array']
    mgf_intensities = mg[x]['intensity array']
    # top 20 intensity
    #indices = sorted(range(len(mg[x]['intensity array'])), key=lambda i: mg[x]['intensity array'][i])[-20:]
    #mgf_mzs = [ mg[x]['m/z array'][i] for i in indices]
    #mgf_intensities = [ mg[x]['intensity array'][i] for i in indices]
    mgf_mzs, mgf_intensities = zip(*sorted(zip(mgf_mzs, mgf_intensities)))
    ts2 = ts[ts['FullUniModPeptideName']==seq]
    if len(ts2) == 0:
        continue
    tsv_mzs = ts2['ProductMz'].values.tolist()
    tsv_mzs = [ round(y,3) for y in tsv_mzs ]
    tsv_intensities = ts2['LibraryIntensity'].values.tolist()
    tsv_intensities = [ round(y,1) for y in tsv_intensities ]

    max_intensities = max(mgf_intensities)
    mgf_intensities = [ (x/max_intensities)*10000 for x in mgf_intensities ]
    mgf_intensities = [ round(y,1) for y in mgf_intensities ]

    tsv_mzs, tsv_intensities = zip(*sorted(zip(tsv_mzs, tsv_intensities)))

    #Make this a visual - is the shape the same on a graph?
    # as a bonus, you could probably calculate the cosine similarity to get a numeric representation
    print(seq)
    print('mgf length: ' + str(len(mgf_mzs)))
    print('tsv length: ' + str(len(tsv_mzs)))
    print('mzs:')
    print(mgf_mzs)
    print(tsv_mzs)
    print('intensities:')
    print(mgf_intensities)
    print(tsv_intensities)
    print()

# Graph
def plot_spec(SPECTRA, COLOR):
    plt.vlines(SPECTRA.columns, np.repeat(0, len(SPECTRA.columns)), SPECTRA, colors=COLOR)
