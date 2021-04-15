import pandas as pd
from pyteomics import mzxml, mgf, mass
from collections import defaultdict
import csodiaq_base_functions as cbf
from random import randint
import re
from Bio import SeqIO
import bisect
from matplotlib import pyplot
from matplotlib_venn import venn2
import numpy as np
import pickle
import bz2
from numba import njit
from timeit import default_timer as timer
import io

pd.set_option('display.max_columns', None)



def fdr_calculation(df, returnType=0): #***NOTE***make return type
    # initializing the two return values at 0
    fdrValues = []
    indices = []
    numDecoys = 0
    df.fillna("nan",inplace=True)
    # for every row in the dataframe
    count = 0


    for index, row in df.iterrows():
        # current criteria for 'decoys' is to have 'decoy' in the protein name. This may change in the future.
        if (returnType==0 or returnType==1): decoy = 'DECOY' in row['Name']
        elif returnType==2: decoy = row['decoy']
        if decoy:
            numDecoys += 1

        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(count+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if len(fdrValues) < 1/0.01:
                if returnType: return [], 0
                else: return 0, 0
            if returnType==1: return fdrValues, numDecoys-1
            if returnType==2: return indices, numDecoys-1
            else: return len(fdrValues), numDecoys-1
        fdrValues.append(curFDR)
        indices.append(index)
        count += 1

    if returnType: return fdrValues, numDecoys-1
    else: return len(fdrValues), numDecoys-1

def return_frag_mzs(peptide, z):
    mzValues = []
    digPat = r'\+\d+\.\d+'
    digs = re.findall(digPat, peptide)
    pepFrags = re.split(digPat, peptide)
    modValues = {}
    seq = ''
    while len(digs) != 0:
        dig = digs.pop(0)
        frag = pepFrags.pop(0)
        seq += frag
        modValues[len(seq)] = float(dig[1:])/z
    seq += pepFrags[0]
    for i in range(1, len(seq)-1):
        mz = mass.fast_mass(sequence=seq[i:], ion_type='y', charge=z)
        mz += sum([modValues[x] for x in modValues if x > i])
        mzValues.append(mz)

    for i in range(len(seq)-1, 1, -1):
        mz = mass.fast_mass(sequence=seq[:i], ion_type='b', charge=z)
        mz += sum([modValues[x] for x in modValues if x <= i])
        mzValues.append(mz)
    return mzValues

def approx(x, y, ppmTol):
    if x==y: return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)

def approx_list(x, l, ppmTol=10):
    for i in range(len(l)):
        if approx(x, l[i], ppmTol): return i
    return -1

#if 'protein' in spec['params']

def clean_mgf_file(file):
    spectra = mgf.read(file)
    fasta = 'C:/Users/ccranney/Desktop/Caleb_Files/data/2019-03-14-td-UP000005640.fasta'
    fDict = {}
    longPep = ''
    for record in SeqIO.parse(open(fasta,'r'),'fasta'):
        fDict[len(longPep)] = record.id
        longPep += str(record.seq) + '.'
    cleaned = []
    count = 0
    pepCount = 0
    for spec in spectra:
        count += 1
        #if count % 40==0: break
        #mzValues = return_frag_mzs(spec['params']['seq'],1)
        #peaks = list(tuple(zip(spec['m/z array'],spec['intensity array'])))
        #for i in range(len(peaks)-1,-1,-1):
        #    if approx_list(peaks[i][0],mzValues)==-1: peaks.pop(i)
        #if len(peaks)==0: continue
        #peaks.sort(key=lambda x:x[0])
        #spec['m/z array'],spec['intensity array'] = map(list,zip(*peaks))
        #'''
        decoy = False
        if 'protein' in spec['params'] and 'DECOY' in spec['params']['protein']: decoy = True
        else:
            seq = re.sub(r'\+\d+\.\d+', '', spec['params']['seq'])
            listOfI = [m.start() for m in re.finditer(seq, longPep)]
            sorted_keys = sorted(fDict.keys())
            proteins = set()
            for i in listOfI:
                insertion_point = bisect.bisect_left(sorted_keys,i)
            # adjust, as bisect returns not exactly what we want
                if insertion_point==len(sorted_keys) or sorted_keys[insertion_point]!=i:
                    insertion_point-=1
                protein = fDict[sorted_keys[insertion_point]]
                proteins.add(fDict[sorted_keys[insertion_point]])
            if len(proteins)==0: proteins.add(spec['params']['seq'])

        if decoy: proteins = ['DECOY_0_'+x for x in proteins]

        protein = str(len(proteins)) + '/' + '/'.join(sorted(proteins))
        spec['params']['protein'] = protein
        if protein != '0/': pepCount += 1
        #'''
        cleaned.append(spec)
        if count % 1000 == 0: print(count); print(pepCount); print(protein)


    cleanedFile = re.sub('(.*).mgf', r'\1_proteinsAdded.mgf', file)
    mgf.write(cleaned, cleanedFile)

print('\n'*30)


'''
scanCount = defaultdict(int)
msplit = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/HeLa-50ppm/MSPLIT_HeLa.txt', sep='\t')
scans = msplit['Scan#']
for scan in scans:
    scanCount[scan]+=1

#hits, decoys = fdr_calculation(msplit)
#print(hits, decoys)
'''
'''
import operator

sorted_tuples = sorted(scanCount.items(), key=operator.itemgetter(1),reverse=True)
temp = []
for x in sorted_tuples[:20]: temp.append(x[0])


#print(msplit[msplit['Scan#']==22501])

#hits, decoys = fdr_calculation(msplit)

#print(hits, decoys)

csodiaq = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/HeLa-50ppm/CsoDIAq-file1_HeLa_160min_DIA_106win_1.csv')
#print(len(csodiaq))
'''
'''
windows = set()
spectralCount = 0
with mzxml.read('C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1.mzXML', use_index=True) as spectra:
#with mzxml.read('C:/Users/ccranney/Desktop/Caleb_Files/data/18302_REP2_500ng_HumanLysate_SWATH_2 (2)..mzXML', use_index=True) as spectra:

    for scan in temp:
        spec = spectra.get_by_id(str(scan))
        print(scan)
        print('PrecursorMz: '+str(spec['precursorMz'][0]['precursorMz']))
        print('Window: '+str(spec['precursorMz'][0]['windowWideness']))

        tempDf = msplit[msplit['Scan#']==scan]
        print('Mz Values:')
        print(tempDf['Mz.1'])
        print('\n')

    for spec in spectra:
        spectralCount += 1
        if spectralCount % 10000 == 0: print(spectralCount)
        if 'precursorMz' in spec:
            windows.add(spec['precursorMz'][0]['windowWideness'])
            if spec['precursorMz'][0]['windowWideness'] == 50.0: print(spec['num'])
for x in windows: print(x)
#'''
'''
file = 'C:/Users/ccranney/Desktop/Caleb_Files/data/18302_REP2_500ng_HumanLysate_SWATH_2 (2)..mzXML'
with open(file, 'rb') as f:
    for _ in range(100): # first 10 lines
        print(f.readline())
#'''
'''
peptide = 'A+42.01057AAAAAGAGPEM+15.9949VR'
#peptide = 'AAAAA'
cbf.return_frag_mzs(peptide, 1)
#'''
'''
file = 'C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf'
fasta = 'C:/Users/ccranney/Desktop/Caleb_Files/data/2019-03-14-td-UP000005640.fasta'
cbf.clean_mgf_file(file, fasta)
#'''
'''
mgf1 = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf')
mgf2 = []
count = 0
print('enter processing')
for spec in mgf1:
    spec['intensity array'] = spec['intensity array'][:10]
    spec['m/z array'] = spec['m/z array'][:10]
    mgf2.append(spec)
    #spec['params']['title'] = count
    #count += 1
    #if count == 10: break
print(mgf2)
print('enter writing')
mgf.write(mgf1, 'C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy3.mgf')
#'''
'''
fasta = 'C:/Users/ccranney/Desktop/Caleb_Files/data/2019-03-14-td-UP000005640.fasta'
fDict = {}
longPep = ''
for record in SeqIO.parse(open(fasta,'r'),'fasta'):
    fDict[len(longPep)] = record.id
    longPep += str(record.seq) + '.'

print(len(fDict))
for i in sorted(fDict)[:3]:
    print(i)
    if i != 0: print(longPep[i-1])
    print(fDict[i])
print(longPep[:200])

pep = 'KMMKRRGINVSDE'
listOfI = [m.start() for m in re.finditer(pep, longPep)]
print(listOfI)
sorted_keys = sorted(fDict.keys())
proteins = set()
for i in listOfI:
    insertion_point = bisect.bisect_left(sorted_keys,i)
# adjust, as bisect returns not exactly what we want
    if insertion_point==len(sorted_keys) or sorted_keys[insertion_point]!=i:
        insertion_point-=1
    proteins.add(fDict[sorted_keys[insertion_point]])
    print(insertion_point)
    print(fDict[sorted_keys[insertion_point]])

proteins = str(len(proteins)) + '/' + '/'.join(sorted(proteins))
print(proteins)
#'''
'''
fasta = 'C:/Users/ccranney/Desktop/Caleb_Files/data/2019-03-14-td-UP000005640.fasta'

record_dict = SeqIO.index(fasta, "fasta")
spectra = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy_cleaned.mgf')
count = 0
num = 5000
for spec in spectra:
    count += 1
    if count % 1000 == 0: print(count)
    if randint(0,num) % num != 0: continue

    seq = re.sub(r'\+\d+\.\d+', '', spec['params']['seq'])
    proteins = spec['params']['protein'].split('/')[1:]
    for protein in proteins:
        if 'DECOY' in protein: continue
        temp = str(record_dict[protein].seq)
        if seq not in temp:
            print('FAIL ' + 'X'*50)
            print(seq)
            print(temp)
            print('\n')
        else:
            print('SUCCESS')
            print('\n')
#'''
'''
peptides = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/MGF-ID1/CsoDIAq-file1_ID1_corrected_peptideFDR.csv')
peptides = peptides.drop_duplicates(subset='peptide', keep='first').reset_index(drop=True)
print(len(peptides))
proteins = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/MGF-ID1/CsoDIAq-file1_ID1_corrected_proteinFDR.csv')
proteins = proteins[proteins['uniquePeptide']==1]
#proteins = proteins.drop_duplicates(subset='leadingProtein', keep='first').reset_index(drop=True)
prots = set()
for index, row in proteins.iterrows():
    prots.update(row['leadingProtein'].split('/'))
print(len(prots))
#'''
'''
# for determining the intensities of an mzxml file - I'm finding they're either all the same or has a lot of zeroes.
count = 0
#with mzxml.read('C:/Users/ccranney/Desktop/wiffFiles/mzxml-round1/18302_REP2_500ng_HumanLysate_SWATH_2.mzxml', use_index=True) as spectra:
with mzxml.read('C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1.mzxml', use_index=True) as spectra:
    for spec in spectra:
        if 'precursorMz' not in spec: continue
        count += 1
        if count > 10: break
        print(list(spec['intensity array']))
#'''
'''
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0.txt',sep='\t')
spectra = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy_proteinsAdded.mgf')
pepDict = {}
count = 0
for spec in spectra:
    count += 1
    if count % 1000 == 0: print(count)
    pepDict[spec['params']['seq']] = spec['params']['protein']

proteins = []
for index, row in df.iterrows():
    proteins.append(pepDict[row['Peptide']])

df['protein'] = proteins
df.to_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded.csv', index=False)
#'''
'''
# Adding MSPLIT proteins to final file
inFile = 'C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded.csv'
specFile = 'C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded_spectralFDR.csv'
pepFile = 'C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded_peptideFDR.csv'
protFile = 'C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded_proteinFDR.csv'
cbf.write_csodiaq_fdr_outputs(inFile, specFile, pepFile, protFile)
#'''
'''
# for comparing msplit with csodiaq
msplit = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_qehfx_09_30_0_proteinsAdded_peptideFDR.csv')
#msplit = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_ID1_02-10-0_temp.csv')
msplit['ID'] = list(zip(msplit['Scan#'].tolist(),msplit['peptide'].tolist(),msplit['z.1'].tolist()))
#msplit = msplit[msplit['uniquePeptide']==1]
#msplit = msplit.drop_duplicates(subset='leadingProtein', keep='first').reset_index(drop=True)
#msplit.set_index('ID',inplace=True)
#msplit = msplit[msplit['protein'].str.contains('DECOY')]


csodiaq = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/WIFF files/CsoDIAq-file1_01_qehfx_lab_SA_R1_corrected_peptideFDR.csv')
#csodiaq = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/WIFF files/CsoDIAq-file1_01_qehfx_lab_SA_R1_corrected_peptideFDR.csv')
csodiaq['ID'] = list(zip(csodiaq['scan'].tolist(),csodiaq['peptide'].tolist(),csodiaq['zLIB'].tolist()))
#csodiaq = csodiaq[csodiaq['uniquePeptide']==1]
#csodiaq = csodiaq.drop_duplicates(subset='leadingProtein', keep='first').reset_index(drop=True)
#csodiaq.set_index('ID',inplace=True)
#csodiaq = csodiaq[csodiaq['protein'].str.contains('DECOY')]

#mValues = set()
#cValues = set()
#for proteinGroup in msplit['leadingProtein']:
#    proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
#    mValues.update([p for p in proteins if 'DECOY' not in p])
#for proteinGroup in csodiaq['leadingProtein']:
#    proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
#    cValues.update([p for p in proteins if 'DECOY' not in p])


mValues = set(msplit['peptide'])
cValues = set(csodiaq['peptide'])

pyplot.figure(figsize=(12,12))
out = venn2([set(cValues),set(mValues)], set_labels = ["CsoDIAq", "MSPLIT-DIA"])
for text in out.set_labels:
    text.set_fontsize(32)
for text in out.subset_labels:
    text.set_fontsize(24)
pyplot.show()

#inM = set(msplit['ID']) - set(csodiaq['ID'])
#inC = set(csodiaq['ID']) - set(msplit['ID'])
#both = set(csodiaq['ID']).intersection(set(msplit['ID']))

#'''
'''
print(len(msplit['ID']))
print('\n')
print(len(inM))
print('\n')
print(len(csodiaq['ID']))
print('\n')
print(len(inC))
print('\n')
print(len(both))

msplit = msplit[msplit['ID'].isin(both)].sort_values('ID').reset_index(drop=True)
csodiaq = csodiaq[csodiaq['ID'].isin(both)].sort_values('ID').reset_index(drop=True)

for i in range(10):
    print(msplit.loc[i])
    print(csodiaq.loc[i])
    print('\n')
#'''
'''
oldMGF = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf')
protMGF = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy_proteinsAdded.mgf')

count = 0
for spec in oldMGF:
    count += 1
    if count % 5 ==0: break#print(count)
    spec2 = protMGF.get_by_id(spec['params']['title'])
    print(spec2)
    if sorted(spec['m/z array']) != sorted(spec2['m/z array']):
        print('m/z array')
        print(len(spec['m/z array']))
        print(len(spec2['m/z array']))
        print('\n')
    if sorted(spec['intensity array']) != sorted(spec2['intensity array']):
        print('intensity array')
        print(len(spec['intensity array']))
        print(len(spec2['intensity array']))
        print('\n')

    if len(spec['charge']) != len(spec2['charge']):
        print('charge')
        print(len(spec['charge']))
        print(len(spec2['charge']))
        print('\n')
    if len(spec['mask']) != len(spec2['mask']):
        print('mask')
        print(len(spec['mask']))
        print(len(spec2['mask']))
        print('\n')

    if spec['params']['title'] != spec2['params']['title']:
        print(spec['params']['title'])
        print(spec2['params']['title'])
        print('\n')
    if spec['params']['charge'] != spec2['params']['charge']:
        print(spec['params']['charge'])
        print(spec2['params']['charge'])
        print('\n')
    if spec['params']['pepmass'] != spec2['params']['pepmass']:
        print(spec['params']['pepmass'])
        print(spec2['params']['pepmass'])
        print('\n')
    if spec['params']['seq'] != spec2['params']['seq']:
        print(spec['params']['seq'])
        print(spec2['params']['seq'])
        print('\n')
    if spec['params']['scan'] != spec2['params']['scan']:
        print(spec['params']['scan'])
        print(spec2['params']['scan'])
        print('\n')
#'''
'''
# simple script for seeing a few scans in mzxml files

#with mzxml.read('C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1_large.mzXML', use_index=True) as spectra:
with mzxml.read('C:/Users/ccranney/Desktop/wiffFiles/K562_10ug_DDA_Top100_r02-K562_Sample10.mzXML', use_index=True) as spectra:
    count = 1
    windows = defaultdict(list)
    for spec in spectra:
        print(list(spec['intensity array']))
        if 'precursorMz' in spec: count += 1; windows[spec['precursorMz'][0]['windowWideness']].append(spec['num'])
        if count % 10 == 0: break
    #for key, value in windows.items():
        #print(key)
        #print(value)
#'''
'''
# simple script for seeing the first few lines of a csv data file
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/docker-shared/consensus.tsv', sep='\t')
#df.sort_values('MaCC_Score', ascending=False, inplace=True)
print(df.head(200))
print(df.columns)
#'''
'''
# checking if there are protein-marked decoys is the pan human library
#df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/docker-shared/human.tsv',sep='\t')
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/lib_tsv.tsv',sep='\t')
#print(df.head(10))
#df = df[df['ProteinName'].str.contains('DECOY')]
df = df[df['decoy']==1]
#df = df.head(100)
#df.to_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/lib_tsv_abbrev.csv')
print(len(df))
#'''
'''
# checking hov FullUniMode differs from PeptideSequence
diffs = set()
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/lib_tsv.tsv',sep='\t')
for index, row in df.iterrows():
    if row['PeptideSequence'] != row['FullUniModPeptideName']:
        diffs.add((row['PeptideSequence'],row['FullUniModPeptideName']))

print(len(diffs))
print(sorted(diffs)[:20])
#'''
'''
# Counting peptides in targetted re-analysis files
files = [
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-30.0.csv',
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-40.0.csv',
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-50.0.csv',
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-60.0.csv',
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-70.0.csv',
    'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CsoDIAq-file1_20210222_DISPA_hela_03_corrected_mostIntenseTargs_-80.0.csv',
]
count = 0
for file in files:
    print(file)
    df = pd.read_csv(file)
    peptides = df['Compound']
    for pep in peptides:
        peps = pep.split('/')
        count += len(peps)
print(count)
#'''
'''
# convert tsv to csv
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/docker-shared/qehfx.tsv', sep='\t')
df.to_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/docker-shared/qehfx.csv')
#'''
'''
#counting MGF spectra in a file
spectra = mgf.read('C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf')
count = 0
for spec in spectra: count += 1
print(count)
'''
'''
#testing table manipulations
with bz2.BZ2File('C:/Users/ccranney/Desktop/Caleb_Files/data/output/compressTest/compressed_processObj_canDeleteAfter/window0', 'rb') as pickleFile: data = pickle.load(pickleFile)
print(data.head(10))

def calc_macc(df):
    lib = df['libIntensity']
    que = df['queIntensity']
    AB = lib.multiply(que).sum()
    A = lib.pow(2).sum()
    B = que.pow(2).sum()
    cosine = cbf.cosine_similarity([AB, A, B])
    return (len(df)**(1/5))*cosine
#test = data.groupby(['libID','queID']).apply(calc_macc)
#print(data.head(10))
tempDict = {}
grouped = data.groupby(['libID','queID'])
print(grouped.ngroups)
data2 = grouped.filter(lambda x: len(x.index) > 2.)
grouped = data2.groupby(['libID','queID'])
tempList = [(21, 42743),
(21, 44234),
(21, 44376),
(21, 44447),
(21, 47003),
(21, 64256),
(21, 64611),
(21, 64966),
(21, 71782),
(21, 71853),
(21, 71924),
(21, 71995),
(21, 72066),
(21, 72208),
(21, 72279),
(21, 72350),
(21, 72421),
(21, 72705),
(21, 72776),
(21, 72847),
(21, 72989),
(21, 73060),
(21, 73131),
(21, 73202),
(21, 73273),
(21, 73344),
(21, 76965),
(21, 77036),
(21, 77107),
(21, 77178),
(21, 77249),
(21, 77320),
(21, 77391),
(21, 77604),
(21, 77817),
(21, 85272),
(41, 85201),
(60, 52825),
(60, 53180),
(60, 65250),
(60, 65321),
(60, 65392),
(60, 69439),
(60, 69510),
(60, 78598),
(60, 78669),
(60, 78811),
(84, 37276),
(84, 37347),
(84, 37489),
(84, 37560),
(86, 36424),
(86, 42388),
(86, 42459),
(86, 42530),
(86, 42601),
(86, 42814),
(86, 42956),
(86, 43382),
(86, 43524),
(86, 49204),
(86, 49346),
(86, 50908),
(86, 59925),
(96, 34862),
(96, 34933),
(96, 35004),
(96, 35075),
(96, 35146),
(96, 37915),
(96, 42246),
(96, 42672),
(96, 42743),
(96, 42814),
(96, 44589),
(96, 44660),
(96, 44731),
(96, 44802),
(96, 44873),
(96, 44944),
(96, 45015),
(96, 45086),
(96, 45228),
(96, 45370),
(96, 45441),
(96, 45583),
(96, 52896),
(96, 56446),
(96, 56517),
(96, 56588),
(96, 56659),
(96, 61203),
(96, 65037),
(96, 78740),
(96, 78953),
(96, 79024),
(96, 79166),
(96, 79237)]
for x in tempList: tempDict[x]=0
#data2['MaCC'] = grouped.apply(lambda x,d: d[x.name])
data2['MaCC'] = data2.groupby(['libID','queID']).apply(lambda x,d: print(d[x.name]),tempDict)
#d[tuple(x.name)]
#count = 0
#for name, group in grouped:
#    print(type(name))
#    count += 1
#    if count==10: break
#print(grouped.ngroups)
#data2['MaCC'] = grouped.apply(calc_macc)
#for name, group in grouped:
#    macc = calc_mass(group)
#print(data2.head(10))
'''
'''
# Writing/testing reduction function
def reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    count = 1
    #maccScores = []
    AB = matchLibIntensities[0]*matchQueIntensities[0]
    A = matchLibIntensities[0]**2
    B = matchQueIntensities[0]**2
    returnLibTags = []
    returnLibIntensities = []
    returnQueTags = []
    returnQueIntensities = []
    returnPpmMatchs = []
    returnMaccScores = []
    length = len(matchLibTags)
    for i in range(1,length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            if count > 2:
                cosine = cbf.cosine_similarity(AB, A, B)
                macc = (count**(1/5))*cosine
                magnitude = (A**0.5) * (B**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
                returnLibTags.extend(matchLibTags[i-count:i])
                returnLibIntensities.extend(matchLibIntensities[i-count:i])
                returnQueTags.extend(matchQueTags[i-count:i])
                returnQueIntensities.extend(matchQueIntensities[i-count:i])
                returnPpmMatchs.extend(ppmMatches[i-count:i])
                returnMaccScores.extend([macc]*count)
            count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
            AB = matchLibIntensities[i]*matchQueIntensities[i]
            A = matchLibIntensities[i]**2
            B = matchQueIntensities[i]**2
        else:
            AB += matchLibIntensities[i]*matchQueIntensities[i]
            A += matchLibIntensities[i]**2
            B += matchQueIntensities[i]**2
            count += 1
    if count > 2:
        cosine = cbf.cosine_similarity(AB, A, B)
        macc = (count**(1/5))*cosine
        returnLibTags.extend(matchLibTags[length-count:])
        returnLibIntensities.extend(matchLibIntensities[length-count:])
        returnQueTags.extend(matchQueTags[length-count:])
        returnQueIntensities.extend(matchQueIntensities[length-count:])
        returnPpmMatchs.extend(ppmMatches[length-count:])
        returnMaccScores.extend([macc]*count)
        count = 1

    return returnLibTags, returnLibIntensities, returnQueTags, returnQueIntensities, returnPpmMatchs, returnMaccScores

libTags = np.repeat(np.arange(0,10),3)
libIntensities = np.arange(1.0,31.0)
queTags = np.repeat(np.arange(0,5),6)
queIntensities = np.arange(31.0,61.0)
ppmMatches = np.arange(61.0,91.0)
#print(libTags)
#print(libIntensities)
#print(queTags)
#print(queIntensities)
#print(ppmMatches)


matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, maccScores = reduce_final_df(libTags, libIntensities, queTags, queIntensities, ppmMatches)

#print(len(libTags))
#print(len(matchLibTags))
#print(len(set(maccScores)))
#print(libTags)
#print(matchLibTags)
print(maccScores)
'''
'''
# Testing last part of "compression" attempt
@njit
def match_score_decoys(matchLibTags, matchQueTags, maccScores, decoys):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    returnMaccs = [maccScores[0]]
    returnDecoys = [decoys[0]]
    length = len(matchLibTags)
    for i in range(1,length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            returnMaccs.append(maccScores[i])
            returnDecoys.append(decoys[i])
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
    #return returnLibTags, returnLibIntensities, returnQueTags, returnQueIntensities, returnPpmMatchs, returnMaccScores
    return returnMaccs, returnDecoys

@njit
def fdr_calculation2(maccs, decoys): #***NOTE***make return type
    # initializing the two return values at 0
    numDecoys = 0

    # for every row in the dataframe
    count = 0
    for i in range(len(maccs)):
        #print(str(maccs[i])+':'+str(decoys[i]))
        if decoys[i]: numDecoys += 1
        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(count+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if count < 1/0.01: return -1
            return maccs[i-1]
        count += 1

    return maccs[-1]

def collect_ppm_values(ppmMatches, maccScores, maccCutoff):
    ppms = []
    length = len(ppmMatches)
    for i in range(1,length):
        if maccScores[i] >= maccCutoff: ppms.append(ppmMatches[i])
    #return returnLibTags, returnLibIntensities, returnQueTags, returnQueIntensities, returnPpmMatchs, returnMaccScores
    return ppms

#test = np.arange(10,0,-1)
#print(type(test.argsort()))
#print(type(np.where((test>2)*(test<8))[0]))
with bz2.BZ2File('C:/Users/ccranney/Desktop/Caleb_Files/data/output/CompressTest/compressed_processObj_canDeleteAfter/windowallFinals', 'rb') as pickleFile: data = pickle.load(pickleFile)
finalLibTags, finalLibIntensities, finalQueTags, finalQueIntensities, finalPpmMatches, finalMaccScores, finalDecoys = data
print(len(finalDecoys))
print(timer())
#temp1, temp2 = match_score_decoys(np.array([],dtype=int),np.array([],dtype=int),np.array([],dtype=float),np.array([],dtype=int))
#print(timer())
maccs, decoys = match_score_decoys(finalLibTags, finalQueTags, finalMaccScores, finalDecoys)
print(timer())
print(len(maccs))
print(len(set(finalMaccScores)))
print(len(set(maccs)))
maccs = np.array(maccs)
decoys = np.array(decoys)
i1 = (-maccs).argsort()
maccs = maccs[i1]
decoys = decoys[i1]
maccCutoff = fdr_calculation2(maccs, decoys)
print(maccCutoff)
cutoff = 1.4357304070871908
ppms = collect_ppm_values(finalPpmMatches, finalMaccScores, maccCutoff)
offset, tolerance = cbf.find_offset_tol(ppms, 'C:/Users/ccranney/Desktop/Caleb_Files/data/output/CompressTest/CsoDIAq-file1_01_qehfx_lab_SA_R1_histogram.png', stdev=0)
print(offset,tolerance)
lowend = offset-tolerance
highend = offset+tolerance
ppmIndices = np.where((finalPpmMatches>lowend)*(finalPpmMatches<highend))[0]
print(len(ppmIndices))
tempTags = finalLibTags[ppmIndices]
print(len(tempTags))
print(ppmIndices[:10])
print(tempTags[:10])
corLibTags, corLibIntensities, corQueTags, corQueIntensities, corDecoys = [x[ppmIndices] for x in [finalLibTags, finalLibIntensities, finalQueTags, finalQueIntensities, finalDecoys]]
corLibTags, corLibIntensities, corQueTags, corQueIntensities, corDecoys, corMaccScores = cbf.reduce_final_df(corLibTags, corLibIntensities, corQueTags, corQueIntensities, corDecoys)
maccs, decoys = match_score_decoys(corLibTags, corQueTags, corMaccScores, corDecoys)
i1 = (-maccs).argsort()
maccs = maccs[i1]
decoys = decoys[i1]
maccCutoff = fdr_calculation2(maccs, decoys)
'''

t1 = np.arange(10)
t2 = np.arange(10,20)
t3 = np.arange(20,30)
compressed_array = io.BytesIO()
compressed_array2 = io.BytesIO()

#np.savez_compressed(compressed_array, t1=t1, t2=t2, t3=t3)
np.savez_compressed(compressed_array, t1=t1)
np.savez_compressed(compressed_array2, t1=t2)
#np.savez_compressed(compressed_array, t3=t3)

compList = [compressed_array, compressed_array2]
for x in compList:
    x.seek(0)
    decompressed_array = np.load(x)
    print(decompressed_array['t1'])
    #print(decompressed_array['t2'])
#print(decompressed_array['t3'])




#java -Xmx2500M -cp C:/Users/ccranney/Desktop/Caleb_Files/MSPLIT-DIAv1.0/MSPLIT-DIAv02102015.jar org.Spectrums.SWATHMSPLITSearch 02 10 0 C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1.mzXML C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_HeLa_02-10-0.tsv
#docker run -it -v C:/Users/ccranney/Desktop/Caleb_Files/data/docker-shared:/data openswath/openswath
#TargetedFileConverter -in phl004_consensus_openms24.TraML -out consensus.tsv
#OpenSwathWorkflow -in 01_qehfx_lab_SA_R1.mzXML -tr consensus.tsv -out_features qehfx.featureXML
#'''
'''
                        decoyList = [idToDecoyDict[x] for x in libTags]
                        #libMzs = np.array([x[0] for x in pooledLibSpectra])
                        queMzs = np.array([x[0] for x in pooledQueSpectra])
                        tempTag = 'identify'
                        returns = initialize_return_values(tempTag)

                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys = spectra_peak_comparison(libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, decoyList, ppmTol, ppmYOffset)
                        #data = pd.DataFrame({'libID':matchLibTags, 'queID':matchQueTags, 'libIntensity':matchLibIntensities, 'queIntensity':matchQueIntensities, 'ppmDiff':ppmMatches})
                        #data = data.groupby(['libID','queID']).filter(lambda x: len(x.index) > 2.)
                        #print(len(data))
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys = [np.array(x) for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys]]
                        i1 = matchQueTags.argsort(kind='mergesort')
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys = [x[i1] for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys]]
                        i2 = matchLibTags.argsort(kind='mergesort')
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys = [x[i2] for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys]]
                        #reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches)
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys, maccScores = reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys)
                        test = [idToDecoyDict[x] for x in matchLibTags]
def reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, matchDecoys):
    curLibTag = matchLibTags[0]
    curQueTag = matchQueTags[0]
    count = 1
    AB = matchLibIntensities[0]*matchQueIntensities[0]
    A = matchLibIntensities[0]**2
    B = matchQueIntensities[0]**2
    returnLibTags = []
    returnLibIntensities = []
    returnQueTags = []
    returnQueIntensities = []
    returnPpmMatches = []
    returnDecoys = []
    returnMaccScores = []
    length = len(matchLibTags)
    for i in range(1,length):
        if matchLibTags[i] != curLibTag or matchQueTags[i] != curQueTag:
            if count > 2:
                test = 0
                cosine = cosine_similarity(AB, A, B)
                macc = (count**(1/5))*cosine
                magnitude = (A**0.5) * (B**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
                returnLibTags.extend(matchLibTags[i-count:i])
                returnLibIntensities.extend(matchLibIntensities[i-count:i])
                returnQueTags.extend(matchQueTags[i-count:i])
                returnQueIntensities.extend(matchQueIntensities[i-count:i])
                returnPpmMatches.extend(ppmMatches[i-count:i])
                returnDecoys.extend(matchDecoys[i-count:i])
                returnMaccScores.extend([macc]*count)

            count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
            AB = matchLibIntensities[i]*matchQueIntensities[i]
            A = matchLibIntensities[i]**2
            B = matchQueIntensities[i]**2
        else:
            AB += matchLibIntensities[i]*matchQueIntensities[i]
            A += matchLibIntensities[i]**2
            B += matchQueIntensities[i]**2
            count += 1
    if count > 2:
        test = 0
        cosine = cosine_similarity(AB, A, B)
        macc = (count**(1/5))*cosine
        returnLibTags.extend(matchLibTags[length-count:])
        returnLibIntensities.extend(matchLibIntensities[length-count:])
        returnQueTags.extend(matchQueTags[length-count:])
        returnQueIntensities.extend(matchQueIntensities[length-count:])
        returnPpmMatches.extend(ppmMatches[length-count:])
        returnDecoys.extend(matchDecoys[length-count:])
        returnMaccScores.extend([macc]*count)

    #return returnLibTags, returnLibIntensities, returnQueTags, returnQueIntensities, returnPpmMatchs, returnMaccScores
    return returnLibTags, returnLibIntensities, returnQueTags, returnQueIntensities, returnPpmMatches, returnDecoys, returnMaccScores

####################################
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches = spectra_peak_comparison(libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, ppmTol, ppmYOffset)

                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches = [np.array(x) for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches]]
                        i1 = matchQueTags.argsort(kind='mergesort')
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches = [x[i1] for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches]]
                        i2 = matchLibTags.argsort(kind='mergesort')
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches = [list(x[i2]) for x in [matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches]]
                        #reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches)
                        matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, maccScores = reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches)
                        print(len(set(maccScores)))
                        #data = pd.DataFrame({'libID':matchLibTags, 'queID':matchQueTags, 'libIntensity':matchLibIntensities, 'queIntensity':matchQueIntensities, 'ppmDiff':ppmMatches})
                        #data = pd.concat(dfs)
                        #data = data.groupby(['libID','queID']).filter(lambda x: len(x.index) > 2.)
                        #data = data.groupby(['libID','queID']).apply(add_macc_column, maccDict)

                        #data, NAlist = reduce_mem_usage(data)
                        #allDfs.append(data)
                        #with bz2.BZ2File(pickleDir+pickleHeader+str(count), 'w') as pickleFile: pickle.dump(data, pickleFile)

                        pooledQueSpectra.clear()
                        del returns


                count += 1
                if count % printCutoff == 0:
                    time = timer()
                    #print('\nNumber of Pooled Experimental Spectra Analyzed: ' + str(count))
                    #print('Number of Spectra in Current Pooled Spectra: ' + str(len(scans)))
                    #print('Time Since Last Checkpoint: ' + str(round(time-prevtime,2)) + ' Seconds', flush=True)
                    print(round(time-prevtime,2),flush=True)
                    prevtime = time


    # Prints the final number of experimental spectra analyzed.
    print('Total Time (seconds): ' + str(timer()))
    print('Count: '+str(count),flush=True)
    if tag=='identify': return ppmList
#''
#@njit
def reduce_final_df(matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches):
    curLibTag = matchLibTags[-1]
    curQueTag = matchQueTags[-1]
    count = 1
    maccScores = []
    AB = matchLibIntensities[-1]*matchQueIntensities[-1]
    A = matchLibIntensities[-1]**2
    B = matchQueIntensities[-1]**2
    for i in reversed(range(len(matchLibTags)-1)):
        if matchLibTags[i] != curLibTag and matchQueTags[i] != curQueTag:
            if count == 1:
                matchLibTags.pop(i-1)
                matchLibIntensities.pop(i-1)
                matchQueTags.pop(i-1)
                matchQueIntensities.pop(i-1)
                ppmMatches.pop(i-1)
            else:
                cosine = cosine_similarity([AB, A, B])
                macc = (count**(1/5))*cosine
                maccScores.extend([macc]*count)
                count = 1
            curLibTag = matchLibTags[i]
            curQueTag = matchQueTags[i]
            AB = matchLibIntensities[i]*matchQueIntensities[i]
            A = matchLibIntensities[i]**2
            B = matchQueIntensities[i]**2
        else:
            AB += matchLibIntensities[-1]*matchQueIntensities[-1]
            A += matchLibIntensities[-1]**2
            B += matchQueIntensities[-1]**2
            count += 1

    return matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches, maccScores[::-1]

def pooled_spectra_analysis(expSpectraFile, outFile, lib, ppmTol, ppmYOffset, queryPooling, spectraKeys=None):

    pickleDir = '/'.join(outFile.split('/')[:-1]) + '/compressed_processObj_canDeleteAfter/'
    if not os.path.exists(pickleDir):
        os.mkdir(pickleDir)
    pickleHeader = 'window'

    # query data file is loaded
    with mzxml.read(expSpectraFile, use_index=True) as spectra:


        # query data is looped over and scans are grouped by mz windows for future pooling.
        #  Additionally, for the second condition, various variables are saved for future reference.
        queScanDict = defaultdict(list)
        queValDict = defaultdict(dict)
        spectralCount = 0

        for spec in spectra:
            if 'precursorMz' not in spec: continue
            queScanDict[spec['precursorMz'][0]['precursorMz'],spec['precursorMz'][0]['windowWideness']].append(spec['num'])
            #if scan=='1199' or scan=='1200': print((spec['precursorMz'][0]['precursorMz'],spec['precursorMz'][0]['windowWideness']))
            peaksCount = spec['peaksCount']
            if 'compensationVoltage' in spec: CV = spec['compensationVoltage']
            else: CV = ''
            queValDict[spec['num']]['peaksCount'] = peaksCount
            queValDict[spec['num']]['CV'] = CV
            spectralCount += 1

        print('Number of Unpooled MS/MS Query Spectra: ' + str(spectralCount))
        print('Number of Pooled MS/MS Query Spectra/Mz Windows: ' + str(len(queScanDict)),flush=True)

        # To enhance the print experience, status prints will be given at intervals tailored to the number of identified windows.
        #  example: if there are 1-99 pooled query spectra, print statements are made after every pooled query spectra analysis is complete.
        #           if there are 100-999, print after every 10 pooled spectra. And so on.
        printCutoff = 100
        while printCutoff < len(queScanDict): printCutoff*=10
        printCutoff /= 100

        # 'lib' dictionary keys are kept as a separate list in this analysis. Note that they are sorted by precursor m/z.
        allLibKeys = lib.keys()
        allLibKeys = sorted(allLibKeys)

        # outfile is opened in advance so results can be written directly to the file as they are produced (mitigating memory use)
        # The second condition of this function returns a list of PPM differences that can be used for correction. Initialized here
        ppmList = []

        # Count variable keeps track of the number of query spectra that have been analyzed for time tracking purposes.
        count = 0

        # Library keys were saved as an integer to save on time and simplify other parts of the algorithm. I believe the purpose is now defunct, but it's harmless, so I'm keeping it.
        idToKeyDict = {}
        for key in allLibKeys: idToKeyDict[lib[key]['ID']] = key

        # tracking time for print statements.
        prevtime = timer()

        print('Enter Pooled Spectra Analysis:')
        print(str(timedelta(seconds=prevtime)), flush=True)

        columns=['libID','queID','libIntensity','queIntensity','ppmDiff']
        allDfs = []
        allData = []
        maccDecoys = []
        decs = []
        # looping through all the windows that the query data corresponds to
        for precMz_win, scans in queScanDict.items():

            # Determining all library spectra that should be pooled for this particular query data window
            top_mz = precMz_win[0] + precMz_win[1] / 2
            bottom_mz = precMz_win[0] - precMz_win[1] / 2
            libKeys = lib_mz_match_query_window( top_mz, bottom_mz, allLibKeys )
            if len(libKeys) == 0: continue
            pooledLibSpectra = pool_lib_spectra(lib, libKeys)

            # begin pooling query spectra
            pooledQueSpectra = []
            returns = initialize_return_values('identify')
            for i in range(len(scans)):

                # adding each scan in the window to the pooled spectra
                scan = scans[i]
                spec = spectra.get_by_id(scan)
                intensity = [x**0.5 for x in spec['intensity array']]
                peakIDs = [int(scan) for x in range(spec['peaksCount'])]
                pooledQueSpectra += list(zip(spec['m/z array'],intensity,peakIDs))

                # to reduce memory use for particularly large files, the user can limit the number of query spectra that are pooled. That's what this conditional statement takes care of.
                if (i % queryPooling == 0 and i!=0) or i == len(scans)-1:
                    pooledQueSpectra.sort()

                    # each peak in the pooled library and query spectra is compared, and the necessary data is extracted from matching peaks (as determined by the ppm tolerance)
                    libMzs = np.array([x[0] for x in pooledLibSpectra])
                    queMzs = np.array([x[0] for x in pooledQueSpectra])

                    pMatch, jMatch, ppmMatch = spectra_peak_comparison(libMzs, queMzs, ppmTol, ppmYOffset)
                    for i in range(len(pMatch)):
                        libPeak = pooledLibSpectra[pMatch[i]]
                        quePeak = pooledQueSpectra[jMatch[i]]

                        update_return_values(returns, pooledLibSpectra[pMatch[i]], pooledQueSpectra[jMatch[i]], pMatch[i], jMatch[i], ppmMatch[i], 'identify')
                        #data.append([libPeak[2], quePeak[2], libPeak[1], quePeak[1], ppmMatch[i]])
                    # output values are saved to the output file

                    pooledQueSpectra.clear()

            cosDict, countDict, ionDict, ppmDict = returns
            data = []
            maccs = []

            for key, value in countDict.items():
                if value > 2:
                    macc = (value**(1/5))*cosine_similarity(cosDict[key])
                    #df = pd.DataFrame(ppmDict[key], columns=columns)
                    #df['MaCC'] = pd.Series(macc, index=df.index)
                    #dfs.append(df)
                    #data += ppmDict[key]
                    #maccDict[key] = macc
                    #data.extend(ppmDict[key])
                    #allData += ppmDict[key]
                    maccs += [macc] * value

                    #libKey = idToKeyDict[key[0]]
                    decoy = 0
                    if 'DECOY' in lib[idToKeyDict[key[0]]]['ProteinName']: decoy = 1
                    #decs += [decoy] * value
                    maccDecoys.append((macc,decoy))

            #data = pd.concat(dfs)
            data = pd.DataFrame(data, columns=['libID','queID','libIntensity','queIntensity','ppmDiff'])
            data['MaCC_Score'] = maccs
            #data = data.groupby(['libID','queID']).filter(lambda x: len(x.index) > 2.)
            #data = data.groupby(['libID','queID']).apply(add_macc_column, maccDict)

            data, NAlist = reduce_mem_usage(data)
            #allDfs.append(data)
            #with bz2.BZ2File(pickleDir+pickleHeader+str(count), 'w') as pickleFile: pickle.dump(data, pickleFile)

            #compressedWindows.append(zlib.compress(pickle.dumps(data)))
            #print('size: '+str(sys.getsizeof(zlib.compress(pickle.dumps(data),zlib.Z_BEST_COMPRESSION))), flush=True)
            #print('size: '+str(sys.getsizeof(data)), flush=True)
            #return []
            #compressedWindows.append(data)
            # print statements for the user to track progress.
            count += 1
            #if count % 1000 == 0:
            if count % printCutoff == 0:
                time = timer()
                #print('\nNumber of Pooled Experimental Spectra Analyzed: ' + str(count))
                #print('Number of Spectra in Current Pooled Spectra: ' + str(len(scans)))
                #print('Time Since Last Checkpoint: ' + str(round(time-prevtime,2)) + ' Seconds', flush=True)
                print(round(time-prevtime,2),flush=True)
                prevtime = time

    pickle.dump(maccDict, open('C:/Users/ccranney/Desktop/Caleb_Files/data/output/CompressTest/maccDict.p', 'wb'))

    #lib = pickle.load(open(args['outDirectory']+'mgf_lib.p', 'rb'))
    #print('size: '+str(sys.getsizeof(compressedWindows)))
    # Prints the final number of experimental spectra analyzed.
    print('Total Time (seconds): ' + str(timer()))
    print('Count: '+str(count),flush=True)
    return ppmList
#'''
