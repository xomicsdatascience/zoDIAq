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
msplit['ID'] = list(zip(msplit['Scan#'].tolist(),msplit['peptide'].tolist(),msplit['z.1'].tolist()))
#msplit = msplit[msplit['uniquePeptide']==1]
#msplit = msplit.drop_duplicates(subset='leadingProtein', keep='first').reset_index(drop=True)
#msplit.set_index('ID',inplace=True)
#msplit = msplit[msplit['protein'].str.contains('DECOY')]


csodiaq = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/WIFF files/CsoDIAq-file1_01_qehfx_lab_SA_R1_corrected_peptideFDR.csv')
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
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/WIFF files/CsoDIAq-file1_HeLa_160min_DIA_106win_1_large_delete.csv')
df.sort_values('MaCC_Score', ascending=False, inplace=True)
print(df.head(200))
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
# checking hov FullUniMode differs from PeptideSequence
diffs = set()
df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/lib_tsv.tsv',sep='\t')
for index, row in df.iterrows():
    if row['PeptideSequence'] != row['FullUniModPeptideName']:
        diffs.add((row['PeptideSequence'],row['FullUniModPeptideName']))

print(len(diffs))
print(sorted(diffs)[:20])


#java -Xmx2500M -cp C:/Users/ccranney/Desktop/Caleb_Files/MSPLIT-DIAv1.0/MSPLIT-DIAv02102015.jar org.Spectrums.SWATHMSPLITSearch 02 10 0 C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1.mzXML C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_HeLa_02-10-0.tsv
