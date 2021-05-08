'''
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
'''
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

df = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/phl004_canonical_sall_osw_decoys.csv', sep='\t')
print(df.head(10))
