from bisect import bisect_left
from Bio import SeqIO
import re
from pyteomics import mgf
from . import spectra_matcher_functions as smf


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

def clean_mgf_file(mgfFile, fasta, ions=False):
    spectra = mgf.read(mgfFile)
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
        #if count % 40 == 0: break
        if ions:
            mzValues = return_frag_mzs(spec['params']['seq'],1)
            peaks = list(tuple(zip(spec['m/z array'],spec['intensity array'])))
            for i in range(len(peaks)-1,-1,-1):
                if smf.approx_list(peaks[i][0],mzValues)==-1: peaks.pop(i)
            if len(peaks)==0: continue
            peaks.sort(key=lambda x:x[0])
            spec['m/z array'],spec['intensity array'] = map(list,zip(*peaks))
        decoy = False
        if 'protein' in spec['params'] and 'DECOY' in spec['params']['protein']: decoy = True
        else:
            seq = re.sub(r'\+\d+\.\d+', '', spec['params']['seq'])
            listOfI = [m.start() for m in re.finditer(seq, longPep)]
            sorted_keys = sorted(fDict.keys())
            proteins = set()
            for i in listOfI:
                insertion_point = bisect_left(sorted_keys,i)
                if insertion_point==len(sorted_keys) or sorted_keys[insertion_point]!=i:
                    insertion_point-=1
                protein = fDict[sorted_keys[insertion_point]]
                proteins.add(fDict[sorted_keys[insertion_point]])
            if len(proteins)==0: proteins.add('protein_not_in_fasta_'+spec['params']['seq'])

        if decoy: proteins = ['DECOY_0_'+x for x in proteins]

        protein = str(len(proteins)) + '/' + '/'.join(sorted(proteins))
        if protein != '0/': spec['params']['protein'] = protein; pepCount += 1
        cleaned.append(spec)
        if count % 1000 == 0: print(count); print(pepCount); print(protein)


    cleanedFile = re.sub('(.*).mgf', r'\1_proteinsAdded.mgf', mgfFile)
    if ions: cleanedFile = re.sub('(.*).mgf', r'\1_YBionsOnly.mgf', cleanedFile)
    mgf.write(cleaned, cleanedFile)
