import pandas as pd
from pyteomics import mzxml, mgf, mass
from collections import defaultdict
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


print('\n'*30)


#'''
scanCount = defaultdict(int)
msplit = pd.read_csv('C:/Users/ccranney/Desktop/Caleb_Files/data/output/HeLa-50ppm/MSPLIT_HeLa.txt', sep='\t')
scans = msplit['Scan#']
for scan in scans:
    scanCount[scan]+=1

#hits, decoys = fdr_calculation(msplit)
#print(hits, decoys)

import operator
'''
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

file = 'C:/Users/ccranney/Desktop/Caleb_Files/data/18302_REP2_500ng_HumanLysate_SWATH_2 (2)..mzXML'
with open(file, 'rb') as f:
    for _ in range(100): # first 10 lines
        print(f.readline())


#java -Xmx2500M -cp C:/Users/ccranney/Desktop/Caleb_Files/MSPLIT-DIAv1.0/MSPLIT-DIAv02102015.jar org.Spectrums.SWATHMSPLITSearch 02 10 0 C:/Users/ccranney/Desktop/Caleb_Files/data/HeLa_160min_DIA_106win_1.mzXML C:/Users/ccranney/Desktop/Caleb_Files/data/human.faims.fixed.decoy.mgf C:/Users/ccranney/Desktop/Caleb_Files/data/output/msplit_HeLa_02-10-0.tsv
