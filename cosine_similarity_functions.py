import pandas as pd
from pyteomics import mzxml
from bisect import bisect
import pickle
from sklearn.metrics.pairwise import cosine_similarity
from timeit import default_timer as timer
from datetime import timedelta
import csv

def dot(A,B):
    return (sum(a*b for a,b in zip(A,B)))

def cosine_similarity_manual(a,b):
    return dot(a,b) / ( (dot(a,a)**.5) * (dot(b,b)**.5) )

def cosine_between_rows(fullSpec,i1,i2):
    return cosine_similarity_manual(list(fullSpec.iloc[i1].values),list(fullSpec.iloc[i2].values))

def approx(x, y, tol=0.005):
    return abs(x-y) <= tol

def cosine_similarity_msplit(libSpectrum, libKeys, exp_mz, exp_i):
    libDict = {k: [0.0, 0.0, 0.0, 0] for k in libKeys} #product, query magnitude, library magnitude, count, #ion_count

    i, j = 0, 0

    while i < len(libSpectrum) and j < len(exp_mz):
        if not approx(libSpectrum[i][0],exp_mz[j]):
            if libSpectrum[i][0] > exp_mz[j]: j += 1; continue
            if libSpectrum[i][0] < exp_mz[j]: i += 1; continue
        p = i + 0
        while (p < len(libSpectrum)):
            libDict[libSpectrum[p][2]][0] += libSpectrum[p][1]*exp_i[j]
            libDict[libSpectrum[p][2]][1] += exp_i[j]**2
            libDict[libSpectrum[p][2]][2] += libSpectrum[p][1]**2
            libDict[libSpectrum[p][2]][3] += 1
            p += 1
            if p==len(libSpectrum) or not approx(libSpectrum[p][0], exp_mz[j]): break
        #ion_count += exp_i[j]
        j += 1
    return([calculateCosine(key, libDict[key]) for key in libDict.keys()])

def calculateCosine(key, row):
    product = row[0]
    magnitude = (row[1]**0.5) * (row[2]**0.5)
    return [key, (product / magnitude if magnitude else 0), row[3]]

def libWindowMassMatch( spec, sorted_lib_keys ):
    temp = sorted_lib_keys[:]
    top_mz = spec['precursorMz'][0]['precursorMz'] + spec['precursorMz'][0]['windowWideness'] / 2
    bottom_mz = spec['precursorMz'][0]['precursorMz'] - spec['precursorMz'][0]['windowWideness'] / 2
    top_key = (top_mz, "z")
    bottom_key = (bottom_mz, "")
    i1 = bisect(temp, bottom_key)
    if i1 == len(temp)-1 and len(temp)!= 1: return []
    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
    temp.insert(i2, top_key)
    return temp[i1+1:i2]

def tramlFileConversionCSV(file_name):
    if file_name.endswith('.tsv'):
        lib_df = pd.read_csv(file_name, sep='\t')
    else:
        lib_df = pd.read_csv(file_name)
    print("Enter library dictionary upload: ")
    print(timedelta(seconds=timer()))
    lib_df = lib_df.loc[:, lib_df.columns.intersection(['PrecursorMz','FullUniModPeptideName','PrecursorCharge','ProductMz','LibraryIntensity','transition_group_id','ProteinName'])]
    lib_df['ID'] = list(zip(lib_df['PrecursorMz'].tolist(),lib_df['FullUniModPeptideName'].tolist(),lib_df['ProteinName'].tolist()))
    mz_dict = lib_df.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensity_dict = lib_df.groupby("ID")['LibraryIntensity'].apply(list).to_dict()
    lib_df.drop_duplicates(subset="ID",inplace=True)
    lib_df.set_index("ID", drop=True, inplace=True)
    lib = lib_df.to_dict(orient="index")
    for key in lib:
        #lib[key]['m/z array'], lib[key]['intensity array'] = (list(t) for t in zip(*sorted(zip(mz_dict[key], intensity_dict[key]))))
        mz, intensity = (list(t) for t in zip(*sorted(zip(mz_dict[key], intensity_dict[key]))))
        keyList = [key for i in range(len(mz))]
        lib[key]['Peaks'] = list(tuple(zip(mz,intensity,keyList)))
    return lib

def extractLibSpectra(lib, libKeys):
    finalList = []
    for key in libKeys: finalList += lib[key]['Peaks']
    return sorted(finalList)

def expSpectraAnalysis( expSpectraFile, outFile, lib ):
    columns = [
        'fileName', # Name of the experimental SWATH spectrum file
        'scan#', # Scan number, corresponding to scans in the experimental SWATH spectrum file
        'MzEXP', # precursor m/z for SWATH scan. Column 'windowWideness' corresponds to this value.
        'zEXP', # precursor charge for SWATH scan. Note that this is likely an estimate, as the exact charge of every compound in the scan is unknown.
        'peptide', # Peptide sequence for the library spectrum corresponding to this row.
        'protein', # Protein name the peptide corresponds to, also derived from the library spectrum corresponding to this row.
        'MzLIB', # precursor m/z for the library spectrum corresponding to this row.
        'zLIB', # precursor charge for the library spectrum corresponding to this row.
        'cosine', # cosine score comparing the library spectrum corresponding to this row with the experimental spectrum.
        'name', # Title - corresponds to the column ""
        '#Peak(Query)', # I assume the number of peaks in the experimental spectrum.
        '#Peaks(Match)', # I assume the number of peaks in the library spectrum.
        'shared', # I assume the number of peaks that matched between experimental spectrum/library spectrum
        'ionCount', # Sum of experimental spectrum intensities, excluding possible duplicates
        'CompensationVoltage', #the compensation voltage of the experimental SWATH spectrum
        'totalWindowWidth' # width of m/z that was captured in the experimental SWATH spectrum. Corresponds to MzSWATH
    ]
    with open(outFile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(columns)
        count = 0
        all_lib_keys = sorted(lib)

        time = timer()
        prevtime = time
        print("Index,TimeElapsed,NumPeaks,ExpPrecursorMz,NumLibrarySpectra")
        with mzxml.read(expSpectraFile) as spectra:
            for spec in spectra:
#                spec['m/z array'], spec['intensity array'] = querySpectrumFilter(spec['m/z array'],spec['intensity array'],3,25.0)
                count += 1
#                time = timer()


                lib_keys = libWindowMassMatch( spec , all_lib_keys)
                if count % 100 == 0:
                    time = timer()
                    print(str(count)+','+str(time-prevtime)+','+str(len(spec['m/z array']))+','+str(spec['precursorMz'][0]['precursorMz'])+','+str(len(lib_keys)))
                    prevtime = time
#                    print(count)
#                    print(timedelta(seconds=timer()))
#                print(str(count)+','+str(time-prevtime)+','+str(len(spec['m/z array']))+','+str(spec['precursorMz'][0]['precursorMz'])+','+str(len(lib_keys)))
#                prevtime = time

                if len(lib_keys) != 0: cosineList = cosine_similarity_msplit(extractLibSpectra(lib, lib_keys), lib_keys, spec['m/z array'], spec['intensity array'])
                else: continue
                for cos in cosineList:
                    if cos[2] > 3:
                        temp = [
                            expSpectraFile, #fileName
                            spec['num'], #scan#
                            spec['precursorMz'][0]['precursorMz'], #MzEXP
                            spec['precursorMz'][0]['precursorCharge'], #zEXP
                            lib[cos[0]]['FullUniModPeptideName'], #peptide
                            lib[cos[0]]['ProteinName'], #protein
                            lib[cos[0]]['PrecursorMz'], #MzLIB
                            lib[cos[0]]['PrecursorCharge'], #zLIB
                            cos[1],#cosine_between_rows(cos_df,0,1), #cosine
                            #cosine_similarity(cos_df[0:1],cos_df[1:])[0][0], #cosine
                            lib[cos[0]]['transition_group_id'], #name
                            len(spec['m/z array']), ##Peak(Query)
                            len(lib[cos[0]]['Peaks']), ##Peaks(Match)
                            cos[2],#len(cos_df.columns), #shared
                            0,#sum(list(cos_df.iloc[1].values)), #ionCount
                            spec['compensationVoltage'], #compensationVoltage
                            spec['precursorMz'][0]['windowWideness'] #totalWindowWidth
                        ]
                        writer.writerow(temp)
        print(count)
    pass

def querySpectrumFilter(specMz, specIntensity, maxNum, ppmTol):
    neighs, remove, left, right = [0], [], 0, 0
#    print()
#    print('original length: '+str(len(specMz)))
    original = len(specMz)
    check = False
    for i in range(len(specMz)):
        for j in range(left,i):
            if not approx(specMz[i], specMz[j], tol=ppmTol):
#                print(neighs)
                neighs.remove(j)
                left = j+1

        if i < len(specMz)-1:
            for j in range(right+1,len(specMz)):
                if not approx(specMz[j], specMz[i], tol=ppmTol):
                    right = j-1
                    break
                neighs.append(j)

#        print('left: '+str(left))
#        print('current: '+str(i))
#        print('right: '+str(right))

        count = 0
        for j in [ i for _,i in sorted(zip(specIntensity,neighs),reverse=True) ]:
            if j == i: break
            if count > maxNum: remove.append(j)
            break
            count += 1
    if len(neighs) > maxNum:
        check = True
    tuples = sorted(zip(specMz,specIntensity))
    if len(remove) > 0:
        print('original list: ')
        for x in tuples: print(x)
        print('to be removed: ')
        for j in remove: print(tuples[j])
    for j in sorted(remove, reverse=True): del tuples[j]
#    print('final length: '+str(len(tuples)))
#    print()
    if check:
        print("here")
    return map(list,zip(*tuples))

def querySpectrumMerge():
    pass

def libSpectrumFilter():
    pass
