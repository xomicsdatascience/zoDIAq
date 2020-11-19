import pandas as pd
from pyteomics import mzxml
from bisect import bisect
import pickle
from sklearn.metrics.pairwise import cosine_similarity
from timeit import default_timer as timer
from datetime import timedelta

def dot(A,B):
    return (sum(a*b for a,b in zip(A,B)))

def cosine_similarity_manual(a,b):
    return dot(a,b) / ( (dot(a,a)**.5) * (dot(b,b)**.5) )

def cosine_between_rows(fullSpec,i1,i2):
    return cosine_similarity_manual(list(fullSpec.iloc[i1].values),list(fullSpec.iloc[i2].values))

def approx(x, y, tol=0.005):
    return abs(x-y) <= tol

def cosine_similarity_msplit(lib_spect, exp_mz, exp_i):
    lib_mz = lib_spect['m/z array']
    lib_i = lib_spect['intensity array']
    lib_mag = lib_spect['magnitude']
    lib_mag = 0.0
    count = 0
    product = 0.0
    exp_mag = 0.0
    i = 0
    j = 0
    while i < len(lib_mz) and j < len(exp_mz):
        if not approx(lib_mz[i],exp_mz[j]):
            if lib_mz[i] > exp_mz[j]: j += 1; continue
            if lib_mz[i] < exp_mz[j]: i += 1; continue
        p = i + 0
        while (p < len(lib_mz)):
            product += lib_i[p]*exp_i[j]
            exp_mag += exp_i[j]**2
            p += 1
            if p==len(lib_mz) or not approx(lib_mz[p], exp_mz[j]): break
        j += 1
    magnitude = lib_mag * (exp_mag**0.5)
    return (product / magnitude if magnitude else 0)

def cosine_similarity_df_generator(lib_mz, lib_i, exp_mz, exp_i):
    df = pd.DataFrame()
    count, i, j = 0, 0, 0
    while i < len(lib_mz) and j < len(exp_mz):
        if not approx(lib_mz[i],exp_mz[j]):
            if lib_mz[i] > exp_mz[j]: j += 1; continue
            if lib_mz[i] < exp_mz[j]: i += 1; continue
        p = i + 0
        while (p < len(lib_mz)):
            df[str(lib_mz[p])+'_'+str(exp_mz[j])] = [ lib_i[p], exp_i[j] ]
            p += 1
            if p==len(lib_mz) or not approx(lib_mz[p], exp_mz[j]): break
        j += 1
    return df

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
    print(timedelta(seconds=timer()))
    lib = {}
    added = set()
    previous_id = 0
    for index, row in lib_df.iterrows():
        id = (row['PrecursorMz'], row['FullUniModPeptideName'])
        if id not in added:
            temp = { 'precursorMz':float(row['PrecursorMz']), 'seq':row['FullUniModPeptideName'], 'precursorCharge':int(row['PrecursorCharge']), 'm/z array':[ float(row['ProductMz']) ], 'intensity array':[ float(row['LibraryIntensity']) ], 'decoy':bool(row['decoy']), 'name':row['transition_group_id'], 'protein':row['ProteinName'] }
            lib[ id ] = temp
            added.add(id)
        else:
            lib[id]['m/z array'].append(float(row['ProductMz']))
            lib[id]['intensity array'].append(float(row['LibraryIntensity']))

    # you can probably do this during the initialization above, but it'll be slightly more complicated. Do this for now.
    for key, value in iter(lib.items()):
        lib[key]['m/z array'], lib[key]['intensity array'] = (list(t) for t in zip(*sorted(zip(lib[key]['m/z array'], lib[key]['intensity array']))))
        lib[key]['magnitude'] = dot(value['intensity array'], value['intensity array'])**.5
    return lib

def cosineDataFrame(content):
    df = pd.DataFrame(content, columns = [
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
    ])
    return df

def expSpectraAnalysis( expSpectraFile, lib ):
    count = 0
    all_lib_keys = sorted(lib)
    final_df_content = []

    with mzxml.read(expSpectraFile) as spectra:
        for spec in spectra:
            count += 1
            if count % 1000 == 0:
                print(count)
                print(timedelta(seconds=timer()))
            lib_keys = libWindowMassMatch( spec , all_lib_keys)
            for x in lib_keys:
                cos_df = cosine_similarity_df_generator(lib[x]['m/z array'], lib[x]['intensity array'], spec['m/z array'][:], spec['intensity array'][:])
                if(len(cos_df)>9):
                    temp = [
                        expSpectraFile, #fileName
                        spec['num'], #scan#
                        spec['precursorMz'][0]['precursorMz'], #MzEXP
                        spec['precursorMz'][0]['precursorCharge'], #zEXP
                        lib[x]['seq'], #peptide
                        lib[x]['protein'], #protein
                        lib[x]['precursorMz'], #MzLIB
                        lib[x]['precursorCharge'], #zLIB
                        cosine_between_rows(cos_df,0,1),
                        #cosine_similarity(cos_df[0:1],cos_df[1:])[0][0], #cosine
                        lib[x]['name'], #name
                        len(spec['m/z array']), ##Peak(Query)
                        len(lib[x]['m/z array']), ##Peaks(Match)
                        len(cos_df.columns), #shared
                        sum(list(cos_df.iloc[1].values)), #ionCount
                        spec['compensationVoltage'], #compensationVoltage
                        spec['precursorMz'][0]['windowWideness'] #totalWindowWidth
                    ]
                    final_df_content.append(temp)
    print(count)
    return(cosineDataFrame(final_df_content))
