import pandas as pd
from pyteomics import mzxml
from bisect import bisect
import pickle
from sklearn.metrics.pairwise import cosine_similarity
pd.set_option("display.max_rows", None, "display.max_columns", None)


'''
# full list
temp_df = pd.read_csv('Part 2/consensus_transitions_lowres_decoys.tsv', sep='\t')

# for sequences
lib_df = temp_df.loc[temp_df['FullUniModPeptideName']=='APIIAVTR']

# for ranges
#lib_df = temp_df.query( 'FullUniModPeptideName == 418 and PrecursorMz < 423' )

# write to file
lib_df.to_csv('Part 2/condensed_APIIAVTR.csv')
'''

def dot(A,B):
    return (sum(a*b for a,b in zip(A,B)))

def cosine_similarity_manual(a,b):
#    print('DF:-----------------------')
#    print('AxB sum: '+str(dot(a,b)))
#    print('A magnitude '+str(dot(a,a)**.5))
#    print('B magnitude: '+str(dot(b,b)**.5))
    return dot(a,b) / ( (dot(a,a)**.5) * (dot(b,b)**.5) )

def cosine_between_rows(fullSpec,i1,i2):
    return cosine_similarity_manual(list(fullSpec.iloc[i1].values),list(fullSpec.iloc[i2].values))

# Tolerance should be a user-determined value
#   They should enter a ppm, which you then multiply by 1e-6
#   This may be the second "mode" referred to in the MSPLIT code
#   "new" tolerance = x*ppm, x = experimental mass
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
#            lib_mag += lib_i[p]**2
            exp_mag += exp_i[j]**2
            p += 1
            if p==len(lib_mz) or not approx(lib_mz[p], exp_mz[j]): break
        j += 1
#    print('Function:-------------------')
#    print('AxB sum: '+str(product))
#    print('A magnitude: '+str(lib_mag))
#    print('B magnitude: '+str(exp_mag**0.5))
    magnitude = lib_mag * (exp_mag**0.5)
#    magnitude = (lib_mag**0.5) * (exp_mag**0.5)
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
    lib = {}
    added = set()
    previous_id = 0
    for index, row in lib_df.iterrows():
        #id = row['transition_group_id']
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

# check to see if the decoy column in the traml corresponds to actual decoys (see it in it's name)
# print how long it takes to reach certain points in the program, esp. 1) loading the full traml df, 2) converting it to a dictionary, and 3) the spectra comparison
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
        'ionCount', #??I presume 'totIonCurrent' from experimental spectra file
        'CompensationVoltage', #the compensation voltage of the experimental SWATH spectrum
        'totalWindowWidth' # width of m/z that was captured in the experimental SWATH spectrum. Corresponds to MzSWATH
    ])
    return df

# Main
lib = tramlFileConversionCSV('Part 2/condensed_APIIAVTR.csv')
all_lib_keys = sorted(lib)

'''
with open('spec_2420.pkl','rb') as f:
    spec = pickle.load(f)
lib_keys = libWindowMassMatch( spec , all_lib_keys)
x = lib_keys[0]
cos = cosine_similarity_msplit(lib[x], spec['m/z array'][:], spec['intensity array'][:])
df2 = cosine_similarity_df_generator(lib[x]['m/z array'], lib[x]['intensity array'], spec['m/z array'][:], spec['intensity array'][:])

cos2 = cosine_between_rows(df2,0,1)
cos2 = cosine_similarity(df2[0:1],df2[1:])

prod = 0
print(cos)
print(cos2[0][0])
print(cos2)
'''
exp_spectra_file = "Part 2/20190411_DI2A_1to16_n1b.mzXML"

def expSpectraAnalysis( expSpectraFile, lib ):
    all_lib_keys = sorted(lib)

    final_df_content = []

    with mzxml.read(expSpectraFile) as spectra:
        for spec in spectra:
            lib_keys = libWindowMassMatch( spec , all_lib_keys)
            for x in lib_keys:
                cos_df = cosine_similarity_df_generator(lib[x]['m/z array'], lib[x]['intensity array'], spec['m/z array'][:], spec['intensity array'][:])
#               print(spec['num'])
#               cos = cosine_between_rows(cos_df,0,1)
#               print(cos)
#               print(cos_df)
#               print(" ")
                temp = [
                    exp_spectra_file, #fileName
                    spec['num'], #scan#
                    spec['precursorMz'][0]['precursorMz'], #MzEXP
                    spec['precursorMz'][0]['precursorCharge'], #zEXP
                    lib[x]['seq'], #peptide
                    lib[x]['protein'], #protein
                    lib[x]['precursorMz'], #MzLIB
                    lib[x]['precursorCharge'], #zLIB
                    cosine_similarity(cos_df[0:1],cos_df[1:])[0][0], #cosine
                    lib[x]['name'], #name
                    len(spec['m/z array']), ##Peak(Query)
                    len(lib[x]['m/z array']), ##Peaks(Match)
                    len(cos_df.columns), #shared
                    sum(list(cos_df.iloc[1].values)), #ionCount
                    spec['compensationVoltage'], #compensationVoltage
                    spec['precursorMz'][0]['windowWideness'] #totalWindowWidth
                ]
                final_df_content.append(temp)
    return(cosineDataFrame(final_df_content))

final_df = expSpectraAnalysis( exp_spectra_file, lib )
final_df.to_csv( 'Part 2/output.csv' )
