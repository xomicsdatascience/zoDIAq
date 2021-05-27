import pandas as pd
from pyteomics import mzxml, mgf, mass
from bisect import bisect, bisect_left
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as pyplot
import numpy as np
from Bio import SeqIO
import idpicker as idp
import re
from collections import defaultdict
from numba import njit
import PooledSpectraMatcher
import QuantSpectraMatcher


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# Template for function descriptions
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''
Function:
Purpose:
Parameters:
Returns:
'''

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# PSM Scoring: Identification
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''
Function: library_file_to_dict()
Purpose: This function reads the file type of the input file (accepts .mgf, .csv and .tsv formats) and runs it through the
            appropriate reading function. Basically, this function is a catch-all used in the corresponding menu function to
            determine how to best read in a particular library spectra file.
Parameters:
    'inFile' - string representing the file path to the desired library spectra file.
Returns:
    dictionary 'lib'
        key - (float, string) tuple.
            float - corresponding to the precursor m/z value of the library spectrum.
            string - corresponding to the full UniMode peptide name/sequence represented by the library
                    spectrum.
        value - dictionary
            'PrecursorCharge': int indicating the charge of the library spectrum.
            'transition_group_id': string indicating library spectrum identifier.
            'ProteinName': string indicating the name of the protein this peptide corresponds to
            'Peaks': list of (float, float, tuple) tuples.
                float 1 - representing the m/z of the peak.
                float 2 - representing the intensity of the peak.
                tuple - equal to the corresponding key in the 'lib' dictionary. This is to identify the
                    library this peak belongs to, as they will be mixed with peaks from other libraries.
'''
def library_file_to_dict(inFile):
    fileType = inFile.split('.')[-1]
    if fileType == 'mgf':
        lib = mgf_library_upload(inFile)
    else:
        lib = traml_library_upload(inFile)
    return lib


'''
Function: traml_library_upload()
Purpose: Given a traml file in .tsv or .csv format representing the library spectra, this function generates a dictionary to
            be used for future analysis. key:value pairings are described in the "Returns" section of this comment. Note that
            the tuple key format was chosen because library spectra will need to be sorted and chosen by precursor mass, so
            identifiers already included in the library spectra were explicitly not used as keys.
Parameters:
    'fileName' - string representing the path to the library spectra file (required traml .tsv or .csv format).
Returns:
    dictionary 'lib' - See dictionary 'lib' key explanation in function library_file_to_dict().
'''
def traml_library_upload(fileName):
    # Library spectra file is read as a pandas dataframe - conditional statement allows for both .tsv and .csv files to be uploaded.
    if fileName.endswith('.tsv'):
        lib_df = pd.read_csv(fileName, sep='\t')
    else:
        lib_df = pd.read_csv(fileName)

    # Print statement for timing the program
    print_milestone('Enter library dictionary upload: ')

    # Pan human and spectraST libraries have different column names. This normalizes the columns.
    headings = traml_column_headings(lib_df.columns)

    # Unneeded columns are removed from the dataframe
    lib_df = lib_df.loc[:, lib_df.columns.intersection([headings['PrecursorMz'],headings['FullUniModPeptideName'],headings['PrecursorCharge'],headings['ProductMz'],headings['LibraryIntensity'],headings['transition_group_id'],headings['ProteinName']])]

    # necessary columns are identified and sorted
    lib_df = lib_df[[headings['PrecursorMz'],headings['FullUniModPeptideName'],headings['PrecursorCharge'],headings['ProductMz'],headings['LibraryIntensity'],headings['transition_group_id'],headings['ProteinName']]]
    lib_df.columns = ['PrecursorMz','FullUniModPeptideName','PrecursorCharge','ProductMz','LibraryIntensity','transition_group_id','ProteinName']

    # Normalize intensities by finding their square root
    lib_df['LibraryIntensity'] = [x**0.5 for x in list(lib_df['LibraryIntensity'])]

    # ID created to become the key of the resulting dictionary
    lib_df['ID'] = list(zip(lib_df['PrecursorMz'].tolist(),lib_df['FullUniModPeptideName'].tolist()))

    # M/z and intensity columns are grouped by ID to be later combined as a list of peaks included in the dictionary
    mz_dict = lib_df.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensity_dict = lib_df.groupby("ID")['LibraryIntensity'].apply(list).to_dict()

    # Dataframe is prepared for and converted to a dictionary
    lib_df.drop_duplicates(subset="ID",inplace=True)
    lib_df = lib_df.loc[:, lib_df.columns.intersection(['ID','PrecursorCharge','transition_group_id','ProteinName'])]
    lib_df.set_index("ID", drop=True, inplace=True)
    lib = lib_df.to_dict(orient="index")

    # pan human library formats are different, including how the peptides are matched to proteins (esp. decoys). This section of code adjusts for this discrepancy.
    if headings['type']=='PanHuman':
        for key, value in lib.items():
            proteins = lib[key]['ProteinName'].split('/')
            num = proteins.pop(0)
            newProteins = [x for x in proteins if 'DECOY' not in x]
            proteinStr = str(len(newProteins))
            for x in newProteins:
                if 'DECOY' in num: proteinStr += ('/DECOY_'+x)
                else: proteinStr += ('/'+x)
            lib[key]['ProteinName'] = proteinStr


    # Peaks list is created and attached to the dictionary
    id = 0
    for key in lib:
        id += 1
        mz, intensity = (list(t) for t in zip(*sorted(zip(mz_dict[key], intensity_dict[key]))))
        keyList = [id for i in range(len(mz))]
        peaks = list(tuple(zip(mz,intensity,keyList)))

        peaks.sort(key=lambda x:x[1],reverse=True)
        if len(peaks) > 10: peaks = peaks[:10]

        peaks.sort(key=lambda x:x[0])
        lib[key]['Peaks'] = peaks
        lib[key]['ID'] = id
        if 'DECOY' in lib[key]['ProteinName']: lib[key]['Decoy'] = 1
        else: lib[key]['Decoy'] = 0
    return lib

'''
Function: traml_column_headings()
Purpose: This function normalizes column headings so the program can accept spectraST and pan human libraries, as they have
            differerent column headings.
Parameters:
    'columns' - a list of strings corresponding to the column names of the inputed library file.
Returns:
    dict - strings with key/value pairings mapping expected values to values found in given library.
'''
def traml_column_headings(columns):
    if 'FullUniModPeptideName' in columns: #SpectraST
        return {
        'type':'SpectraST',
        'PrecursorMz':'PrecursorMz',
        'FullUniModPeptideName':'FullUniModPeptideName',
        'PrecursorCharge':'PrecursorCharge',
        'ProductMz':'ProductMz',
        'LibraryIntensity':'LibraryIntensity',
        'transition_group_id':'transition_group_id',
        'ProteinName':'ProteinName',
        }
    else: #pan human
        return {
        'type':'PanHuman',
        'PrecursorMz':'PrecursorMz',
        'FullUniModPeptideName':'ModifiedPeptideSequence',
        'PrecursorCharge':'PrecursorCharge',
        'ProductMz':'ProductMz',
        'LibraryIntensity':'LibraryIntensity',
        'transition_group_id':'TransitionGroupId',
        'ProteinName':'ProteinId',
        }


'''
Function: mgf_library_upload()
Purpose: Same purpose as traml_library_upload_csv() function, but for mgf files.
Parameters: see traml_library_upload_csv() parameters.
Returns: see traml_library_upload_csv() return values.
'''
def mgf_library_upload(fileName):

    # mgf file is read in using the pyteomics mgf module
    libMGF = mgf.read(fileName)

    # Print statement for timing the program
    print_milestone('Enter library dictionary upload: ')

    # return value is initialized
    lib = {}

    # each spectrum in the mgf file
    id = 0
    for spec in libMGF:
        id += 1

        # key for the final dictionary is initialized
        key = (spec['params']['pepmass'][0], spec['params']['seq'])

        # final dictionary value values are created
        charge = int(re.sub('[+-]','',str(spec['params']['charge'][0])))
        name = spec['params']['title']

        # note that the protein value is optional in MGF files - if it's nonexistent, protein value initialized as an empty string
        if 'protein' in spec['params']: protein = spec['params']['protein']
        else: protein = ''

        if 'DECOY' in name: decoy=1
        else: decoy=0

        # peaks of the library file are intialized.
        mz = spec['m/z array']
        intensity = spec['intensity array']
        intensity = [x**0.5 for x in intensity]
        keyList = [id for x in mz]
        peaks = list(tuple(zip(mz,intensity,keyList)))
        peaks.sort(key=lambda x:x[1],reverse=True)
        if len(peaks) > 10: peaks = peaks[:10]

        peaks.sort(key=lambda x:x[0])

        # final dictionary value is created
        tempDict = {
            'PrecursorCharge':charge,
            'transition_group_id':name,
            'ProteinName':protein,
            'Peaks':peaks,
            'ID':id,
            'Decoy':decoy,

        }

        # entry placed in final dictionary
        lib[key] = tempDict

    return lib


def perform_spectra_pooling_and_analysis(querySpectraFile, outFile, lib, tolerance, maxQuerySpectraToPool, corrected, histFile):

    print_milestone('Begin Grouping Scans by m/z Windows:')
    queWindowDict, queScanValuesDict = group_scans_by_mz_windows(querySpectraFile)

    print('Number of Unpooled MS/MS Query Spectra: ' + str(len(queScanValuesDict)))
    print('Number of Pooled MS/MS Query Spectra/Mz Windows: ' + str(len(queWindowDict)),flush=True)

    # To enhance the print experience, status prints will be given at intervals tailored to the number of identified windows.
    #  example: if there are 1-99 pooled query spectra, print statements are made after every pooled query spectra analysis is complete.
    #           if there are 100-999, print after every 10 pooled spectra. And so on.
    printFriendlyCounter = 100
    while printFriendlyCounter < len(queWindowDict): printFriendlyCounter*=10
    printFriendlyCounter /= 100

    allLibKeys, libIdToKeyDict, libIdToDecoyDict = gather_library_metadata(lib)
    allSpectraMatches = PooledSpectraMatcher.PooledSpectraMatcher()
    numWindowsAnalyzed = 0

    prevtime = timer()
    print_milestone('Begin Pooled Spectra Analysis:')
    with mzxml.read(querySpectraFile, use_index=True) as spectra:

        for precMz_win, scans in queWindowDict.items():
            top_mz = precMz_win[0] + precMz_win[1] / 2
            bottom_mz = precMz_win[0] - precMz_win[1] / 2
            libKeys = identify_lib_spectra_in_window( top_mz, bottom_mz, allLibKeys )
            if len(libKeys) == 0: continue
            pooledLibSpectra = pool_lib_spectra(lib, libKeys)
            pooledQueSpectra = []

            for i in range(len(scans)):
                scanNumber = scans[i]
                queSpectrum = spectra.get_by_id(scanNumber)
                pooledQueSpectra += format_spectra_for_pooling(queSpectrum, scanNumber)

                if (i % maxQuerySpectraToPool == 0 and i!=0) or i == len(scans)-1:
                    pooledQueSpectra.sort()
                    windowSpectraMatches = PooledSpectraMatcher.PooledSpectraMatcher()
                    windowSpectraMatches.compare_spectra(pooledLibSpectra, pooledQueSpectra, tolerance, libIdToDecoyDict)
                    allSpectraMatches.extend_all_spectra(windowSpectraMatches)
                    pooledQueSpectra.clear()

            numWindowsAnalyzed += 1
            if numWindowsAnalyzed % printFriendlyCounter == 0:
                time = timer()
                print('\nNumber of Pooled Experimental Spectra Analyzed: ' + str(numWindowsAnalyzed))
                print('Number of Spectra in Current Pooled Spectra: ' + str(len(scans)))
                print('Time Since Last Checkpoint: ' + str(round(time-prevtime,2)) + ' Seconds', flush=True)
                prevtime = time

    print_milestone('Begin FDR Analysis:')
    maccCutoff = allSpectraMatches.find_score_fdr_cutoff()

    if corrected != -1:
        print_milestone('Begin Correction Process:')
        allSpectraMatches.filter_by_corrected_ppm_window(corrected, maccCutoff, histFile)

        print_milestone('Begin Corrected FDR Analysis:')
        maccCutoff = allSpectraMatches.find_score_fdr_cutoff()

    print_milestone('\nBegin Writing to File: ')
    allSpectraMatches.write_output(outFile, querySpectraFile, maccCutoff, queScanValuesDict, libIdToKeyDict, lib)

def identify_lib_spectra_in_window( top_mz, bottom_mz, sortedLibKeys ):
    temp = sortedLibKeys[:]
    top_key = (top_mz, "z")
    bottom_key = (bottom_mz, "")

    i1 = bisect(temp, bottom_key)
    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
    if i2-i1==1: return []
    temp.insert(i2, top_key)

    return temp[i1+1:i2]

def pool_lib_spectra(lib, libKeys):
    finalList = []
    for key in libKeys: finalList += lib[key]['Peaks']
    return sorted(finalList)


def group_scans_by_mz_windows(querySpectraFile):
    queWindowDict = defaultdict(list)
    queScanValuesDict = defaultdict(dict)

    with mzxml.read(querySpectraFile, use_index=True) as spectra:
        for spec in spectra:

            if 'precursorMz' not in spec: continue
            scan = spec['num']
            precMz = spec['precursorMz'][0]['precursorMz']
            windowWidth = spec['precursorMz'][0]['windowWideness']
            queWindowDict[precMz,windowWidth].append(scan)

            queScanValuesDict[scan]['precursorMz'] = precMz
            queScanValuesDict[scan]['windowWideness'] = windowWidth
            queScanValuesDict[scan]['peaksCount'] = spec['peaksCount']
            if 'compensationVoltage' in spec: CV = spec['compensationVoltage']
            else: CV = ''
            queScanValuesDict[scan]['CV'] = CV

    return queWindowDict, queScanValuesDict

def gather_library_metadata(lib):
    allLibKeys = sorted(lib.keys())
    libIdToKeyDict = {}
    libIdToDecoyDict = {}
    for key in allLibKeys: libIdToKeyDict[lib[key]['ID']] = key; libIdToDecoyDict[lib[key]['ID']] = lib[key]['Decoy']
    return allLibKeys, libIdToKeyDict, libIdToDecoyDict

def format_spectra_for_pooling(spectrum, scanNumber, sqrt=True):
    scanNumber = int(scanNumber)
    if sqrt: intensity = [x**0.5 for x in spectrum['intensity array']]
    else: intensity = [x for x in spectrum['intensity array']]
    peakIDs = [scanNumber for x in range(spectrum['peaksCount'])]
    return list(zip(spectrum['m/z array'],intensity,peakIDs))

def print_milestone(text):
    print(text)
    print(str(timedelta(seconds=timer())),flush=True)

def write_fdr_outputs(inFile, specFile, pepFile, protFile):

    print_milestone('Generating FDR Analysis Files:')
    overallDf = pd.read_csv(inFile).sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
    spectralDf = add_fdr_to_csodiaq_output(overallDf)
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    spectralDf.to_csv(specFile, index=False)
    peptideDf.to_csv(pepFile, index=False)

    if protFile:
        peptideProteinConnections = format_peptide_protein_connections(peptideDf)
        verifiedProteinDict = idp.find_valid_proteins(peptideProteinConnections)
        proteinDf = add_leading_protein_column(peptideDf, verifiedProteinDict)
        tempProtDf = add_fdr_to_csodiaq_output(proteinDf, filterType='leadingProtein')
        proteinMetaInfoDict = tempProtDf.set_index('leadingProtein').T.to_dict()
        proteinDf = remove_invalid_peptides_and_add_metadata(proteinDf, proteinMetaInfoDict)
        proteinDf = mark_peptides_unique_to_proteins(proteinDf)
        proteinDf.to_csv(protFile, index=False)

def format_peptide_protein_connections(peptideDf):
    peptideProteinConnections = []

    for i in range(len(peptideDf)):
        peptide = peptideDf['peptide'].loc[i]
        proteinGroup = peptideDf['protein'].loc[i]

        for protein in proteinGroup.split('/')[1:]:
            peptideProteinConnections.append((peptide,protein))
    return peptideProteinConnections

def remove_invalid_peptides_and_add_metadata(proteinDf, proteinMetaInfoDict):
    removables, proteinCosine, proteinFDR = [], [], []

    for i in range(len(proteinDf)):
        protein = proteinDf['leadingProtein'].loc[i]
        if protein in proteinMetaInfoDict: #dict only contains proteins with FDR < 0.01
            proteinCosine.append(proteinMetaInfoDict[protein]['cosine'])
            proteinFDR.append(proteinMetaInfoDict[protein]['leadingProteinFDR'])
        else:
            removables.append(i)
    proteinDf = proteinDf.drop(proteinDf.index[removables]).reset_index(drop=True)
    proteinDf['proteinCosine'] = proteinCosine
    proteinDf['leadingProteinFDR'] = proteinFDR
    return proteinDf


def mark_peptides_unique_to_proteins(proteinDf):
    proteinDf = proteinDf.sort_values(['proteinCosine','leadingProtein', 'cosine'], ascending=[False,False,False]).reset_index(drop=True)
    uniquePepsDict = defaultdict(set)
    for i in range(len(proteinDf)):
        uniquePepsDict[proteinDf.loc[i]['peptide']].add(proteinDf.loc[i]['leadingProtein'])

    uniquePeps = []
    for i in range(len(proteinDf)):
        p = proteinDf.loc[i]['peptide']
        if len(uniquePepsDict[proteinDf.loc[i]['peptide']]) == 1: uniquePeps.append(1)
        else: uniquePeps.append(0)

    proteinDf['uniquePeptide'] = uniquePeps
    return proteinDf

def add_fdr_to_csodiaq_output(df, filterType='spectral', bestMatchNum=0):

    finalDf = df.copy()

    # final dataframe is initialized as a filtered version of the paraemter input 'df' - filters are minimum number of
    #   allowed matches and highest-scoring unique value of a column type if given
    if filterType != 'spectral': finalDf = finalDf.drop_duplicates(subset=filterType, keep='first').reset_index(drop=True)

    # FDR is calculated
    fdrList, decoyNum = fdr_calculation(finalDf, returnType=1)

    # All rows below the FDR cutoff are removed
    finalDf = finalDf.truncate(after=len(fdrList)-1)

    # FDR column is added to the dataframe
    finalDf[filterType + 'FDR'] = fdrList
    finalDf = finalDf.reset_index(drop=True)
    return finalDf

def add_leading_protein_column(df, verifiedProteinDict):
    # final dataframe is initialized - empty, but with columns matching the 'df' parameter input
    finalDf = pd.DataFrame(columns = df.columns)

    # This list will become the 'leadingProtein' column of the output
    leadingProteins = []

    # for each peptide in the 'df' parameter input
    for i in range(len(df)):

        # proteins connected to this peptide are separated into a list
        proteinGroup = df['protein'].loc[i]
        #proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
        proteins = proteinGroup.split('/')[1:]
        # for each protein connected to this peptide
        for pro in proteins:

            # this ensures that decoys are tagged while non-decoys are not
            #protein = pro[0] + pro[1]
            protein = pro

            # if the protein is one of the proteins verified by the IDPicker algorithm, a new row is added to the input
            #   with the (IDPicker-determined) protein group added in the new column 'leadingProtein'.
            if protein in verifiedProteinDict:
                finalDf = finalDf.append(df.loc[i])
                leadingProteins.append(verifiedProteinDict[protein])

    # new leading protein column set
    finalDf['leadingProtein'] = leadingProteins

    finalDf = finalDf.reset_index(drop=True)

    # For proteins that were part of an IDPicker-determined protein group the same row was added for every protein in the group. This is dropping those duplicate rows.
    finalDf = finalDf.drop_duplicates(keep='first').reset_index(drop=True)

    return finalDf

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# DISPA: Targeted Re-Analysis
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

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
        if (returnType==0 or returnType==1): decoy = 'DECOY' in row['protein']
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

def find_offset_tol(data, histFile, stdev, mean=True):
    if len(data)==0: return 0, 10
    hist, bins = np.histogram(data, bins=200)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    # offset is calculated as the mean or median value of the provided data, though this is overwritten if no corrected standard deviation is provided
    if mean: offset = sum(data)/len(data)
    else: offset = data[len(data)//2]

    # If a corrected standard deviation is provided, it is used to set the tolerance. If not, see below conditional statement
    tolerance = np.std(data)*stdev

    # If no corrected standard deviation is provided, one is customized.
    #  Offset is considered the highest point of the histogram,
    #  tolerance the range necessary before the peaks are the same size as the noise on either end.
    if not stdev:
        index_max = max(range(len(hist)), key=hist.__getitem__)
        histList = list(hist)
        noise = np.mean(hist[:10] + hist[-10:])
        min_i = 0
        max_i = len(hist)-1
        for i in range(index_max, 0, -1):
            if hist[i] < noise: min_i = i; break
        for i in range(index_max, len(hist)):
            if hist[i] < noise: max_i = i; break

        offset = center[index_max]
        if index_max - min_i >= max_i - index_max: tolerance = offset - center[min_i]
        else: tolerance = center[max_i] - offset

    # if a histogram file is provided, it is created, with offset (black) and tolerance (red) lines drawn for reference
    if histFile:
        pyplot.clf()
        pyplot.bar(center, hist, align='center', width=width)
        pyplot.axvline(x=offset, color='black', linestyle = 'dashed', linewidth=4)
        pyplot.axvline(x=offset-tolerance, color='red', linestyle = 'dashed', linewidth=4)
        pyplot.axvline(x=offset+tolerance, color='red', linestyle = 'dashed', linewidth=4)
        pyplot.suptitle('offset: '+str(offset) + ', tolerance: '+str(tolerance))
        pyplot.savefig(histFile)
    return offset, tolerance


'''
Function: calc_heavy_mz()
Purpose: In order to quantify peptides or proteins, we need to re-run DISPA in a targetted fashion on peptides of interest
            identified in previous experiments. Furthermore, we need to target the "heavy" version of the peptide, namely the
            version of the peptide that have heavy lysine and arginine instead of their typical, lighter counterparts. This
            function calculates the expected precursor m/z value for the heavy peptide.
Parameters:
    'seq' - string representing the amino acid composition of the peptide.
    'mz' - float representing the precursor mz for the peptide.
    'z' - int representing the charge of the peptide.
Returns:
    'heavyMz' - float representing the heavy calculation of the precursor mz for the peptide.
'''
def calc_heavy_mz(seq, mz, z):
    hK=8.014199 ## mass of heavy lysine
    hR=10.00827 ## mass of heavy arg

    nK = seq.count('K')
    nR = seq.count('R')
    heavyMz = mz + (nK*hK)/z + (nR*hR)/z
    return heavyMz

def filter_fdr_output_for_targeted_reanalysis(fdrFile, proteins, heavy):
    print_milestone('Generate DISPA Targeted Reanalysis Files:')
    fdrDf = pd.read_csv(fdrFile)
    fdrDf = fdrDf[~fdrDf['protein'].str.contains('DECOY',na=False)].reset_index(drop=True)
    if heavy: fdrDf = fdrDf[fdrDf['peptide'].str.endswith('R') | fdrDf['peptide'].str.endswith('K')].reset_index(drop=True)
    if proteins:
        fdrDf = fdrDf[fdrDf['uniquePeptide']==1].sort_values('ionCount', ascending=False).reset_index(drop=True)
        fdrDf = fdrDf.groupby(['leadingProtein']).head(proteins).reset_index()
    return fdrDf

def gather_all_possible_cv_values(fdrDf):
    CVs = set(fdrDf['CompensationVoltage'])
    def notNan(x): return ~np.isnan(x) # get list of all CVs, or CV==0 if there are no CVs
    CVs = set(filter(notNan, CVs))
    if len(CVs)==0: CVs.add('')
    return CVs

def calculate_all_light_heavy_mzs(cvDf):
    lightMzs = []
    heavyMzs = []
    peptides = []
    for index, row in cvDf.iterrows():
        peptide = row['peptide']
        lightMz = float(row['MzLIB'])
        charge = row['zLIB']
        heavyMz = calc_heavy_mz(peptide, lightMz, charge)
        lightMzs.append(lightMz)
        heavyMzs.append(heavyMz)
        peptides.append(peptide)
    return lightMzs, heavyMzs, peptides

def compile_reanalysis_files_with_bins(cvDf, heavy):
    lightMzs, heavyMzs, peptides = calculate_all_light_heavy_mzs(cvDf)
    binWidth = 1.0
    lv, lBins = bin_assignment(lightMzs, binWidth)
    hv, hBins = bin_assignment(heavyMzs, binWidth)
    scanLightMzs = [lBins[lv[i]] for i in range(len(cvDf))]
    scanHeavyMzs = [hBins[hv[i]] for i in range(len(cvDf))]

    cvDf['scanLightMzs'] = [x-(binWidth/2) for x in scanLightMzs]
    cvDf['scanHeavyMzs'] = [x-(binWidth/2) for x in scanHeavyMzs]

    # preparing the allCSV output
    binDict = defaultdict(list)
    for i in range(len(cvDf)):
        if heavy: binDict[scanLightMzs[i],scanHeavyMzs[i]].append(peptides[i])
        else: binDict[scanLightMzs[i],0].append(peptides[i])

    # preparing the files that are fed into an MS machine
    data = []
    binKeys = sorted(binDict)
    for i in range(len(binKeys)):
        data.append(['/'.join(binDict[binKeys[i]]), '', '(no adduct)', binKeys[i][0]-(binWidth/2), 2, i+1])
        if heavy: data.append(['/'.join(binDict[binKeys[i]]), '', '(no adduct)', binKeys[i][1]-(binWidth/2), 2, i+1])
    return cvDf, data

def compile_reanalysis_files_without_bins(cvDf, heavy):
    data = []
    scanLightMzs = []
    scanHeavyMzs = []
    for i in range(len(cvDf)):
        compound = cvDf.loc[i]['peptide']
        formula = ''
        adduct = '(no adduct)'
        lightMz = float(cvDf.loc[i]['MzLIB'])
        charge = cvDf.loc[i]['zLIB']
        heavyMz = calc_heavy_mz(compound, lightMz, charge)
        MSXID = i+1
        scanLightMzs.append(round(lightMz, ndigits = 2))
        scanHeavyMzs.append(round(heavyMz, ndigits = 2))
        data.append([compound, formula, adduct, scanLightMzs[-1], charge, MSXID])
        if heavy: data.append([compound, formula, adduct, scanHeavyMzs[-1], charge, MSXID])

    cvDf['scanLightMzs'] = scanLightMzs
    cvDf['scanHeavyMzs'] = scanHeavyMzs
    return cvDf, data


def write_targeted_reanalysis_outputs(header, fdrDf, heavy):
    bins = True
    CVs = gather_all_possible_cv_values(fdrDf)
    allDfs = []
    for CV in CVs:
        if CV: cvDf = fdrDf[fdrDf['CompensationVoltage']==CV].reset_index(drop=True)
        else: cvDf = fdrDf
        if bins: cvDf, data = compile_reanalysis_files_with_bins(cvDf, heavy)
        else: cvDf, data = compile_reanalysis_files_without_bins(cvDf, heavy)
        allDfs.append(cvDf)
        finalDf = pd.DataFrame(data, columns=['Compound','Formula','Adduct','m.z','z','MSXID'])
        outFile = header
        if CV: outFile += '_' + str(CV)
        outFile += '.txt'
        finalDf.to_csv(outFile, sep='\t', index=False)

    #final print
    compDf = pd.concat(allDfs)
    if bins: binMark = '_withBins_'
    else: binMark = '_withoutBins_'
    compDf.to_csv(header+binMark+'allCVs.csv', index=False)




def bin_assignment(mzValues, binWidths):
    maxi = max(mzValues)+binWidths
    mini = min(mzValues)
    bins = np.arange(mini, maxi, binWidths)
    bins = [round(x,2) for x in bins]
    values = np.digitize(mzValues, bins)
    return values, bins


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# Light:Heavy Quantification
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
@njit
def approx(x, y, ppmTol):
    if x==y: return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)


'''
Function: ppm_offset()
Purpose: Given a peak mz value and a ppm offset value, this function calculates the mz offset that would need to be
            subtracted from the mz value to align the final spread of ppm differences around 0.
Parameters:
    'mz' - Mz value of a given peak (generally a query peak)
    'ppm' - Ppm offset value. This is generally calculated from an uncorrected csodiaq output. During the uncorrected
        csodiaq output generation, this value is generally 0, which in turn makes the returned mz offset 0. NOTE: the
        ppm offset is made negative to apply appropriately to the dataset.
Returns:
    float - Mz offset value.
'''
@njit
def ppm_offset(mz, ppm):
    return (mz*-ppm)/1000000


'''
Function: pool_library_spectra_by_scan()
Purpose: This function reads the file type of the input file (accepts .mgf, .csv and .tsv formats) and runs it through the
            appropriate reading function. Basically, this function is a catch-all used in the corresponding menu function to
            determine how to best read in a particular library spectra file. This is similar to the library_file_to_dict()
            function, but is specific to quantification library prep.
Parameters:
    'libFile' - string representing the path to the library spectra, ideally the same library that was used in steps leading
        up to FDR analysis.
    dictionary 'libScanDict'
        key - (float, string) tuple.
            float - corresponding to the precursor m/z value of the library spectrum.
            string - corresponding to the full UniMode peptide name/sequence represented by the library
                    spectrum.
        value - string representing the scan number of the DISPA re-analysis. Scan numbers are consistent across datasets
            from the same run.
    'fragDf' - Pandas dataframe, created from 'inFile' parameter of connect_mzxml_csodiaq_library() function.
    'maxPeaks' - maximum number of peaks to be included in each library (light and heavy combined will have at most double
        this number)
Returns:
    dictionary - see parameter 'libDict' in heavy_light_quantification() function.
'''
def pool_library_spectra_by_scan(libFile, libScanDict, fragDf, maxPeaks):
    fileType = libFile.split('.')[-1]
    if fileType == 'mgf':
        uniquePeps = set(fragDf['peptide'])
        digPat = r'\+\d+\.\d+'
        uniqueDigs = set()
        for pep in uniquePeps: uniqueDigs.update(re.findall(digPat, pep))
        digDict = dict(zip(uniqueDigs, [str(x) for x in list(range(len(uniqueDigs)))]))
        customAAmass = dict(mass.std_aa_mass)
        for key in digDict: customAAmass[digDict[key]] = float(key[1:])
        return mgf_library_upload_quant(libFile, libScanDict, digDict, customAAmass, maxPeaks)
    else:
        return traml_library_upload_quant(libFile, libScanDict, maxPeaks)


'''
Function: mgf_library_upload_quant()
Purpose: This function creates a dictionary with scan:peak entries for identifying library spectra based on scan number. This
            function is specific to MGF library files, but is specific to quantification, distinguishing it from the
            mgf_library_upload() function.
Parameters:
    'fileName' - string representing the path to the library spectra, ideally the same library that was used in steps leading
        up to FDR analysis.
    dictionary 'scanDict' - see parameter 'libScanDict' in pool_library_spectra_by_scan() function.
    dictionary 'digDict'
        key - numeric, single-digit string representing a modification found in a peptide sequence.
        value - string representing the mass of the modification (float).
    dictionary 'aaDict' - same as mass.std_aa_mass, but with values from digDict added.
    'maxPeaks' - maximum number of peaks to be included in each library (light and heavy combined will have at most double
        this number)
Returns:
    dictionary 'lib' - see parameter 'libDict' in heavy_light_quantification() function.
'''
def mgf_library_upload_quant(fileName, scanDict, digDict, aaDict, maxPeaks):

    # mgf file is read in using the pyteomics mgf module
    libMGF = mgf.read(fileName)

    # return value is initialized
    lib = defaultdict(list)

    keyList = sorted(list(scanDict.keys()))
    # each spectrum in the mgf file
    for spec in libMGF:

        seq = spec['params']['seq']
        precMz = spec['params']['pepmass'][0]

        key = (round(precMz,2), seq)
        if key not in scanDict: continue

        # Decimal values are replaced with numeric placeholders to be included in the analysis.
        sequence = re.sub(r'\+\d+\.\d+', lambda m: digDict.get(m.group()), seq)

        # peaks of the library file are intialized
        mz = list(spec['m/z array'])
        intensity = [x for x in list(spec['intensity array'])]
        z = spec['params']['charge'][0]

        # The y-ion mz value for each fragment of the peptide is calculated. If it is in the library, it and it's intensity are stored in a list
        # NOTE: y-ions are singled out because they should have at least one lysine or arginine, so will have a heavy counterpart that can show up. B-ions don't have that guarantee.
        fragList = []
        for x in range(1, len(sequence)-1):
            fragseq = sequence[x:]
            lightfragmz = mass.fast_mass(sequence=sequence[x:], ion_type='y', charge=1, aa_mass = aaDict) # Do I need to use different possible charges?
            i = approx_list(lightfragmz, mz)
            if i==-1: continue
            fragList.append((intensity[i], lightfragmz, fragseq))

        # y-ion peaks are sorted by intensity, and lower-intensity peaks are filtered out.
        fragList.sort(reverse=True)
        if maxPeaks !=0 and len(fragList) >= maxPeaks: fragList = fragList[:maxPeaks]

        # heavy counterpart mz is calculated. Light and heavy pairs are additionally tagged by their intensity rank and included in the final output.
        peaks = []
        for i in range(len(fragList)):
            fragMz = fragList[i][1]
            fragInt = fragList[i][0]
            peaks.append((fragMz, fragInt, (0, i, seq)))
            peaks.append((calc_heavy_mz(fragList[i][2], fragMz, 1), fragInt, (1, i, seq)))

        # peaks are sorted by mz value (as they appear in a spectrum)
        peaks.sort(key=lambda x:x[0])

        # @DEBUG: error 1
        # entry placed in final dictionary
        #lib[scanDict[keyMatch]] = peaks
        lib[scanDict[key]] += peaks
    return lib


'''
Function: traml_library_upload_quant()
Purpose: This function creates a dictionary with scan:peak entries for identifying library spectra based on scan number. This
            function is specific to TraML library files, but is specific to quantification, distinguishing it from the
            traml_library_upload() function.
Parameters:
    'fileName' - string representing the path to the library spectra, ideally the same library that was used in steps leading
        up to FDR analysis.
    dictionary 'scanDict' - see parameter 'libScanDict' in pool_library_spectra_by_scan() function.
    'maxPeaks' - maximum number of peaks to be included in each library (light and heavy combined will have at most double
        this number)
Returns:
'''
def traml_library_upload_quant(fileName, scanDict, maxPeaks):

    # Library spectra file is read as a pandas dataframe - conditional statement allows for both .tsv and .csv files to be uploaded.
    if fileName.endswith('.tsv'):
        lib_df = pd.read_csv(fileName, sep='\t')
    else:
        lib_df = pd.read_csv(fileName)

    # to save on time, an initial filtering step removes all non-relevant peptides from consideration
    peptides = set([x[1] for x in scanDict])
    lib_df = lib_df[lib_df['FullUniModPeptideName'].isin(peptides)]

    # Unneeded columns are removed from the dataframe
    lib_df = lib_df.loc[:, lib_df.columns.intersection(['PrecursorMz','PeptideSequence','FullUniModPeptideName','ProductMz','LibraryIntensity','FragmentCharge','FragmentType','FragmentSeriesNumber'])]

    # Rounding to two decimal places to match the re-analysis files
    lib_df['PrecursorMz'] = [round(x,2) for x in lib_df['PrecursorMz']]

    # ID created to become the key of the resulting dictionary
    lib_df['ID'] = list(zip(lib_df['PrecursorMz'].tolist(),lib_df['FullUniModPeptideName'].tolist()))

    # Dataframe filters out rows with an ID that is not found in the scanDict dictionary
    lib_df = lib_df[lib_df['ID'].isin(scanDict)]

    # b-ions are removed from consideration.
    # NOTE: y-ions are singled out because they should have at least one lysine or arginine, so will have a heavy counterpart that can show up. B-ions don't have that guarantee.
    lib_df = lib_df[lib_df['FragmentType']=='y']

    # M/z and intensity columns are grouped by ID to be later combined as a list of peaks included in the dictionary
    mz_dict = lib_df.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensity_dict = lib_df.groupby("ID")['LibraryIntensity'].apply(list).to_dict()

    # Dataframe is prepared for and converted to a dictionary
    lib_df.drop_duplicates(subset="ID",inplace=True)
    lib_df = lib_df.loc[:, lib_df.columns.intersection(['ID', 'PeptideSequence', 'FragmentCharge', 'FragmentSeriesNumber'])]
    lib_df.set_index("ID", drop=True, inplace=True)
    lib = lib_df.to_dict(orient="index")

    # Peaks list is created and attached to the dictionary
    finalDict = defaultdict(list)
    for key in lib:

        # y-ion peaks are sorted by intensity, and lower-intensity peaks are filtered out.
        fragList = sorted(list(tuple(zip(intensity_dict[key], mz_dict[key]))), reverse=True)
        if maxPeaks !=0 and len(fragList) >= maxPeaks: fragList = fragList[:maxPeaks]

        # heavy counterpart mz is calculated. Light and heavy pairs are additionally tagged by their intensity rank and included in the final output.
        peaks = []
        for i in range(len(fragList)):
            fragMz = fragList[i][1]
            fragInt = fragList[i][0]
            peaks.append((fragMz, fragInt, (0,i,key[1])))
            fragSeq = lib[key]['PeptideSequence'][-lib[key]['FragmentSeriesNumber']:]
            heavyMz = calc_heavy_mz(fragSeq, fragMz, lib[key]['FragmentCharge'])
            peaks.append((heavyMz, fragInt, (1,i,key[1])))

        # entry placed in final dictionary
        finalDict[scanDict[key]] += peaks

    return finalDict


def create_mzxml_to_csodiaq_dict(mzxmlFile):
    mzxmlToCsodiaqDict = {}
    with mzxml.read(mzxmlFile) as spectra:
        for x in spectra:
            key = ( round(x['precursorMz'][0]['precursorMz'], 2),
                    round(x['precursorMz'][1]['precursorMz'], 2),
                    x['compensationVoltage']
                    )
            mzxmlToCsodiaqDict[key] = x['num']
    return mzxmlToCsodiaqDict

def connect_csodiaq_data_to_scans(idFile, mzxmlToCsodiaqDict, fragDf):
    fragVarDict = defaultdict(list)
    libScanDict = {}


    for i in range(len(fragDf)):

        # @DEBUG: Top line for my ouput, bottom for Jesse's old output.
        seq, mz, z, CV = fragDf.loc[i]['peptide'], fragDf.loc[i]['MzLIB'], fragDf.loc[i]['zLIB'], fragDf.loc[i]['CompensationVoltage']
        #seq, mz, z, CV = fragDf.loc[i]['Peptide'], fragDf.loc[i]['prec_light_mz'], fragDf.loc[i]['z'], fragDf.loc[i]['CV']
        lightMz = round(mz, 2)

        key = (round(fragDf.loc[i]['scanLightMzs'],2), round(fragDf.loc[i]['scanHeavyMzs'],2), CV)
        #key = (lightMz, heavyMz, CV)
        if key in mzxmlToCsodiaqDict:
            scan = mzxmlToCsodiaqDict[key]
            libScanDict[lightMz, seq] = scan
            fragVarDict[scan].append({'seq':seq, 'mz':mz, 'z':z, 'CV':CV})
    return fragVarDict, libScanDict

def connect_mzxml_to_csodiaq_and_library(idFile, libFile, mzxmlFiles, maxPeaks):
    print_milestone('Preparing Quantification Dictionaries:')
    metadataToScanDict = create_mzxml_to_csodiaq_dict(mzxmlFiles[0])
    fileType = idFile.split('.')[-1]
    if fileType == 'csv': fragDf = pd.read_csv(idFile)
    else: fragDf = pd.read_csv(idFile, sep='\t')

    scanToCsodiaqDict, libMetadataToScanDict = connect_csodiaq_data_to_scans(idFile, metadataToScanDict, fragDf)
    scanToLibPeaksDict = pool_library_spectra_by_scan(libFile,libMetadataToScanDict, fragDf, maxPeaks)
    return scanToCsodiaqDict, scanToLibPeaksDict

@njit()
def spectra_peak_comparison_quant(libMzs, queMzs, ppmTol, ppmOffset):

    # initializing various values to be used in the function, including the return values.
    lenLib = len(libMzs)
    lenQue = len(queMzs)
    pMatches = []
    jMatches = []
    ppmMatches = []

    # By tracking the indices of the current library/query peaks we reduce the time complexity of the algorithm
    i, j = 0, 0
    expPeakMz = queMzs[j] - ppm_offset(queMzs[j], ppmOffset)
    while i < lenLib and j < lenQue:

        # If the m/z of the peaks are not within the given ppm tolerance, the indices of the smaller of the two is incremented
        #   and we are returned to the top of the while loop.
        if not approx(libMzs[i],expPeakMz, ppmTol):
            if libMzs[i] > expPeakMz:
                j += 1
                if j < lenQue: expPeakMz = queMzs[j] - ppm_offset(queMzs[j], ppmOffset)
                continue
            if libMzs[i] < expPeakMz: i += 1; continue

        # To account for the matching of one query peak to multiple library peaks, library peaks are looped over
        #   after the initial match. Every matching peak contributes to the various variables in the final returned dictionary.
        p = i + 0
        while (p < lenLib):
            ppm = approx(libMzs[p], expPeakMz, ppmTol)
            if p==lenLib or not ppm: break
            pMatches.append(p)
            jMatches.append(j)
            ppmMatches.append(ppm)
            p += 1

        # Note that the possibility of one library peak matching to multiple query peaks is automatically accounted for
        #   by the fact that the query peak is the next default increment after all match calculations have been made.
        j += 1
        if j < lenQue: expPeakMz = queMzs[j] - ppm_offset(queMzs[j], ppmOffset)
    return pMatches, jMatches, ppmMatches

'''
Function: heavy_light_quantification()
Purpose: Following DISPA targeted re-analysis, this function calculates changes in quantity of the peptides identified in
            previous calculations relative to the given (heavy) control. The heavy control represents a culture grown with
            heavy lysine and arginine, the two amino acids at which trypsin cuts peptides. By mixing the heavy control with
            a light counterpart representing a change in environment in equal quantities, one can determine the impact of the
            environment change on the peptide expession. If the heavy and light cultures were grown in the same conditions,
            it would be expected that the relative ratio for each peptide would be 1:1.
Parameters:
    dictionary 'fragDict'
        key - string representing the scan number of the DISPA re-analysis. Scan numbers are consistent across datasets from
            the same run.
        value - list of dictionaries, dictionary values below.
            'seq': sequence of a peptide corresponding to the scan of interest.
            'mz': precursor mass of the peptide.
            'z': charge of the peptide.
            'CV': Compensation Voltage setting where the peptide was originally found.
    dictionary 'libDict'
        key - string representing the scan number of the DISPA re-analysis. Scan numbers are consistent across datasets from
            the same run. NOTE: this library dictionary ONLY represents peptides of interest.
        value - See dictionary 'lib' value 'Peaks' explanation in function library_file_to_dict(), with one difference.
            Rather than being be tagged by a library-specific tuple, the tags are tuples of (str, int, str) composition. The
            string differentiates between 'light' and 'heavy' peaks, and the int represents the rank of peak intensity (0
            being high), and the last string represents the peptide represented by the library.
    'mzxmlFiles' - list of strings representing the path to files from the DISPA re-analysis.
Returns:
    'finalDf' - pandas dataframe. The first two columns of the dataframe are the scans and peptides, respectively. Each
        subsequent column represents the ratios for those peptides derived from each data file of the DISPA re-analysis.

'''

def initialize_quantification_output(fragDict, libDict):
    finalDf = pd.DataFrame()

    tempKeys = []
    for scan in libDict:
        for tempDict in fragDict[scan]:
            tempKeys.append((scan, tempDict['seq']))

    scans = []
    peptides = []
    for x in sorted(tempKeys):
        scans.append(x[0])
        peptides.append(x[1])

    finalDf['scan'] = scans
    finalDf['peptide'] = peptides
    return finalDf

def heavy_light_quantification(fragDict, libDict, mzxmlFiles, outDir, massTol, minMatch, ratioType, correction, hist):

    finalDf = initialize_quantification_output(fragDict, libDict)
    def initialize_ratio_dict_values(): return np.nan

    # Heavy:Light Ratio for each peptide is calculated for every DISPA reanalysis file. Looping over the files begins here.
    for f in mzxmlFiles:
        # In order to apply a ppm correction, optimal offset and tolerance for the data.
        ppmDiffs = []
        allSpectraMatch = QuantSpectraMatcher.QuantSpectraMatcher()
        scanToNoiseIntensityCutoffDict = dict()
        with mzxml.read(f, use_index =True) as file:
            for scan in sorted(libDict.keys()):

                # experimental spectrum is formatted for peak comparison (see 'libDict' parameter value for peak description)
                spec = file.get_by_id(scan)

                scanToNoiseIntensityCutoffDict[int(scan)] = np.mean(sorted(spec['intensity array'])[:10])/2

                expSpectrum = format_spectra_for_pooling(spec, scan, sqrt=False)
                expSpectrum.sort()

                libSpectra = sorted(libDict[scan])

                quantSpectraMatch = QuantSpectraMatcher.QuantSpectraMatcher()
                quantSpectraMatch.compare_spectra(libSpectra, expSpectrum, massTol, minMatch)
                allSpectraMatch.extend_all_spectra(quantSpectraMatch)

        if correction != -1: allSpectraMatch.filter_by_corrected_ppm_window(correction, hist, minMatch)
        ratioDict = defaultdict(initialize_ratio_dict_values)
        if len(allSpectraMatch.libraryIntensities) != 0: ratioDict = allSpectraMatch.determine_ratios(ratioDict, scanToNoiseIntensityCutoffDict, ratioType, minMatch)


        finalDf[f] = [ratioDict[(int(row['scan']),row['peptide'])] for index, row in finalDf.iterrows()]
    print_milestone('Finish SILAC Quantification')
    return finalDf

def make_lib_dict_quant(libFile, libScanDict, fragDf, maxPeaks):
    fileType = libFile.split('.')[-1]
    if fileType == 'mgf':

        # all code before the return statement assumes that the peptides have decimal values corresponding to non-standard modifications.
        uniquePeps = set(fragDf['peptide'])
        digPat = r'\+\d+\.\d+'
        uniqueDigs = set()
        for pep in uniquePeps: uniqueDigs.update(re.findall(digPat, pep))
        digDict = dict(zip(uniqueDigs, [str(x) for x in list(range(len(uniqueDigs)))]))

        #NOTE: If you have more than 10 custom entries you'll have problems
        customAAmass = dict(mass.std_aa_mass)
        for key in digDict: customAAmass[digDict[key]] = float(key[1:])
        return mgf_library_upload_quant(libFile, libScanDict, digDict, customAAmass, maxPeaks)
    else:
        return traml_library_upload_quant(libFile, libScanDict, maxPeaks)

def initialize_return_values(tag):
    if tag=='qPPM':
        light = defaultdict(list)
        heavy = defaultdict(list)
        return [light, heavy]
    elif tag=='ratio':
        light = defaultdict(dict)
        heavy = defaultdict(dict)
        return [light, heavy]

def update_return_values(returns, peak1, peak2, i1, i2, ppm, tag):
    if tag=='qPPM':
        if peak1[2][0] == 0: returns[0][peak1[2][2]].append(ppm)
        if peak1[2][0] == 1: returns[1][peak1[2][2]].append(ppm)
    elif tag=='ratio':
        if peak1[2][0] == 0: returns[0][peak1[2][2]][peak1[2][1]] = peak2[1]
        if peak1[2][0] == 1: returns[1][peak1[2][2]][peak1[2][1]] = peak2[1]


def match_score(keys):
    if len(keys) == 0: return False
    if min(keys) < 3: return True
    return False


'''
Function: return_ratio()
Purpose: This function may or may not have use in the future. There are often several ratios calculated for any given
            peptide, one for each light:heavy peak pair identified. This function con return a different value depending on
            the type of ratio requested. See comments in the function itself for specifics.
Parameters:
    NOTE: both dictionaries have the same keys, but different values.
    dictionary 'lightDict'
        key - int representing the rank of intensity for an identified light:heavy peak pair.
        value - intensity of the light peak of the pair.
    dictionary 'lightDict'
        key - int representing the rank of intensity for an identified light:heavy peak pair.
        value - intensity of the heavy peak of the pair.
    'type' - string representing the ratio-deciding method of choice.
Returns:
    float representing the ratio calculated. Returns string 'invalidType' if an invalid type is provided in the 'type'
        parameter.
'''
def return_ratio(lightDict, heavyDict, type):
    keys = lightDict.keys()

    # intensity: the ratio from the peak pair with the highest intensity is returned.
    minKey = min(keys)
    if type=='intensity': return np.log2((heavyDict[minKey]/lightDict[minKey]))

    lightInt = [lightDict[key] for key in keys]
    heavyInt = [heavyDict[key] for key in keys]
    log2Ratios = [np.log2(x/y) for x,y in zip(heavyInt, lightInt)]

    # mean: the mean value of log2(ratio) is returned.
    if type=='mean': return np.mean(log2Ratios)

    # median: the median value of log2(ratio) is returned.
    if type=='median': return np.median(log2Ratios)

    # weighted: a weighted mean calculation. Peak pairs with a higher intensity contribute more to the average.
    if type=='weighted':
        ratio = 0.0
        sumInt = sum(lightInt + heavyInt)
        for i in range(len(lightInt)): ratio += ((lightInt[i]+heavyInt[i])/sumInt)*log2Ratios[i]
        return ratio

    # ratio_median: the median value of un-normalized, raw ratios is returned.
    normRatios = [x/y for x,y in zip(heavyInt, lightInt)]
    if type=='ratio_median': return np.median(normRatios)

    # ratio_mean: the mean value of un-normalized, raw ratios is returned.
    if type=='ratio_mean': return np.mean(normRatios)

    # if the type doesn't match one of the above options, an error message is returned instead.
    return 'invalidType'


'''
Function: approx_list()
Purpose: This function determines if a value is within a given ppm tolerance of a value in a list.
Parameters:
    'x' - float value.
    'l' - list of float values.
    'ppmTol' - the ppm tolerance for float values to be considered a match. Default is 10ppm.
Returns:
    int - the index of the match. If there is no match, returns -1 instead.
'''
def approx_list(x, l, ppmTol=10):
    for i in range(len(l)):
        if approx(x, l[i], ppmTol): return i
    return -1



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
                if approx_list(peaks[i][0],mzValues)==-1: peaks.pop(i)
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
