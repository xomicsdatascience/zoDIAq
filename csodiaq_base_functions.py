import pandas as pd
from pyteomics import mzxml, mgf, mass
from bisect import bisect, bisect_left
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as pyplot
import numpy as np
import idpicker as idp
import re
from collections import defaultdict
import linecache
from Bio import SeqIO
from numba import njit, jit, typeof
from numba.core import types
from numba.typed import List

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
    print('\nEnter library dictionary upload: ',flush=True)
    print(str(timedelta(seconds=timer())),flush=True)

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
    if 'FullUniModPeptideName' in columns:
        return {
        'PrecursorMz':'PrecursorMz',
        'FullUniModPeptideName':'FullUniModPeptideName',
        'PrecursorCharge':'PrecursorCharge',
        'ProductMz':'ProductMz',
        'LibraryIntensity':'LibraryIntensity',
        'transition_group_id':'transition_group_id',
        'ProteinName':'ProteinName',
        'type':'SpectraST'
        }
    else:
        return {
        'PrecursorMz':'PrecursorMz',
        'FullUniModPeptideName':'ModifiedPeptideSequence',
        'PrecursorCharge':'PrecursorCharge',
        'ProductMz':'ProductMz',
        'LibraryIntensity':'LibraryIntensity',
        'transition_group_id':'TransitionGroupId',
        'ProteinName':'ProteinId',
        'type':'PanHuman'
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
    print('\nEnter library dictionary upload: ',flush=True)
    print(str(timedelta(seconds=timer())),flush=True)

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
        }

        # entry placed in final dictionary
        lib[key] = tempDict



    return lib


'''
Function: lib_mz_match_query_window()
Purpose: This function determines which library spectra fit within the baseline precursor mass window of the
            experimental spectrum. This significantly reduces the time complexity of the overall analysis by eliminating
            all irrelevant library spectra from consideration, the vast majority of them.
Parameters:
    'spec' - dictionary corresponding to pyteomics.mzxml.read() values. This variable contains all data corresponding
        to the query spectrum in question.
    'sortedLibKeys' - a list of (float, string) tuples. List represents all keys of the dictionary 'lib' from
        library_file_to_dict(). See dictionary 'lib' key explanation in function library_file_to_dict().
Returns:
    A list of (float, string) tuples. Represents the section of values from the 'sortedLibKeys' parameter that
        are relevant to the query spectrum represented by the 'spec' parameter.
'''
def lib_mz_match_query_window( top_mz, bottom_mz, sortedLibKeys ):
    # Values will be added to sortedLibKeys parameter, so a copy of that list is created here.
    temp = sortedLibKeys[:]

    # A fake 'key' value is created representing the upper limit of keys that correspond to the query spectrum.
    top_key = (top_mz, "z")

    # A fake 'key' value is created representing the lower limit of keys that correspond to the query spectrum.
    bottom_key = (bottom_mz, "")

    # The fake keys created above are inserted in an O(log(n)) fashion to the sorted library keys.
    i1 = bisect(temp, bottom_key)

    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
    if i2-i1==1: return []

    temp.insert(i2, top_key)

    # All keys between the upper and lower limit fake 'keys' are returned.
    return temp[i1+1:i2]


'''
Function: pool_lib_spectra()
Purpose: Given a list of library spectra (represented by a list of library keys), a comprehensive library spectrum
            is generated and returned, representing the pooled peaks of several library spectra. This variable
            signifincantly reduces the time complexity of the algorithm.
Parameters:
    'lib' - dictionary as returned by the library_file_to_dict() function.
    'libKeys' - list of keys corresponding to the 'lib' parameter as returned by the lib_mz_match_query_window()
        function.
Returns:
    A list of (float, float, tuple) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict().
'''
def pool_lib_spectra(lib, libKeys):
    finalList = []
    for key in libKeys: finalList += lib[key]['Peaks']
    return sorted(finalList)


'''
Function: approx()
Purpose: Given two float values, this function returns a float value indicating the difference between two peaks m/z.
            Note that ppm is proportional, not static. For example, the ppm calculated for a difference between 1 and 2
            will be much higher than a ppm calculated for a difference between 1,000,000 and 1,000,001, even though
            statically the difference is identical.
Parameters:
    'x' - first float value, corresponding to the theoretical mass (library mz).
    'y' - second float value, corresponding to the experimental mass (query mz).
    'ppmTol' - ppm used to determine tolerance.
Returns:
    float - Returns the difference (in ppm) of the two values if the difference is less than the provided tolerance.
            If the difference is outside the provided tolerance, 0 is returned instead for conditional purposes. If the
            values are exactly equal, a value close to 0 is returned instead (1x10e-7).
'''
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
Function: spectra_peak_comparison()
Purpose: This function determines peaks that are within a sufficiently close tolerance to one another. Various functions
            requiring peak comparison then take place inside this function, including compiling data for a cosine
            similarity score. The specifics of the various functions can be found in the return values below.
Parameters:
    'libMzs' - numpy float array. Represents library spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict().
    'queMzs' - numpy float array. Represents query spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict().
    'ppmTol' - see 'ppmTol' parameter description for function approx().
    'ppmYOffset' - The value in ppm that should be subtracted from the query peak m/z. This is the average ppm difference
        of an uncorrected csodiaq output, being applied to the experimental spectra to generate a corrected
        csodiaq calculation output.
Returns:
    tag 'identify':
        dictionary 'cosDict'
            key - ((float, string),int) tuple. See dictionary 'lib' key explanation in function library_file_to_dict() for
                the first part of the tuple. The second part is the scan number of a query spectrum.
            value - list of floats
                float - sum of products of intensities between matched library and query peaks. In context of the cosine
                    similarity score algorithm, this is the "A*B" sum.
                float - sum of squared query peak intensities when matched to a library peaks. In context of the cosine
                    similarity score algorithm, this is the "B^2" sum.
                float - sum of squared library peak intensities when matched to a query peaks. In context of the cosine
                    similarity score algorithm, this is the "A^2" sum.
        dictionary 'countDict'
            key - ((float, string),int) tuple. See dictionary 'cosDict' key explanation in this function.
            value - int
                Representing the number of matched peaks. This value will be the 'shared' column value in the output file.
        dictionary 'ionDict'
            key - ((float, string),int) tuple. See dictionary 'cosDict' key explanation in this function.
            value - set representing the indices of the query peaks that matched the library spectrum represented in this
                    dictionary value. The sum of these query peak intensities will be returned as the 'ionCount' column
                    value in the output file.
        dictionary 'ppmDict'
            key - ((float, string),int) tuple. See dictionary 'cosDict' key explanation in this function.
            value - list of floats represents the ppm difference of all peak matches that were within the given ppm tolerance. This is
                    most directly used for determining the ppm offset and standard deviation of an uncorrected
                    csodiaq output to be used in generating a corrected csodiaq output. For a corrected csodiaq output,
                    it is primarily used for figure generation to compare with uncorrected csodiaq output.
    tag 'qPPM':
        dictionary 'lightPPM' -
            'key' - string representing the peptide of interest.
            'value' - list of ppm differences between light library peak and a query peak.
        dictionary 'heavyPPM' -
            'key' - string representing the peptide of interest.
            'value' - list of ppm differences between heavy library peak and a query peak.
    tag 'ratio':
        dictionary 'lightRatio'
            'key' - string representing the peptide of interest.
            'value' - dictionary 'intensityRank'
                key - int representing intensity rank of the light library peak that matched a query peak.
                value - float representing the intensity of the matched query peak.
        dictionary 'heavyRatio'
            'key' - string representing the peptide of interest.
            'value' - dictionary 'intensityRank'
                key - int representing intensity rank of the heavy library peak that matched a query peak.
                value - float representing the intensity of the matched query peak.
'''
@njit()
def spectra_peak_comparison(libMzs, queMzs, ppmTol, ppmOffset):

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
Function: initialize_return_values()
Purpose: Initializes the return values for the peak comparison function. Because the same function is used several times for
            different purposes, this and the update_return_values() function allow for customized use of the same function.
            Explanations for return values are written under the spectra_peak_comparison() function, separated by tag.
Parameters:
    'tag' - string representing the condition under which the peak comparison is happening. This changes the return values.
Returns:
    list - contents depend on the value of 'tag'
'''
def initialize_return_values(tag):
    if tag=='identify':
        def cosDictValues(): return [0.0, 0.0, 0.0]
        cosDict = defaultdict(cosDictValues)
        countDict = defaultdict(int)
        ionDict = defaultdict(set)
        ppmDict = defaultdict(list)
        return [cosDict, countDict, ionDict, ppmDict]
    elif tag=='iPPM':
        def cosDictValues(): return [0.0, 0.0, 0.0]
        cosDict = defaultdict(cosDictValues)
        countDict = defaultdict(int)
        return [cosDict, countDict]
    elif tag=='qPPM':
        light = defaultdict(list)
        heavy = defaultdict(list)
        return [light, heavy]
    elif tag=='ratio':
        light = defaultdict(dict)
        heavy = defaultdict(dict)
        return [light, heavy]


'''
Function: update_return_values()
Purpose: Updates the return values for the peak comparison function. Because the same function is used several times for
            different purposes, this and the initialize_return_values() function allow for customized use of the same
            function. Explanations for each set of return values are written under the spectra_peak_comparison() function.
Parameters:
    'returns' - list of values. See initialize_return_values() for origin and spectra_peak_comparison() for contents.
    'peak1' - tuple in (float, float, tuple) format, representing a library peak. Contents of tuple are (mz, intensity,
        identifier).
    'peak2' - tuple in (float, float, tuple) format, representing a query peak. Contents of tuple are (mz, intensity,
        identifier).
    'i1' - index of the library peak.
    'i2' - index of the query peak.
    'ppm' - float representing the ppm difference between the mz values of the peaks.
    'tag' - string representing the condition under which the peak comparison is happening. This changes the return values.
Returns:
    no return value. Values in the 'return' parameter are updated.
'''
def update_return_values(returns, peak1, peak2, i1, i2, ppm, tag):
    if tag=='identify':
        key = (peak1[2],peak2[2])

        # Data required to calculate the cosine score
        returns[0][key][0] += peak1[1]*peak2[1]
        returns[0][key][1] += peak2[1]**2
        returns[0][key][2] += peak1[1]**2

        # Number of matching peaks
        returns[1][key] += 1

        # Ion Count
        returns[2][key].add(i2)

        # ppm differences in matched peaks
        returns[3][key].append(ppm)
    elif tag=='iPPM':
        key = (peak1[2],peak2[2])

        # Data required to calculate the cosine score
        returns[0][key][0] += peak1[1]*peak2[1]
        returns[0][key][1] += peak2[1]**2
        returns[0][key][2] += peak1[1]**2

        # Number of matching peaks
        returns[1][key] += 1

    elif tag=='qPPM':
        if peak1[2][0] == 'light': returns[0][peak1[2][2]].append(ppm)
        if peak1[2][0] == 'heavy': returns[1][peak1[2][2]].append(ppm)
    elif tag=='ratio':
        if peak1[2][0] == 'light': returns[0][peak1[2][2]][peak1[2][1]] = peak2[1]
        if peak1[2][0] == 'heavy': returns[1][peak1[2][2]][peak1[2][1]] = peak2[1]


'''
Function: cosine_similarity()
Purpose: Provided a value of the 'cosDict' dictionary output of the spectra_peak_comparison() function, this function
            explicitly calculates the cosine similarity score and returns it.
Parameters:
    'row' - list of various values. See dictionary 'cosDict' value explanation in function spectra_peak_comparison().
Returns:
    float - represents the cosine similarity score of a library and query spectrum comparison.
'''
def cosine_similarity(row):
    magnitude = (row[1]**0.5) * (row[2]**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
    return (row[0] / magnitude if magnitude else 0)


def generate_valid_FDR_spectra_keys(inFile):
    print('enter valid FDR spectra calculation:', flush=True)
    print(str(timedelta(seconds=timer())), flush=True)
    #if corrected: inFile = re.sub('(.*).csv', r'\1_corrected.csv', inFile)
    fields = ['decoy','MaCC_Score']
    df = pd.read_csv(inFile, sep=',', usecols=fields).sort_values('MaCC_Score', ascending=False)#.reset_index(drop=True)
    hits, decoys = fdr_calculation(df, returnType=2)
    spectraKeys = defaultdict(list)
    for h in hits:
        line = linecache.getline(inFile,int(h)+2).strip().split(',')
        spectraKeys[line[0]].append((float(line[2]),line[1]))
    return spectraKeys


'''
Function: pooled_library_query_spectra_analysis()
Purpose: This function loops through all query spectra and calculates the cosine similarity score and other
            values between it and every library spectra with a precursor m/z value within its designated window.
            For an explanation of each column in the output (references as the csodiaq output file in other
            comments), see the comments alongside the 'column' value in the function.
            NOTE: This function is used in two conditions to reduce the memory use. Conditions are determined
            by the presence of the spectralKeys parameter. The first condition (no spectralKeys parameter)
            determines the PSM values above an FDR cutoff of 0.01%, while the second calculates the various
            values for those PSMs that will be in the final output.
Parameters:
    'expSpectraFile' - string representing the path to the query spectra file (required .mzXML format).
    'lib' - dictionary as returned by the library_file_to_dict() function.
    'ppmTol' - see 'ppmTol' parameter description for function approx().
    'ppmYOffset' - see 'ppmYOffset' parameter description for function spectra_peak_comparison().
    'queryPooling' - int determining the maximum number of query spectra that can be pooled at any given time.
    'spectraKeys' - a dictionary that indicates the library/query spectra of interest as created by the
Returns:
    No Return Value. Results are written directly to the output file and ppm file (outFile and ppmFile,
        respectively). The description of specific columns of output file are provided in the function comments.
'''
def pooled_spectra_analysis(expSpectraFile, outFile, lib, ppmTol, ppmYOffset, queryPooling, spectraKeys=None):

    #
    if spectraKeys: tag='identify'
    else: tag='iPPM'

    # Column headers for the output file are initialized.
    if tag=='iPPM':
        columns = [
            'scan', # Scan number, corresponding to scans in the query spectra file.
            'peptide', # Peptide sequence for the library spectrum corresponding to this row.
            'mzLIB', # precursor charge for the library spectrum corresponding to this row.
            'decoy', # Protein name the peptide corresponds to, also derived from the library spectrum corresponding to this row.
            'MaCC_Score' # score unique to CsoDIAq, the fifith root of the number of matches ('shared') multiplied by the cosine score ('cosine')
        ]
    elif tag=='identify':
        columns = [
            'fileName', # Name of the query spectra file.
            'scan', # Scan number, corresponding to scans in the query spectra file.
            'MzEXP', # precursor m/z for query spectrum. Column 'windowWideness' corresponds to this value.
            'peptide', # Peptide sequence for the library spectrum corresponding to this row.
            'protein', # Protein name the peptide corresponds to, also derived from the library spectrum corresponding to this row.
            'MzLIB', # precursor m/z for the library spectrum corresponding to this row.
            'zLIB', # precursor charge for the library spectrum corresponding to this row.
            'cosine', # Cosine score comparing the library spectrum corresponding to this row with the query spectrum.
            'name', # Title - corresponds to the column "transition_group_id," a library spectrum identifier.
            'Peak(Query)', # The number of peaks in the query spectrum.
            'Peaks(Library)', # The number of peaks in the library spectrum.
            'shared', # The number of peaks that matched between query spectrum/library spectrum.
            'ionCount', # Sum of query spectrum intensities, excluding possible duplicates - currently uncalculated, set to 0.
            'CompensationVoltage', # The compensation voltage of the query spectrum.
            'totalWindowWidth', # width of m/z that was captured in the query spectrum. Corresponds to MzEXP.
            'MaCC_Score',# score unique to CsoDIAq, the fifith root of the number of matches ('shared') multiplied by the cosine score ('cosine')
        ]

    # query data file is loaded
    with mzxml.read(expSpectraFile, use_index=True) as spectra:


        # query data is looped over and scans are grouped by mz windows for future pooling.
        #  Additionally, for the second condition, various variables are saved for future reference.
        queScanDict = defaultdict(list)
        queValDict = defaultdict(dict)
        allLibKeys = set()
        spectralCount = 0
        if tag=='iPPM':
            for spec in spectra:
                if 'precursorMz' not in spec: continue
                queScanDict[spec['precursorMz'][0]['precursorMz'],spec['precursorMz'][0]['windowWideness']].append(spec['num'])
                spectralCount += 1
        elif tag=='identify':
            for scan in spectraKeys.keys():
                allLibKeys.update(spectraKeys[scan])
                spec = spectra.get_by_id(scan)
                queScanDict[spec['precursorMz'][0]['precursorMz'],spec['precursorMz'][0]['windowWideness']].append(spec['num'])
                #if scan=='1199' or scan=='1200': print((spec['precursorMz'][0]['precursorMz'],spec['precursorMz'][0]['windowWideness']))
                peaksCount = spec['peaksCount']
                if 'compensationVoltage' in spec: CV = spec['compensationVoltage']
                else: CV = ''
                queValDict[scan]['peaksCount'] = peaksCount
                queValDict[scan]['CV'] = CV
                spectralCount += 1

        print('Number of Unpooled MS/MS Query Spectra: ' + str(spectralCount))
        print('Number of Pooled MS/MS Query Spectra/Mz Windows: ' + str(len(queScanDict)),flush=True)

        # To enhance the print experience, status prints will be given at intervals tailored to the number of identified windows.
        #  example: if there are 1-99 pooled query spectra, print statements are made after every pooled query spectra analysis is complete.
        #           if there are 100-999, print after every 10 pooled spectra. And so on.
        printCutoff = 100
        while printCutoff < len(queScanDict): printCutoff*=10
        printCutoff /= 100

        # outfile is opened in advance so results can be written directly to the file as they are produced (mitigating memory use)
        with open(outFile, 'w', newline='') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerow(columns)

            # The second condition of this function returns a list of PPM differences that can be used for correction. Initialized here
            ppmList = []

            # Count variable keeps track of the number of query spectra that have been analyzed for time tracking purposes.
            count = 0

            # 'lib' dictionary keys are kept as a separate list in this analysis. Note that they are sorted by precursor m/z.
            if tag=='iPPM': allLibKeys = lib.keys()
            allLibKeys = sorted(allLibKeys)

            # Library keys were saved as an integer to save on time and simplify other parts of the algorithm. I believe the purpose is now defunct, but it's harmless, so I'm keeping it.
            idToKeyDict = {}
            for key in allLibKeys: idToKeyDict[lib[key]['ID']] = key

            # tracking time for print statements.
            prevtime = timer()

            print('Enter Pooled Spectra Analysis:')
            print(str(timedelta(seconds=prevtime)), flush=True)

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
                for i in range(len(scans)):

                    # adding each scan in the window to the pooled spectra
                    scan = scans[i]
                    spec = spectra.get_by_id(scan)
                    intensity = [x**0.5 for x in spec['intensity array']]
                    peakIDs = [scan for x in range(spec['peaksCount'])]
                    pooledQueSpectra += list(zip(spec['m/z array'],intensity,peakIDs))

                    # to reduce memory use for particularly large files, the user can limit the number of query spectra that are pooled. That's what this conditional statement takes care of.
                    if (i % queryPooling == 0 and i!=0) or i == len(scans)-1:
                        pooledQueSpectra.sort()

                        # each peak in the pooled library and query spectra is compared, and the necessary data is extracted from matching peaks (as determined by the ppm tolerance)
                        libMzs = np.array([x[0] for x in pooledLibSpectra])
                        queMzs = np.array([x[0] for x in pooledQueSpectra])
                        returns = initialize_return_values(tag)
                        pMatch, jMatch, ppmMatch = spectra_peak_comparison(libMzs, queMzs, ppmTol, ppmYOffset)
                        for i in range(len(pMatch)):
                            update_return_values(returns, pooledLibSpectra[pMatch[i]], pooledQueSpectra[jMatch[i]], pMatch[i], jMatch[i], ppmMatch[i], tag)

                        # output values are saved to the output file
                        if tag=='iPPM':
                            cosDict, countDict = returns
                            for key,value in cosDict.items():
                                # Library spectra that had too few matching peaks are excluded. numPeakMatch variable determines the threshold.
                                if countDict[key] > 2:
                                    cosine = cosine_similarity(cosDict[key])
                                    libKey = idToKeyDict[key[0]]
                                    scan = str(int(key[1]))
                                    if 'DECOY' in lib[libKey]['ProteinName']: decoy=1
                                    else: decoy=0
                                    temp = [
                                        scan, #scan
                                        libKey[1], #peptide
                                        libKey[0], #MzLIB
                                        decoy,
                                        (countDict[key]**(1/5))*cosine
                                    ]
                                    writer.writerow(temp)
                        if tag=='identify':
                            cosDict, countDict, ionDict, ppmDict = returns
                            for key,value in cosDict.items():
                                libKey = idToKeyDict[key[0]]

                                # Library spectra that had too few matching peaks are excluded. numPeakMatch variable determines the threshold.
                                if (libKey[0],libKey[1]) in spectraKeys[key[1]]:
                                    scan = str(int(key[1]))
                                    #print(scan)
                                    #print(sorted(queValDict)[:10])
                                    cosine = cosine_similarity(value)
                                    ionCount = sum([ pooledQueSpectra[j][1]+ppm_offset(pooledQueSpectra[j][1],ppmYOffset) for j in ionDict[key] ])
                                    temp = [
                                        expSpectraFile, #fileName
                                        scan, #scan
                                        precMz_win[0], #MzEXP
                                        libKey[1], #peptide
                                        lib[libKey]['ProteinName'], #protein
                                        libKey[0], #MzLIB
                                        lib[libKey]['PrecursorCharge'], #zLIB
                                        cosine, #cosine
                                        lib[libKey]['transition_group_id'], #name
                                        queValDict[scan]['peaksCount'], #Peaks(Query)
                                        len(lib[libKey]['Peaks']), #Peaks(Library)
                                        countDict[key], #shared
                                        ionCount, #ionCount
                                        queValDict[scan]['CV'], #compensationVoltage
                                        precMz_win[1], #totalWindowWidth
                                        (countDict[key]**(1/5))*cosine
                                    ]
                                    writer.writerow(temp)
                                    ppmList += ppmDict[key]

                        pooledQueSpectra.clear()
                        del returns

                # print statements for the user to track progress.
                count += 1
                if count % printCutoff == 0:
                    time = timer()
                    print('\nNumber of Pooled Experimental Spectra Analyzed: ' + str(count))
                    print('Number of Spectra in Current Pooled Spectrum: ' + str(len(scans)))
                    print('Time Since Last Checkpoint: ' + str(round(time-prevtime,2)) + ' Seconds', flush=True)
                    prevtime = time

    # Prints the final number of experimental spectra analyzed.
    print('Total Time (seconds): ' + str(timer()))
    print('Count: '+str(count),flush=True)
    if tag=='identify': return ppmList


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# PSM Scoring: Scoring
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''
Function:read_ppm_file_to_dict() *****DEFUNCT*****
Purpose: After filtering the csodiaq output for the optimal minimum allowed peak match number, the returned
            dictionary is then used to create a comprehensive list of ppm differences of the output, which is
            then in turn used to calculate the ppm difference mean (offset) and ppm standard deviation (ppm
            tolerance) for future calculations, as well as generate relevant figures.
Parameters:
    'ppmFile' - string corresponding to the 'ppmFile' parameter of the query_spectra_analysis() function.
Returns:
    'ppmDict' - dictionary representing the ppm differences of a given row in the csodiaq output.

    Defunct Reasons: I no longer write ppm to an outfile. This was eliminated when I broke
    pooled_spectra_analysis() into two parts.
'''
def read_ppm_file_to_dict(ppmFile):
    ppmDict = {}
    with open(ppmFile, 'r') as r:
        csvReader = csv.reader(r)
        for row in csvReader:

            # The first three rows correspond to elements of a 'key', and everything after that is ppm differences.
            ppmDict[(int(row[0]),row[1],int(row[2]))] = row[3:]
    return ppmDict


'''
Function: write_ppm_spread() *****DEFUNCT*****
Purpose: Given an csodiaq output file and a ppm output file (both products of the query_spectra_analysis()
            function), the comprehensive list of ppm differences or "ppm spread" list is created and written
            to a file. This accounts for filtering for the optimal number of allowed matched peaks.
Parameters:
    'cosFile' - corresponds to 'outFile' parameter of query_spectra_analysis() function, AFTER having been
        filtered for the optimal minimum allowed number of peak matches.
    'ppmFile' - corresponds to 'ppmFile' parameter of query_spectra_analysis() function.
    'outFile' - a string representing the path to the ppm spread file.
Returns:
    No Return Value. Results are written directly to the output file ('outFile').

    Defunct Reason: See read_ppm_file_to_dict() defunct reason
'''
def return_ppm_spread(df, ppmFile):
    # Data is read in.
    df.fillna('',inplace=True)

    ppmDict = read_ppm_file_to_dict(ppmFile)

    # list of keys corresponding to ppmDict are generated from the csodiaq data frame.
    listOfKeys = [(df['scan'].loc[i],df['peptide'].loc[i],df['zLIB'].loc[i]) for i in range(len(df))]

    # all values from the ppmFile corresponding to those keys are saved into a single list.
    ppmList = []
    for key in listOfKeys: ppmList += [float(x) for x in ppmDict[key]]

    return ppmList


'''
Function: find_offset_tol()
Purpose: Given a list of ppm differences from peaks that matched in a cosine similarity score analysis, this function
            determines the offset and tolerance that can be used in correcting a second analysis. Generally when seeing the
            list of ppm differences from an uncorrected analysis as a histogram you will see the main peak off-center from
            "0", suggesting all query spectra results were almost identically biased due to machine error. Uncorrected
            analyses have no offset (offset=0) and a wide tolerance to determine the "peak" in this histogram for corrective
            analysis.
        Also, if provided a string corresponding to the file path and file name of a desired histogram, the function will
            also write a histogram for visual confirmation.
        NOTE: We used the mean and second standard deviation as the offset and tolerance, respectively. However, other
            means can be used to determine all values that lie in the "peak" shown in the histogram.
Parameters:
    'data' - a list of floats. Each float corresponds to a library-query peak match in the previous analysis, depicting the
        ppm difference between the two (which was within the ppm tolerance of the analysis).
    'histFile' - string corresponding to the file path and file name of a desired histogram. If histFile == 0, no histogram
        will be written.
    'stdev' - The number of standard deviations of data to use as a tolerance. If 0, a custom tolerance is made of the
        histogram peak formed from the data is returned.
    'mean' - determines if the offset is determined by mean or median. I think this is defunct, and is overruled as
        custom value in the case of custom tolerance anyways.
Returns:
    'offset' - float corresponding to the calculated ppm offset to be used in future corrected analyses.
    'tolerance' - float corresponding to the calculated ppm tolerance to be used in future corrected analyses.
'''
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
        pyplot.savefig(histFile)
    return offset, tolerance


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# FDR Filtering
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''
Function: fdr_calculation()
Purpose: Given the output of the query_spectra_analysis() function as contained in a csodiaq output file and an
            FDR cutoff point, this function returns the number of rows in the dataframe above the FDR cutoff
            point. Note that this function is usually applied to various subsets of the actual csodiaq output
            file after filtering for minimum allowed peak matches. A higher return value from this function
            indicates a more optimal minimum allowed peak match number.
Parameters:
    'df' - Pandas Data Frame representing the contents of a csodiaq output file.
    'FDRCutoff' - float representing the minumum allowed FDR of the csodiaq output.
Returns:
    'fdr' - int representing number of rows in dataframe 'df' that are above the 'FDRCutoff' FDR cutoff point.
    'numDecoys' - int representing the number of rows represented by decoys specifically. NOTE: This return
        value is not currently used by the program, but is being kept here for possible future use.
'''
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


'''
Function: write_csodiaq_fdr_outputs()
Purpose: This function takes the output of the query_spectra_analysis() function, calculates the spectral, peptide, and
            protein FDR (removing values below an FDR threshold of 0.01), then writes the output of each to it's own file.
            File names are provided as parameters to the function. Each output has the same format as the
            query_spectra_analysis() function output, but with columns for FDR values added and rows below the FDR cutoff
            removed.
Parameters:
    'inFile' - string representing file path and name to the query_spectra_analysis() function output. Data from said file is
        read in and analyzed in this function.
    'specFile' - string representing file path and name to expected output of the spectral FDR. Spectral FDR column is
        included.
    'pepFile' - string representing file path and name to expected output of the peptide FDR. Peptide FDR column is included.
    'protFile' - string representing file path and name to expected output of the protein FDR. Protein FDR, Protein cosine
        score and peptide FDR columns are included. The peptide FDR column in this file is not expected to match what is
        written to 'pepFile', as the minimum number of allowed peaks is based on the protein FDR calculation rather than the
        optimal peptide FDR calculation used for the 'pepFile' output.
Returns:
    No return value. Data is written directly to files provided as string parameters.
'''
def write_csodiaq_fdr_outputs(inFile, specFile, pepFile, protFile):

    # cosine score data is read in
    overallDf = pd.read_csv(inFile).sort_values('MaCC_Score', ascending=False).reset_index(drop=True)
    #overallDf = pd.read_csv(inFile).sort_values('cosine', ascending=False).reset_index(drop=True)

    # spectral FDR is calculated and written to dataframe 'spectralDf'
    spectralDf = add_fdr_to_csodiaq_output(overallDf)

    # peptide FDR is calculated and written to dataframe 'peptideDf'
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    # Data from all of the above dataframes are written to their respective files.
    spectralDf.to_csv(specFile, index=False)
    peptideDf.to_csv(pepFile, index=False)

    if protFile:
        # Connections from peptides-proteins are listed in a file as (string, string) tuples.
        peptideProteinConnections = []

        for i in range(len(peptideDf)):
            peptide = peptideDf['peptide'].loc[i]

            # Notably, the protein group from the query_spectra_analysis() function is essentially a list of proteins the peptide is connected to.
            #   Thus, a connection is added for every protein in these protein groups.
            proteinGroup = peptideDf['protein'].loc[i]

            # a regular expression is used to separate proteins in the protein group
            #proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
            proteins = proteinGroup.split('/')[1:]
            for pro in proteins:

                # For decoys, 'pro[0]' will be added as the decoy tag. For non-decoys, 'pro[0]' is blank and therefore adds nothing to the protein name.
                #protein = pro[0] + pro[1]
                #peptideProteinConnections.append((peptide,protein))
                peptideProteinConnections.append((peptide,pro))

        # valid proteins are identified using the IDPicker algorithm
        verifiedProteinDict = idp.find_valid_proteins(peptideProteinConnections)

        # for each protein in the verified list, add all connected peptides found above the peptideFDR cutoff as a new dataframe.
        #   Note that this means the peptide can appear multiple times if found in more than one protein group provided by the IDPicker algorithm.
        proteinDf = add_leading_protein_column(peptideDf, verifiedProteinDict)

        tempProtFile = re.sub('(.*).csv', r'\1_delete.csv', protFile)
        proteinDf.to_csv(tempProtFile)

        # Protein FDR is calculated using the highest-scoring peptide for each protein group.
        tempProtDf = add_fdr_to_csodiaq_output(proteinDf, filterType='leadingProtein')
        proteinDict = tempProtDf.set_index('leadingProtein').T.to_dict()

        # Peptides that don't map to a protein above the FDR cutoff are excluded (indices of rows to be removed are added to this list)
        removables = []

        # protein cosine scores are included as a new column in the output.
        proteinCosine = []

        # protein FDR scores are included as a new column in the output.
        proteinFDR = []

         # Loops for every peptide in the recalculated peptide FDR dataframe.
        for i in range(len(proteinDf)):

            # if the leading protein group is one of the protein groups above the FDR cutoff point, protein cosine and FDR are added.
            protein = proteinDf['leadingProtein'].loc[i]
            if protein in proteinDict:
                proteinCosine.append(proteinDict[protein]['cosine'])
                proteinFDR.append(proteinDict[protein]['leadingProteinFDR'])

            # if leading protein group is NOT one of the protein groups above the FDR cutoff point, it is marked to be removed.
            else:
                removables.append(i)

        # invalid peptides are removed
        proteinDf = proteinDf.drop(proteinDf.index[removables]).reset_index(drop=True)

        # protein cosine/FDR scores are added as new columns
        proteinDf['proteinCosine'] = proteinCosine
        proteinDf['leadingProteinFDR'] = proteinFDR

        # for readability, the output is sorted by peptide cosine score, then leading protein, then protein cosine score
        #   In this way, you see a dataframe that is primarily sorted by proteins, and peptides inside the protein are ordered by cosine score
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
        proteinDf.to_csv(protFile, index=False)



'''
Function: add_fdr_to_csodiaq_output()
Purpose: This function calculates the FDR value for each row in a dataframe, filters out rows below an FDR rate of 0.01, and
            adds the final FDR values as a new column.
Parameters:
    'df' - Pandas dataframe with the same basic format as the query_spectra_analysis() function output
    'filterType' - string that determines if the dataframe should only use the first unique value (highest cosine score) from
     a given column. Also determines the name of the FDR column (filterType + 'FDR'). Default is set to 'spectral', where no
     filtering would be applied and the new column name will be 'spectralFDR'.
    'bestMatchNum' - int that determines the minimum number of matches allowed in calculating the FDR. Generally no value is
        provided, in which case the optimal minimum number of matches allowed is calculated in-function. This parameter
        exists primarily to allow for setting the number when calculating the protein-specific peptide FDR.
Returns:
    'finalDf' - Pandas dataframe with rows below an FDR rate of 0.01 removed and with a new FDR column added.
'''
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


'''
Function: add_leading_protein_column()
Purpose: Given a list of proteins verified by the IDPicker algorithm, this function adds the a new column 'leadingProtein'
            that indicates the protein group connected to the peptide represented by the row. Note that, for peptides that
            are found in multiple valid peptide groups, a new row will be added for each protein group connection.
Parameters:
    'df' - Pandas dataframe with the same basic format as the query_spectra_analysis() function output, though columns such
        as FDR may have been added.
    'verifiedProteinDict' - see idpicker.py find_valid_proteins() function return value 'finalDict' description.
Returns:
    'finalDf' - Dataframe with the new 'leadingProtein' column added.
'''
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


'''
Function: prepare_DIA_rerun_data()
Purpose: This function takes a protein FDR dataframe and converts it into a format that an MS machine can read for targetted
            DISPA re-analysis.
Parameters:
    'df' - Pandas dataframe, representing the protein FDR output from write_csodiaq_fdr_outputs().
    'trypsin' - boolean variable that indicates if identified peptides should be filtered for those that resulted from trypsin
        treatment.
Returns:
    'finalDf' - Pandas dataframe. 'finalDf' can be written to a file for direct input to an MS machine.
'''
def return_DISPA_targeted_reanalysis_dfs(header, inFile, proteins, trypsin, heavy):
    df = pd.read_csv(inFile)
    df = df[~df['protein'].str.contains('DECOY',na=False)].reset_index(drop=True)
    # Necessary?
    if trypsin: df = df[df['peptide'].str.endswith('R') | df['peptide'].str.endswith('K')].reset_index(drop=True)
    if 'uniquePeptide' in df.columns: df = df[df['uniquePeptide']==1].sort_values('ionCount', ascending=False).reset_index(drop=True)
    if proteins: df = df.groupby(['leadingProtein']).head(proteins).reset_index()

    CVs = set(df['CompensationVoltage'])
    def notNan(x): return ~np.isnan(x)
    CVs = set(filter(notNan, CVs))

    allDfs = []
    if len(CVs)==0: CVs.add('')

    for CV in CVs:
        if CV: tempDf = df[df['CompensationVoltage']==CV].reset_index(drop=True)
        else: tempDf = df
        #if proteins: tempDf = tempDf.groupby(['leadingProtein']).head(proteins).reset_index()
        lightMzs = []
        heavyMzs = []
        peptides = []
        rows = len(tempDf)
        for index, row in tempDf.iterrows():
            peptide = row['peptide']
            lightMz = float(row['MzLIB'])
            charge = row['zLIB']
            heavyMz = calc_heavy_mz(peptide, lightMz, charge)
            lightMzs.append(lightMz)
            heavyMzs.append(heavyMz)
            peptides.append(peptide)

        binWidth = 1.0
        lv, lBins = bin_assignment(lightMzs, binWidth)
        hv, hBins = bin_assignment(heavyMzs, binWidth)
        scanLightMzs = [lBins[lv[i]] for i in range(rows)]
        scanHeavyMzs = [hBins[hv[i]] for i in range(rows)]

        tempDf['scanLightMzs'] = [x-(binWidth/2) for x in scanLightMzs]
        tempDf['scanHeavyMzs'] = [x-(binWidth/2) for x in scanHeavyMzs]

        #tempDf.to_csv('/Users/calebcranney/Desktop/0_DataFiles/delete'+str(CV)+'.csv',index=False)
        allDfs.append(tempDf)

        binDict = defaultdict(list)
        for i in range(rows):
            if heavy: binDict[scanLightMzs[i],scanHeavyMzs[i]].append(peptides[i])
            else: binDict[scanLightMzs[i],0].append(peptides[i])

        data = []
        binKeys = sorted(binDict)
        for i in range(len(binKeys)):
            data.append(['/'.join(binDict[binKeys[i]]), '', '(no adduct)', binKeys[i][0]-(binWidth/2), 2, i+1])
            if heavy: data.append(['/'.join(binDict[binKeys[i]]), '', '(no adduct)', binKeys[i][1]-(binWidth/2), 2, i+1])

        finalDf = pd.DataFrame(data, columns=['Compound','Formula','Adduct','m.z','z','MSXID'])
        outFile = header
        if CV: outFile += '_' + str(CV)
        outFile += '.txt'
        finalDf.to_csv(outFile, sep='\t', index=False)
    #if len(CVs) > 1:
    compDf = pd.concat(allDfs)
    compDf.to_csv(header+'allCVs.csv', index=False)




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
            Rather than being be tagged by a library-specific tuple, the tags are tuples of (str, int) composition. The
            string differentiates between 'light' and 'heavy' peaks, and the int represents the rank of peak intensity (0
            being high).
    'mzxmlFiles' - list of strings representing the path to files from the DISPA re-analysis.
    'var' - @DEBUG: It is expected that this parameter will be removed by the end of this project. It currently holds
        variables that will be constant, but for testing/automating purposes are putting then in this list.
Returns:
    'finalDf' - pandas dataframe. The first two columns of the dataframe are the scans and peptides, respectively. Each
        subsequent column represents the ratios for those peptides derived from each data file of the DISPA re-analysis.

'''
def heavy_light_quantification(fragDict, libDict, mzxmlFiles, outDir, massTol, minMatch, ratioType, correction, hist):

    # final return value initialized, and the first two columns set.
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

    # Heavy:Light Ratio for each peptide is calculated for every DISPA reanalysis file. Looping over the files begins here.
    for f in mzxmlFiles:
        if correction == -1: offset, tolerance = 0, massTol
        else:
            # In order to apply a ppm correction, optimal offset and tolerance for the data.
            ppmDiffs = []
            with mzxml.read(f, use_index =True) as file:
                for scan in sorted(libDict.keys()):

                    # experimental spectrum is formatted for peak comparison (see 'libDict' parameter value for peak description)
                    spec = file.get_by_id(scan)
                    spec['intensity array'] = [x for x in spec['intensity array']]
                    peakIDs = [scan for x in range(len(spec['m/z array']))]
                    expSpectrum = list(tuple(zip(spec['m/z array'],spec['intensity array'],peakIDs)))
                    expSpectrum.sort(key=lambda x:x[0])

                    libSpectra = sorted(libDict[scan])
                    libMzs = np.array([x[0] for x in libSpectra])
                    queMzs = np.array([x[0] for x in expSpectrum])
                    tag = 'qPPM'
                    returns = initialize_return_values(tag)
                    pMatch, jMatch, ppmMatch = spectra_peak_comparison(libMzs, queMzs, massTol, 0)
                    for i in range(len(pMatch)):
                        update_return_values(returns, libSpectra[pMatch[i]], expSpectrum[jMatch[i]], pMatch[i], jMatch[i], ppmMatch[i], tag)
                    l, h = returns

                    # all PPM differences from matches of interest are included for consideration
                    #l, h = spectra_peak_comparison(sorted(libDict[scan]), expSpectrum, massTol, 0, tag='qPPM')
                    for pep in l:
                        if len(l[pep]) >= minMatch: ppmDiffs += l[pep]
                    for pep in h:
                        if len(h[pep]) >= minMatch: ppmDiffs += h[pep]


            offset, tolerance = find_offset_tol(sorted(ppmDiffs), hist, stdev=correction)

        # Now that an offset and tolerance has been established specific to the data, the program runs the algorithm to calculate heavy:light ratios
        data = []
        ratioDict = {}
        ppmDiffs = []
        with mzxml.read(f, use_index =True) as file:
            for scan in sorted(libDict.keys()):

                # as above, the scan spectrum is formatted appropriately.
                spec = file.get_by_id(scan)
                spec['intensity array'] = [x for x in spec['intensity array']]
                peakIDs = [scan for x in range(len(spec['m/z array']))]
                expSpectrum = list(tuple(zip(spec['m/z array'],spec['intensity array'],peakIDs)))
                expSpectrum.sort(key=lambda x:x[0])
                libSpectra = sorted(libDict[scan])
                libMzs = np.array([x[0] for x in libSpectra])
                queMzs = np.array([x[0] for x in expSpectrum])
                tag = 'ratio'
                returns = initialize_return_values(tag)
                pMatch, jMatch, ppmMatch = spectra_peak_comparison(libMzs, tolerance, offset, 0)
                for i in range(len(pMatch)):
                    update_return_values(returns, libSpectra[pMatch[i]], expSpectrum[jMatch[i]], pMatch[i], jMatch[i], ppmMatch[i], tag)

                lightPeps, heavyPeps = returns

                #lightPeps, heavyPeps = spectra_peak_comparison(sorted(libDict[scan]), expSpectrum, tolerance, offset, tag='ratio')
                peps = [x['seq'] for x in fragDict[scan]]
                for pep in peps:
                    l = lightPeps[pep]
                    h = heavyPeps[pep]
                    # for light or heavy peaks that are identified without a match, a psuedo peak is created
                    #   with an intensity equal to the average of 10 least intense peaks divided by two.

                    smallPeakIntensity = np.mean(sorted(spec['intensity array'])[:10])/2

                    # creating the psuedo peaks
                    lightSet = set(l.keys())
                    heavySet = set(h.keys())
                    onlyLight = lightSet - heavySet
                    for key in onlyLight: h[key] = smallPeakIntensity
                    onlyHeavy = heavySet - lightSet
                    for key in onlyHeavy: l[key] = smallPeakIntensity

                    # ratio is calculated if there are enough peaks to warrant a match.
                    if minMatch: check = (len(l) >= minMatch)
                    else: check = match_score(list(l.keys()))
                    if check:
                        ratio = return_ratio(l, h, ratioType)

                    # if there are not enough peaks to warrant a match, the ratio is returned as NaN.
                    else:
                        ratio=np.nan

                    # This dictionary maps scans to calculated ratios.
                    ratioDict[scan, pep] = ratio

        # column of ratios from the just-analyzed file is added to the output
        finalDf[f] = [ratioDict[key] for key in sorted(ratioDict.keys())]
    return finalDf


def match_score(keys):
    if len(keys) == 0: return False
    if min(keys) < 3: return True
    return False

'''
Function: make_quant_dicts()
Purpose: This function prepares the dictionaries used by the heavy_light_quantification() function.
Parameters:
    'inFile' - string representing the path to output from the FDR calculation specific to the data that was used to target
        the DISPA re-analysis.
    'libFile' - string representing the path to the library spectra, ideally the same library that was used in steps leading
        up to FDR analysis.
    'mzxmlFiles' - list of strings representing the path to files from the DISPA re-analysis.
    'maxPeaks' - maximum number of peaks to be included in each library (light and heavy combined will have at most double
        this number)
Returns:
    'fragVarDict' - see parameter 'fragDict' in heavy_light_quantification() function.
    'libPeakDict' - see parameter 'libDict' in heavy_light_quantification() function.
'''
def make_quant_dicts(inFile, libFile, mzxmlFiles, maxPeaks):
    fileType = inFile.split('.')[-1]
    if fileType == 'csv': fragDf = pd.read_csv(inFile)
    else: fragDf = pd.read_csv(inFile, sep='\t')

    # Dictionary is created that can match FDR calculation variables to the scan number of the DISPA re-analysis.
    # NOTE: all DISPA re-analysis files should have corresponding scan numbers. The first of the list is arbitrarily chosen there.
    fragScanDict = {}

    with mzxml.read(mzxmlFiles[0]) as spectra:
        for x in spectra:

            # @MAKE NOTE
            key = ( round(x['precursorMz'][0]['precursorMz'], 2),
                    round(x['precursorMz'][1]['precursorMz'], 2),
                    x['compensationVoltage']
                    )
            fragScanDict[key] = x['num']

    fragVarDict = defaultdict(list)
    libScanDict = {}
    print('total scans: ' + str(len(fragScanDict)))
    # each peptide of interest should have a corresponding scan. This runs over each peptide.
    for i in range(len(fragDf)):

        # @DEBUG: Top line for my ouput, bottom for Jesse's old output.
        seq, mz, z, CV = fragDf.loc[i]['peptide'], fragDf.loc[i]['MzLIB'], fragDf.loc[i]['zLIB'], fragDf.loc[i]['CompensationVoltage']
        #seq, mz, z, CV = fragDf.loc[i]['Peptide'], fragDf.loc[i]['prec_light_mz'], fragDf.loc[i]['z'], fragDf.loc[i]['CV']
        lightMz = round(mz, 2)

        key = (round(fragDf.loc[i]['scanLightMzs'],2), round(fragDf.loc[i]['scanHeavyMzs'],2), CV)
        #key = (lightMz, heavyMz, CV)
        if key in fragScanDict:
            scan = fragScanDict[key]
            libScanDict[lightMz, seq] = scan
            fragVarDict[scan].append({'seq':seq, 'mz':mz, 'z':z, 'CV':CV})


    print('scans identified: '+str(len(libScanDict)))
    # dictionary with a scan:peaks entries is created.
    libPeakDict = make_lib_dict_quant(libFile, libScanDict, fragDf, maxPeaks)
    print('libraries identified: '+str(len(libPeakDict)))
    return fragVarDict, libPeakDict


'''
Function: tuple_key_match() *****DEFUNCT*****
Purpose: Given a list of tuple keys and a key that may be in the list, this function determines if an approximately equal key
            exists. String values must be exact, but float values must be withit 10 ppm to be considered a match.
Parameters:
    'keyList' - a list of tuples. Composition of the tuples in the list should be consistent with each other.
    'key' - a tuple of the same composition as the tuples in keyList.
Returns:
    boolean - if 'key' is found in 'keyList', this function returns true, else false.

Defunct Reasons: We've opted to just round float values to 2 decimal places rather than finding values within a ppm
                    tolerance. At the time we chose not to use this function, I was debugging instances where 1) a different
                    CV match was being chosen first and 2) the wrong heavy-mz values were being identified, possibly because
                    the light-mz value (ordered before heavy) was different.
'''
def tuple_key_match(keyList, key):
    i = bisect(keyList, key)
    if i != 0 and approx_tuple(key, keyList[i-1]): return keyList[i-1]
    if i != len(keyList) and approx_tuple(key, keyList[i]): return keyList[i]
    if i<3 and i < len(keyList)-4:
        print('i: '+str(i))
        for j in range(i-2,i+2):
            print(str(j) + str(keyList[i-1]))
        print(key)
    return False


'''
Function: approx_tuple() *****DEFUNCT*****
Purpose: This function compares each value of two tuples to determine if they are approximately (in the case of float values)
            or exactly (in the case of string or int values) the same.
Parameters:
    't1' - a tuple.
    't2' - a tuple.
    'ppmTol' - the ppm tolerance for float values to be considered a match. Default is 10ppm.
Returns:
    boolean - if each value of tuples match this function returns true, else false.

Defunct Reasons: We've opted to just round float values to 2 decimal places rather than finding values within a ppm
                    tolerance.
'''
def approx_tuple(t1, t2, ppmTol=10):

    # if the tuples aren't the same length, return False. This shouldn't happen, I'm just being careful.
    if len(t1) != len(t2): return False
    l = []
    for i in range(len(t1)):
        # if the tuple values aren't the same type, return False. This also shouldn't happen.
        if type(t1[i]) != type(t2[i]): return False

        # if the values are strings and don't match, return False.
        if type(t1[i]) is str and t1[i] != t2[i]: return False

        # If the values are numeric and not within the ppm tolerance, return False.
        if not approx(t1[i], t2[i], ppmTol): return False

    return True

'''
Function: make_lib_dict_quant()
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
    'fragDf' - Pandas dataframe, created from 'inFile' parameter of make_quant_dicts() function.
    'maxPeaks' - maximum number of peaks to be included in each library (light and heavy combined will have at most double
        this number)
Returns:
    dictionary - see parameter 'libDict' in heavy_light_quantification() function.

'''
def make_lib_dict_quant(libFile, libScanDict, fragDf, maxPeaks):
    fileType = libFile.split('.')[-1]
    if fileType == 'mgf':

        # all code before the return statement assumes that the peptides have decimal values corresponding to non-standard modifications.
        uniquePeps = set(fragDf['Peptide'])
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


'''
Function: mgf_library_upload_quant()
Purpose: This function creates a dictionary with scan:peak entries for identifying library spectra based on scan number. This
            function is specific to MGF library files, but is specific to quantification, distinguishing it from the
            mgf_library_upload() function.
Parameters:
    'fileName' - string representing the path to the library spectra, ideally the same library that was used in steps leading
        up to FDR analysis.
    dictionary 'scanDict' - see parameter 'libScanDict' in make_lib_dict_quant() function.
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
            peaks.append((fragMz, fragInt, ('light', i, seq)))
            peaks.append((calc_heavy_mz(fragList[i][2], fragMz, 1), fragInt, ('heavy', i, seq)))

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
    dictionary 'scanDict' - see parameter 'libScanDict' in make_lib_dict_quant() function.
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

    # @DEBUG: for printing out results specific to a scan. To be removed by the end.
    scanEx = ''

    # Peaks list is created and attached to the dictionary
    finalDict = defaultdict(list)
    for key in lib:

        # y-ion peaks are sorted by intensity, and lower-intensity peaks are filtered out.
        fragList = sorted(list(tuple(zip(intensity_dict[key], mz_dict[key]))), reverse=True)
        # @DEBUG: see above
        if scanDict[key]==scanEx:
            for x in fragList: print(x)
        if maxPeaks !=0 and len(fragList) >= maxPeaks: fragList = fragList[:maxPeaks]

        # heavy counterpart mz is calculated. Light and heavy pairs are additionally tagged by their intensity rank and included in the final output.
        peaks = []
        for i in range(len(fragList)):
            fragMz = fragList[i][1]
            fragInt = fragList[i][0]
            peaks.append((fragMz, fragInt, ('light',i,key[1])))
            fragSeq = lib[key]['PeptideSequence'][-lib[key]['FragmentSeriesNumber']:]
            heavyMz = calc_heavy_mz(fragSeq, fragMz, lib[key]['FragmentCharge'])
            peaks.append((heavyMz, fragInt, ('heavy',i,key[1])))
        # @DEBUG: see above
        if scanDict[key]==scanEx:
            for x in peaks: print(x)

        # entry placed in final dictionary
        finalDict[scanDict[key]] += peaks

    return finalDict

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



#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#
# Miscellaneous Functions
#
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
'''
Function: set_plot_settings()
Purpose: Sets parameters for histogram graphics.
Parameters:
    'xlabel' - string indicating the x-axis label.
    'ylabel' - string indicating the y-axis label.
    'wide' - boolean for two different dimension possibilities
Returns:
    No return value. Just sets parameters using pyplot.
'''
def set_plot_settings(xlabel, ylabel, wide=True):
    if wide: pyplot.figure(figsize=(18,12)).add_axes([0.11, 0.1, 0.85, 0.85])
    else: pyplot.figure(figsize=(12,12))
    pyplot.axhline(linewidth=4, color='black')
    pyplot.axvline(linewidth=4, x=-10, color='black')
    pyplot.xlim(-10,10)
    pyplot.xlabel(xlabel, fontsize = 36, weight='bold')
    pyplot.ylabel(ylabel, fontsize = 36, weight='bold')
    pyplot.tick_params(axis="x", labelsize=36)
    pyplot.tick_params(axis="y", labelsize=36)

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
