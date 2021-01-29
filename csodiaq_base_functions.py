import pandas as pd
from pyteomics import mzxml, mgf, mass
from bisect import bisect
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as pyplot
import numpy as np
import idpicker as idp
import re
from collections import defaultdict
from os import listdir
from os.path import isfile, join


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
        lib = traml_library_upload_csv(inFile)
    return lib


'''
Function: traml_library_upload_csv()
Purpose: Given a traml file in .tsv or .csv format representing the library spectra, this function generates a dictionary to
            be used for future analysis. key:value pairings are described in the "Returns" section of this comment. Note that
            the tuple key format was chosen because library spectra will need to be sorted and chosen by precursor mass, so
            identifiers already included in the library spectra were explicitly not used as keys.
Parameters:
    'fileName' - string representing the path to the library spectra file (required traml .tsv or .csv format).
Returns:
    dictionary 'lib' - See dictionary 'lib' key explanation in function library_file_to_dict().
'''
def traml_library_upload_csv(fileName):
    # Library spectra file is read as a pandas dataframe - conditional statement allows for both .tsv and .csv files to be uploaded.
    if fileName.endswith('.tsv'):
        lib_df = pd.read_csv(fileName, sep='\t')
    else:
        lib_df = pd.read_csv(fileName)

    # Print statement for timing the program
    print("#Enter library dictionary upload: ")
    print('#'+str(timedelta(seconds=timer())))

    # Unneeded columns are removed from the dataframe
    lib_df = lib_df.loc[:, lib_df.columns.intersection(['PrecursorMz','FullUniModPeptideName','PrecursorCharge','ProductMz','LibraryIntensity','transition_group_id','ProteinName'])]

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

    # Peaks list is created and attached to the dictionary
    for key in lib:
        mz, intensity = (list(t) for t in zip(*sorted(zip(mz_dict[key], intensity_dict[key]))))
        keyList = [key for i in range(len(mz))]
        peaks = list(tuple(zip(mz,intensity,keyList)))
        peaks.sort(key=lambda x:x[0])
        lib[key]['Peaks'] = peaks
    return lib


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
    print('#Enter library dictionary upload: ')
    print('#'+str(timedelta(seconds=timer())))

    # return value is initialized
    lib = {}

    # each spectrum in the mgf file
    for spec in libMGF:

        # key for the final dictionary is initialized
        key = (spec['params']['pepmass'][0], spec['params']['seq'])

        # final dictionary value values are created
        charge = spec['params']['charge'][0]
        name = spec['params']['title']

        # note that the protein value is optional in MGF files - if it's nonexistent, protein value initialized as an empty string
        if 'protein' in spec['params']: protein = spec['params']['protein']
        else: protein = ''

        # peaks of the library file are intialized.
        mz = spec['m/z array']
        intensity = spec['intensity array']
        intensity = [x**0.5 for x in intensity]
        keyList = [key for x in mz]
        peaks = list(tuple(zip(mz,intensity,keyList)))
        peaks.sort(key=lambda x:x[0])

        # final dictionary value is created
        tempDict = {
            'PrecursorCharge':charge,
            'transition_group_id':name,
            'ProteinName':protein,
            'Peaks':peaks
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
    if i1 == len(temp)-1 and len(temp)!= 1: return []
    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
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
            If the difference is outside the provided tolerance, 0 is returned instead for conditional purposes.
'''
def approx(x, y, ppmTol):
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)


'''
Function: ppm_offset()
Purpose: Given a peak mz value and a ppm offset value, this function calculates the mz offset that would need to be
            subtracted from the mz value to align the final spread of ppm differences around 0.
Parameters:
    'mz' - Mz value of a given peak (generally a query peak)
    'ppm' - Ppm offset value. This is generally calculated from an uncorrected csodiaq output. During the uncorrected
        csodiaq output generation, this value is generally 0, which in turn makes the returned mz offset 0.
Returns:
    float - Mz offset value.
'''
def ppm_offset(mz, ppm):
    return (mz*ppm)/1000000


'''
Function: pooled_all_peak_comparison()
Purpose: This function determines peaks that are within a sufficiently close tolerance to one another. Various functions
            requiring peak comparison then take place inside this function, including compiling data for a cosine
            similarity score. The specifics of the various functions can be found in the return values below.
Parameters:
    'libSpectrum' - list of (float, float, tuple) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict().
    'expSpectrum' - list of (float, float, int) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict(). The only variation is that the int represents the scan number.
    'ppmTol' - see 'ppmTol' parameter description for function approx().
    'ppmYOffset' - The value in ppm that should be subtracted from the query peak m/z. This is the average ppm difference
        of an uncorrected csodiaq output, being applied to the experimental spectra to generate a corrected
        csodiaq calculation output.
Returns:
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
'''
def spectra_peak_comparison(libSpectrum, expSpectrum, ppmTol, ppmYOffset, tag='identify'):
    # final dictionary returned is initialized. See function description for details on contents.
    returns = initialize_return_values(tag)

    # By tracking the indices of the current library/query peaks we reduce the time complexity of the algorithm
    i, j = 0, 0
    expPeakMz = expSpectrum[j][0] - ppm_offset(expSpectrum[j][0], ppmYOffset)
    while i < len(libSpectrum) and j < len(expSpectrum):

        # If the m/z of the peaks are not within the given ppm tolerance, the indices of the smaller of the two is incremented
        #   and we are returned to the top of the while loop.
        if not approx(libSpectrum[i][0],expPeakMz, ppmTol):
            if libSpectrum[i][0] > expPeakMz:
                j += 1
                if j < len(expSpectrum): expPeakMz = expSpectrum[j][0] - ppm_offset(expSpectrum[j][0], ppmYOffset)
                continue
            if libSpectrum[i][0] < expPeakMz: i += 1; continue

        # To account for the matching of one query peak to multiple library peaks, library peaks are looped over
        #   after the initial match. Every matching peak contributes to the various variables in the final returned dictionary.
        p = i + 0
        while (p < len(libSpectrum)):
            ppm = approx(libSpectrum[p][0], expPeakMz, ppmTol)
            if p==len(libSpectrum) or not ppm: break
            update_return_values(returns, libSpectrum[p], expSpectrum[j], p, j, ppm, tag)

            p += 1

        # Note that the possibility of one library peak matching to multiple query peaks is automatically accounted for
        #   by the fact that the query peak is the next default increment after all match calculations have been made.
        j += 1
        if j < len(expSpectrum): expPeakMz = expSpectrum[j][0] - ppm_offset(expSpectrum[j][0], ppmYOffset)
    return returns


def initialize_return_values(tag):
    if tag=='identify':
        def cosDictValues(): return [0.0, 0.0, 0.0]
        cosDict = defaultdict(cosDictValues)
        countDict = defaultdict(int)
        ionDict = defaultdict(set)
        ppmDict = defaultdict(list)
        return [cosDict, countDict, ionDict, ppmDict]
    elif tag=='qPPM': return[[],[]]
    elif tag=='ratio':
        light = defaultdict(int)
        heavy = defaultdict(int)
        return [light, heavy]


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

    elif tag=='qPPM':
#        if peak1[2][0] == 'light': returns[0].append((peak2[1],peak2[0]))
#        if peak1[2][0] == 'heavy': returns[1].append((peak2[1],peak2[0]))
        if peak1[2][0] == 'light': returns[0].append(ppm)
        if peak1[2][0] == 'heavy': returns[1].append(ppm)
    elif tag=='ratio':
        if peak1[2][0] == 'light': returns[0][peak1[2][1]] = peak2[1]
        if peak1[2][0] == 'heavy': returns[1][peak1[2][1]] = peak2[1]


'''
Function: cosine_similarity()
Purpose: Provided a value of the 'cosDict' dictionary output of the peak_comparison() function, this function
            explicitly calculates the cosine similarity score and returns it.
Parameters:
    'row' - list of various values. See dictionary 'cosDict' value explanation in function peak_comparison().
Returns:
    float - represents the cosine similarity score of a library and query spectrum comparison.
'''
def cosine_similarity(row):
    magnitude = (row[1]**0.5) * (row[2]**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
    return (row[0] / magnitude if magnitude else 0)

'''
Function: pooled_library_query_spectra_analysis()
Purpose: This function loops through all query spectra and calculates the cosine similarity score and other
            values between it and every library spectra with a precursor m/z value within its designated window.
            For an explanation of each column in the output (references as the csodiaq output file in other
            comments), see the comments alongside the 'column' value in the function.
Parameters:
    'expSpectraFile' - string representing the path to the query spectra file (required .mzXML format).
    'outFile' - string representing the path to the output/results file.
    'ppmFile' - a string representing the path to the ppm difference results file. The first three columns
        are used for creating a key that corresponds to rows in the outFile output. This ppm file is used
        to compile a comprehensive list of ppm differences used in calculating the ppm offset and ppm
        standard deviation of an uncorrected csodiaq output. Such calculations are done after the optimal
        minimum peak matching number is determined and applied as a filter to the data.
    'lib' - dictionary as returned by the library_file_to_dict() function.
    'numPeakMatch' - int representing the minimum number of allowed peak matches. This number is generally 3
        to catch a wide number of matches that can later be filtered for the optimal number of minimum
        allowed peak matches.
    'ppmTol' - see 'ppmTol' parameter description for function approx().
    'ppmYOffset' - see 'ppmYOffset' parameter description for function peak_comparison().
Returns:
    No Return Value. Results are written directly to the output file and ppm file (outFile and ppmFile,
        respectively). The description of specific columns of output file are provided in the function comments.
'''
def pooled_all_query_spectra_analysis(expSpectraFile, outFile, ppmFile, lib, ppmTol, ppmYOffset):
    # Column headers for the output file are initialized.
    columns = [
        'fileName', # Name of the query spectra file.
        'scan', # Scan number, corresponding to scans in the query spectra file.
        'MzEXP', # precursor m/z for query spectrum. Column 'windowWideness' corresponds to this value.
        'zEXP', # precursor charge for query spectrum. Note that this is likely an estimate, as the exact charge of every compound in the scan is unknown.
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
        'MaCC_Score'# score unique to CsoDIAq, the fifith root of the number of matches ('shared') multiplied by the cosine score ('cosine')
    ]

    # Output file is opened and column headers are written as the first row.
    with open(outFile, 'w', newline='') as csvFile, open(ppmFile, 'w', newline='') as ppmFile:

        writer = csv.writer(csvFile)
        writer.writerow(columns)

        ppmWriter = csv.writer(ppmFile)

        # Count variable keeps track of the number of query spectra that have been analyzed for time tracking purposes.
        count = 0

        # 'lib' dictionary keys are kept as a separate list in this analysis. Note that they are sorted by precursor m/z.
        allLibKeys = sorted(lib)

        quePeakDict = defaultdict(list)

        queValDict = {}

        with mzxml.read(expSpectraFile) as spectra:
            for spec in spectra:
#                if int(spec['num']) < 3601:
                num = spec['num']
                precursorMz = spec['precursorMz'][0]['precursorMz']
                precursorCharge = spec['precursorMz'][0]['precursorCharge']
                peakCount = len(spec['intensity array'])
                CV = spec['compensationVoltage']
                window = spec['precursorMz'][0]['windowWideness']

                spec['intensity array'] = [x**0.5 for x in spec['intensity array']]
                peakIDs = [num for x in range(len(spec['m/z array']))]
                top_mz = precursorMz + window / 2
                bottom_mz = precursorMz - window / 2
                quePeakDict[(top_mz, bottom_mz)] += zip(spec['m/z array'],spec['intensity array'],peakIDs)

                queValDict[num] = [ precursorMz, precursorCharge, peakCount, CV, window ]

        time = timer()
        prevtime = time
        for w in quePeakDict:
            # Printing time taken to analyze every 100 spectra.
            count += 1
            if count % 100 == 0:
                time = timer()
                print(str(count)+','+str(time-prevtime)+','+str(len(spec['m/z array']))+','+str(spec['precursorMz'][0]['precursorMz'])+','+str(len(libKeys))+','+outFile)
                prevtime = time

            quePeakDict[w] = sorted(quePeakDict[w])
            libKeys = lib_mz_match_query_window( w[0], w[1], allLibKeys )

            if len(libKeys) != 0:
                cosDict, countDict, ionDict, ppmDict = spectra_peak_comparison(pool_lib_spectra(lib, libKeys), quePeakDict[w], ppmTol, ppmYOffset)
            else: continue
            for key in cosDict:
                # Library spectra that had too few matching peaks are excluded. numPeakMatch variable determines the threshold.

                if countDict[key] > 2:
                    cosine = cosine_similarity(cosDict[key])
                    ionCount = sum([ quePeakDict[w][j][1]+ppm_offset(quePeakDict[w][j][1],ppmYOffset) for j in ionDict[key] ])
                    temp = [
                        expSpectraFile, #fileName
                        key[1], #scan
                        queValDict[key[1]][0], #MzEXP
                        queValDict[key[1]][1], #zEXP
                        key[0][1], #peptide
                        lib[key[0]]['ProteinName'], #protein
                        key[0][0], #MzLIB
                        lib[key[0]]['PrecursorCharge'], #zLIB
                        cosine, #cosine
                        lib[key[0]]['transition_group_id'], #name
                        queValDict[key[1]][2], #Peaks(Query)
                        len(lib[key[0]]['Peaks']), #Peaks(Library)
                        countDict[key], #shared
                        ionCount, #ionCount
                        queValDict[key[1]][3], #compensationVoltage
                        queValDict[key[1]][4], #totalWindowWidth
                        (countDict[key]**(1/5))*cosine
                    ]
                    writer.writerow(temp)
                    ppmWriter.writerow([key[1], key[0][1], lib[key[0]]['ProteinName']] + sorted(ppmDict[key]))

    # Prints the final number of experimental spectra analyzed.
    print('#'+str(count))


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
Function:read_ppm_file_to_dict()
Purpose: After filtering the csodiaq output for the optimal minimum allowed peak match number, the returned
            dictionary is then used to create a comprehensive list of ppm differences of the output, which is
            then in turn used to calculate the ppm difference mean (offset) and ppm standard deviation (ppm
            tolerance) for future calculations, as well as generate relevant figures.
Parameters:
    'ppmFile' - string corresponding to the 'ppmFile' parameter of the query_spectra_analysis() function.
Returns:
    'ppmDict' - dictionary representing the ppm differences of a given row in the csodiaq output.
'''
def read_ppm_file_to_dict(ppmFile):
    ppmDict = {}
    with open(ppmFile, 'r') as r:
        csvReader = csv.reader(r)
        for row in csvReader:

            # The first three rows correspond to elements of a 'key', and everything after that is ppm differences.
            ppmDict[(int(row[0]),row[1],row[2])] = row[3:]
    return ppmDict


'''
Function: write_ppm_spread()
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
'''
def return_ppm_spread(df, ppmFile):
    # Data is read in.
    df.fillna('',inplace=True)

    ppmDict = read_ppm_file_to_dict(ppmFile)

    # list of keys corresponding to ppmDict are generated from the csodiaq data frame.
    listOfKeys = [(df['scan'].loc[i],df['peptide'].loc[i],df['protein'].loc[i]) for i in range(len(df))]

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
Returns:
    'offset' - float corresponding to the calculated ppm offset to be used in future corrected analyses.
    'tolerance' - float corresponding to the calculated ppm tolerance to be used in future corrected analyses.
'''
def find_offset_tol(data, histFile, stdev=2, mean=True):

    # offset is calculated as the mean value of the provided data
    if mean: offset = sum(data)/len(data)
    else: offset = data[len(data)//2]
    # tolerance is calculated as the second standard deviation of the provided data
    tolerance = statistics.pstdev(data)*stdev

    # if a histogram file is provided, it is created, with offset (black) and tolerance (red) lines drawn for reference
    if histFile:
        hist, bins = np.histogram(data, bins=200)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        pyplot.clf()
#        set_plot_settings('Difference between Matched Peaks (PPM)','Frequency')

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
def fdr_calculation(df, fdrList=False):
    # initializing the two return values at 0
    fdrValues = []
    numDecoys = 0
    df.fillna("nan",inplace=True)
    # for every row in the dataframe
    for i in range(len(df)):

        # current criteria for 'decoys' is to have 'decoy' in the protein name. This may change in the future.
        if 'DECOY' in df.loc[i]['protein']:
            numDecoys += 1

        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(i+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if len(fdrValues) < 1/0.01:
                if fdrList: return [], 0
                else: return 0, 0
            if fdrList: return fdrValues, numDecoys-1
            else: return len(fdrValues), numDecoys-1
        fdrValues.append(curFDR)
    if fdrList: return fdrValues, numDecoys-1
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

    # spectral FDR is calculated and written to dataframe 'spectralDf'
    spectralDf = add_fdr_to_csodiaq_output(overallDf)

    # peptide FDR is calculated and written to dataframe 'peptideDf'
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    # Connections from peptides-proteins are listed in a file as (string, string) tuples.
    peptideProteinConnections = []

    for i in range(len(peptideDf)):
        peptide = peptideDf['peptide'].loc[i]

        # Notably, the protein group from the query_spectra_analysis() function is essentially a list of proteins the peptide is connected to.
        #   Thus, a connection is added for every protein in these protein groups.
        proteinGroup = peptideDf['protein'].loc[i]

        # a regular expression is used to separate proteins in the protein group
        proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
        for pro in proteins:

            # For decoys, 'pro[0]' will be added as the decoy tag. For non-decoys, 'pro[0]' is blank and therefore adds nothing to the protein name.
            protein = pro[0] + pro[1]
            peptideProteinConnections.append((peptide,protein))

    # valid proteins are identified using the IDPicker algorithm
    verifiedProteinDict = idp.find_valid_proteins(peptideProteinConnections)

    # for each protein in the verified list, add all connected peptides found above the peptideFDR cutoff as a new dataframe.
    #   Note that this means the peptide can appear multiple times if found in more than one protein group provided by the IDPicker algorithm.
    proteinDf = add_leading_protein_column(peptideDf, verifiedProteinDict)

    # Protein FDR is calculated using the highest-scoring peptide for each protein group.
    tempProtDf = add_fdr_to_csodiaq_output(proteinDf, filterType='leadingProtein')

    tempProtDf.to_csv('Data/oldOutput/protein.csv')

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
        if len(uniquePepsDict[proteinDf.loc[i]['peptide']]) == 1: uniquePeps.append(1)
        else: uniquePeps.append(0)

    proteinDf['uniquePeptide'] = uniquePeps

    # Data from all of the above dataframes are written to their respective files.
    spectralDf.to_csv(specFile, index=False)
    peptideDf.to_csv(pepFile, index=False)
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
    fdrList, decoyNum = fdr_calculation(finalDf, fdrList=True)

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
        proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)

        # for each protein connected to this peptide
        for pro in proteins:

            # this ensures that decoys are tagged while non-decoys are not
            protein = pro[0] + pro[1]

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
    #nK = len(seq) - len(re.sub('K','',seq))
    #nR = len(seq) - len(re.sub('R','',seq))

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
def return_DISPA_targeted_reanalysis_dfs(header, inFile, proteins, trypsin):
    df = pd.read_csv(inFile)
    df = df[~df['protein'].str.contains('DECOY')].reset_index(drop=True)
    # Necessary?
    if trypsin: df = df[df['peptide'].str.endswith('R') | df['peptide'].str.endswith('K')].reset_index(drop=True)
    df = df[df['uniquePeptide']==1].sort_values('ionCount', ascending=False).reset_index(drop=True)

    CVs = set(df['CompensationVoltage'])
    for CV in CVs:
        tempDf = df[df['CompensationVoltage']==CV].reset_index(drop=True)
        if proteins: tempDf.groupby('leadingProtein').head(proteins).reset_index()

        data = []
        for i in range(len(tempDf)):
            compound = tempDf.loc[i]['peptide']
            formula = ''
            adduct = '(no adduct)'
            lightMz = float(tempDf.loc[i]['MzLIB'])
            charge = tempDf.loc[i]['zLIB']
            heavyMz = calc_heavy_mz(compound, lightMz, charge)
            MSXID = i+1
            data.append([compound, formula, adduct, round(lightMz, ndigits = 2), charge, MSXID])
            data.append([compound, formula, adduct, round(heavyMz, ndigits = 2), charge, MSXID])

        finalDf = pd.DataFrame(data, columns=['Compound','Formula','Adduct','m.z','z','MSXID'])
        finalDf.to_csv(header+str(CV)+'.txt', sep='\t', index=False)
    pass


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
Function:
Purpose:
Parameters:
Returns:
'''
import pickle
def quantify(inFile, libFile, mzxmlFiles, var):
    '''
    print(1)
################ TEMP ################
    fragDf = pd.read_csv('Data/Input/TempHold/mostintense_quantmzlist.txt', sep='\t')
#    libDict = mgf_library_upload_quant_delete('Data/Input/human.faims.fixed.decoy.mgf')

#    pickle.dump(libDict, open( "Data/Input/TempHold/lib.p", "wb" ))
#    libDict = pickle.load(open( "Data/Input/TempHold/lib.p", "rb" ))

################ TEMP ################
    print(5)

#    fragDf = pd.read_csv(inFile)
#    libDf = pd.read_csv(inFile)
    tmp_prec_dict = {}
    with mzxml.read(mzxmlFiles[0]) as spectra:
        for x in spectra:
            tmp_prec_dict[round(x['precursorMz'][0]['precursorMz'], 2),
                          round(x['precursorMz'][1]['precursorMz'], 2),
                          x['compensationVoltage']] = x['num']
#    pickle.dump(tmp_prec_dict, open( "Data/Input/TempHold/dict.p", "wb" ))

    print(10)

    peakDict = {}
    fragDict = {}
    tempLibDict = {}
    for i in range(len(fragDf)):
        #seq, mz, z, CV = fragDf['peptide'], fragDf['MzLIB'], fragDf['zLIB'], fragDf['CompensationVoltage']
        seq, mz, z, CV = fragDf.loc[i]['Peptide'], fragDf.loc[i]['prec_light_mz'], fragDf.loc[i]['z'], fragDf.loc[i]['CV']
        heavyMz = round(calc_heavy_mz(seq, mz, z), 2)
        lightMz = round(mz, 2)
        key = (lightMz, heavyMz, CV)
        scan = tmp_prec_dict[key]
        tempLibDict[mz, seq] = scan
        fragDict[scan] = {'seq':seq, 'mz':mz, 'z':z, 'CV':CV}

    print(15)

    uniquePeps = set(fragDf['Peptide'])
    digPat = r'\+\d+\.\d+'
    uniqueDigs = set()
    for pep in uniquePeps: uniqueDigs.update(re.findall(digPat, pep))
    digDict = dict(zip(uniqueDigs, [str(x) for x in list(range(len(uniqueDigs)))]))
    #NOTE: If you have more than 10 custom entries you'll have problems
    customAAcomp = dict(mass.std_aa_mass)
    for key in digDict: customAAcomp[digDict[key]] = float(key[1:]); print(customAAcomp[digDict[key]])

    count = 0
    for pep in uniquePeps: count += pep.count('C')
    libDict = mgf_library_upload_quant_delete('Data/Input/human.faims.fixed.decoy.mgf', tempLibDict, digDict, customAAcomp)

    pickle.dump(libDict, open( "Data/Input/TempHold/lib_3.p", "wb" ))
    pickle.dump(fragDict, open( "Data/Input/TempHold/frag.p", "wb" ))
    '''
    libDict = pickle.load(open( "Data/Input/TempHold/lib_3.p", "rb" ))
    fragDict = pickle.load(open( "Data/Input/TempHold/frag.p", "rb" ))

    minMatch = var[0]
    m = var[1]
    stdev = var[2]


    finalDf = pd.DataFrame(index=sorted(libDict.keys()))
    finalDf['peptide'] = [fragDict[key]['seq'] for key in sorted(fragDict.keys())]
    check = False
    for f in sorted(mzxmlFiles):
        #print(f)
        if f == 'Data/Input/jesse/20190503_DI2A_tMS2_OTmostint_A549_1to1_01.mzXML': check = True
        ppmDiffs = []
        with mzxml.read(f, use_index =True) as file:
            for scan in sorted(libDict.keys()):

                spec = file.get_by_id(scan)
                spec['intensity array'] = [x**0.5 for x in spec['intensity array']]
                peakIDs = [scan for x in range(len(spec['m/z array']))]
                expSpectrum = list(tuple(zip(spec['m/z array'],spec['intensity array'],peakIDs)))
                expSpectrum.sort(key=lambda x:x[0])
                l, h = spectra_peak_comparison(libDict[scan], expSpectrum, 50, 0, tag='qPPM')

                if len(l) > minMatch and len(h) > minMatch:
                    ppmDiffs += l; ppmDiffs += h


        offset, tolerance = find_offset_tol(sorted(ppmDiffs), 'Data/oldOutput/hist.png', stdev=stdev, mean=False)
        data = []

        ratioDict = {}
        with mzxml.read(f, use_index =True) as file:
            for scan in sorted(libDict.keys()):
                spec = file.get_by_id(scan)
                spec['intensity array'] = [x**0.5 for x in spec['intensity array']]
                peakIDs = [scan for x in range(len(spec['m/z array']))]
                expSpectrum = list(tuple(zip(spec['m/z array'],spec['intensity array'],peakIDs)))
                expSpectrum.sort(key=lambda x:x[0])
                l, h = spectra_peak_comparison(libDict[scan], expSpectrum, tolerance, offset, tag='ratio')
                lightSet = set(l.keys())
                heavySet = set(h.keys())


                smallPeakIntensity = np.mean(sorted(spec['intensity array'])[:10])
                onlyLight = lightSet - heavySet
                for key in onlyLight: h[key] = smallPeakIntensity
                onlyHeavy = heavySet - lightSet
                for key in onlyHeavy: l[key] = smallPeakIntensity

                keys = l.keys()
                if len(keys) > minMatch:
                    minKey = min(keys)
                    light = [l[key] for key in keys]
                    heavy = [h[key] for key in keys]
                    allInt = sorted(light + heavy)
                    mini = allInt[0] + allInt[1]
                    maxi = allInt[-1] + allInt[-2]



#                    mini = min(light+heavy)
#                    maxi = max(light+heavy)
                    ratios = [np.log2(x/y)*weight_ratio(x,y,mini,maxi) for x,y in zip(heavy, light)]
                    if m == 'mean': ratio = np.mean(ratios)
                    elif m=='median': ratio = np.median(ratios)
                    elif m=='intensity': ratio = np.log2((h[minKey]/l[minKey]))

                else:
                    ratio=np.nan

                '''
                print(l)
                print(h)
                print('-'*20)
                print(' '*20)
                '''
                ratioDict[scan] = ratio
        finalDf[f] = [ratioDict[key] for key in sorted(ratioDict.keys())]
    if var[2]==0.5: var[2]='P5'
    var = [str(x) for x in var]

    #finalDf.to_csv('Data/oldOutput/test_'+'_'.join(var)+'.csv')
    return finalDf

def weight_ratio(x, y, mini, maxi, scale=10):
    return 1
    #top = max([x,y])
    #top = np.mean([x,y])
    #top = x + y
    #return ((top-mini)/maxi)*scale

#        print()
'''
                    if check:
                        print('Scan number: ' + str(scan))
                        print('fragment info: ' + str(fragDict[scan]))
                        print('matched light: ' + str(sorted(l,reverse=True)))
                        print('matched heavy: ' + str(sorted(h,reverse=True)))
                        lightMz = [x[0] for x in libDict[scan] if x[2][0]=='light']
                        heavyMz = [x[0] for x in libDict[scan] if x[2][0]=='heavy']
                        lightInt = [x[1]**0.5 for x in libDict[scan] if x[2][0]=='light']
                        heavyInt = [x[1]**0.5 for x in libDict[scan] if x[2][0]=='heavy']
                        expMz = [x[0] for x in expSpectrum]
                        expInt = [x[1] for x in expSpectrum]

                        print('light mz: ' + str([(i,lightMz[i]) for i in range(len(lightMz))]))
                        print('heavy mz: ' + str([(i,heavyMz[i]) for i in range(len(heavyMz))]))
                        printI = [(i,lightInt[i]) for i in range(len(lightMz))]
                        printI.sort(key=lambda x:x[1], reverse=True)

                        print('intensity: ' + str(printI))


                        light = pd.DataFrame([lightInt], columns=lightMz)
                        heavy = pd.DataFrame([heavyInt], columns=heavyMz)
                        exp = pd.DataFrame([expInt], columns=expMz)

                        plot_spec(exp, 'black')
                        plot_spec(-light, 'cyan')
                        plot_spec(-heavy, 'blue')
                        pyplot.savefig('Data/oldOutput/check.png')
                        check = False
'''

def mgf_library_upload_quant_delete(fileName, d, digDict, aaDict):

    # mgf file is read in using the pyteomics mgf module
    libMGF = mgf.read(fileName)

    # Print statement for timing the program
    print('#Enter library dictionary upload: ')
    print('#'+str(timedelta(seconds=timer())))

    # return value is initialized
    lib = {}
    count=0
    time = timer()
    prevtime = time
    minPeaks = 1

    def repl(m):
        return digDict[m]
    # each spectrum in the mgf file
    for spec in libMGF:
        count += 1
        if count % 5000 == 0:
            time = timer()
            if len(peaks) > 0: print(str(count)+', '+str(time-prevtime) + ', ' + str(peaks[0]))
            prevtime = time
        seq = spec['params']['seq']
        precMz = spec['params']['pepmass'][0]
        key = (precMz, seq)
        if key not in d: continue
        sequence = re.sub(r'\+\d+\.\d+', lambda m: digDict.get(m.group()), seq)

        # peaks of the library file are intialized.
        mz = list(spec['m/z array'])
        intensity = list(spec['intensity array'])
        z = spec['params']['charge'][0]

        fragList = []
        for x in range(1, len(sequence)-1):
            fragseq = sequence[x:]
            lightfragmz = mass.fast_mass(sequence=sequence[x:], ion_type='y', charge=1, aa_mass = aaDict) # Do I need to use different possible charges?
            i = mz.index(round(lightfragmz, 3)) if round(lightfragmz, 3) in mz else -1
            if i==-1: continue
            fragList.append((intensity[i], lightfragmz, fragseq))


        fragList.sort(reverse=True)
        if len(fragList) >= minPeaks: fragList = fragList[:minPeaks]
        peaks = []
        for i in range(len(fragList)):
            fragMz = fragList[i][1]
            fragInt = fragList[i][0]
            peaks.append((fragMz, fragInt, ('light',i)))
            peaks.append((calc_heavy_mz(fragList[i][2], fragMz, 1), fragInt, ('heavy', i)))

        '''

        #sorted, then filter, then calculate heavy
        print(fragList)

        heavyMzs = []
        for i in range(1,z):
            heavyMzs.append([round(calc_heavy_mz(sequence, lightfragmz, i),2) for x in mz])

        #print(mz)
        #for x in fragList: print(x)




        mz = [round(x,2) for x in mz]
        id = list(range(len(mz)))
        lightTag = list(tuple(zip(['light' for x in mz],id)))
        heavyTag = list(tuple(zip(['heavy' for x in mz],id)))

        peaks = list(tuple(zip(mz,intensity,lightTag)))
        peaks.sort(key=lambda x:x[1])
        peaks = peaks[-minPeaks:]
        for heavyMz in heavyMzs:
            heavyPeak = list(tuple(zip(heavyMz,intensity,heavyTag)))
            heavyPeak.sort(key=lambda x:x[1])
            heavyPeak = heavyPeak[-minPeaks:]
            peaks += heavyPeak
        '''
        peaks.sort(key=lambda x:x[0])

        # entry placed in final dictionary
        lib[d[key]] = peaks


    return lib

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

def plot_spec(SPECTRA, COLOR):
    pyplot.vlines(SPECTRA.columns, np.repeat(0, len(SPECTRA.columns)), SPECTRA, colors=COLOR)
#    pyplot.ylim(-8000,8000)
#    ax = pyplot.gca()
#    ax.axes.xaxis.set_visible(False)
#    ax.axes.yaxis.set_visible(False)
