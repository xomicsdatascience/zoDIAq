import pandas as pd
from pyteomics import mzxml, mgf
from bisect import bisect
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as plt
import numpy as np
import idpicker as idp
import re
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# Template for function descriptions
'''
Function:
Purpose:
Parameters:
Returns:
'''

'''
Function: library_file_to_dict()
Purpose: This function reads the file type of the input file (accepts .mgf, .csv and .tsv formats) and runs it through the
            appropriate reading function. Basically, this function is a catch-all used in the corresponding menu function to
            determine how to best read in a particular library spectra file.
Parameters:
    'inFile' - string representing the file path to the desired library spectra file.
    'numLibPeaks' - int representing the maximum number of peaks allowed per library spectrum. The [numLibPeaks] peaks with
        the highest intensity are included, everything else is excluded. Default maximum number of peaks is 31, as determined
        in the corresponding menu function.
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
def library_file_to_dict(inFile, numLibPeaks):
    fileType = inFile.split('.')[-1]
    if fileType == 'mgf':
        lib = mgf_library_upload(inFile, numLibPeaks)
    else:
        lib = traml_library_upload_csv(inFile, numLibPeaks)
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
def traml_library_upload_csv(fileName, numLibPeaks):
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
        peaks.sort(key=lambda x:x[1],reverse=True)
        peaks = peaks[:numLibPeaks]
        peaks.sort(key=lambda x:x[0])
        lib[key]['Peaks'] = peaks
    return lib


'''
Function: mgf_library_upload()
Purpose: Same purpose as traml_library_upload_csv() function, but for mgf files.
Parameters: see traml_library_upload_csv() parameters.
Returns: see traml_library_upload_csv() return values.
'''
def mgf_library_upload(fileName, numLibPeaks):

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

        # low-intensity peaks are filtered out, based on value of 'numLibPeaks'
        peaks.sort(key=lambda x:x[1],reverse=True)
        peaks = peaks[:numLibPeaks]
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
def lib_mz_match_query_window( spec, sortedLibKeys ):
    # Values will be added to sortedLibKeys parameter, so a copy of that list is created here.
    temp = sortedLibKeys[:]

    # A fake 'key' value is created representing the upper limit of keys that correspond to the query spectrum.
    top_mz = spec['precursorMz'][0]['precursorMz'] + spec['precursorMz'][0]['windowWideness'] / 2
    top_key = (top_mz, "z")

    # A fake 'key' value is created representing the lower limit of keys that correspond to the query spectrum.
    bottom_mz = spec['precursorMz'][0]['precursorMz'] - spec['precursorMz'][0]['windowWideness'] / 2
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
Function: peak_comparison()
Purpose: This function determines peaks that are within a sufficiently close tolerance to one another. Various functions
            requiring peak comparison then take place inside this function, including compiling data for a cosine
            similarity score. The specifics of the various functions can be found in the return values below.
Parameters:
    'libSpectrum' - list of (float, float, tuple) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        library_file_to_dict().
    'libKeys' - (float, string) tuple. See dictionary 'lib' key explanation in function library_file_to_dict().
    'expMz' - list of floats. Corresponds to m/z values from the experimental spectrum being analyzed. Matched by index
        with the 'expIntensity' parameter.
    'expIntensity' -list of floats. Corresponds to intensity values from the experimental spectrum being analyzed. Matched
        by index with the 'expMz' parameter.
    'ppmTol' - see 'ppmTol' parameter description for function approx().
    'ppmYOffset' - The value in ppm that should be subtracted from the query peak m/z. This is the average ppm difference
        of an uncorrected csodiaq output, being applied to the experimental spectra to generate a corrected
        csodiaq calculation output.
Returns:
    dictionary 'cosDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function library_file_to_dict().
        value - list of floats
            float - sum of products of intensities between matched library and query peaks. In context of the cosine
                similarity score algorithm, this is the "A*B" sum.
            float - sum of squared query peak intensities when matched to a library peaks. In context of the cosine
                similarity score algorithm, this is the "B^2" sum.
            float - sum of squared library peak intensities when matched to a query peaks. In context of the cosine
                similarity score algorithm, this is the "A^2" sum.
    dictionary 'countDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function library_file_to_dict().
        value - int
            Representing the number of matched peaks. This value will be the 'shared' column value in the output file.
    dictionary 'ionDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function library_file_to_dict().
        value - set representing the indices of the query peaks that matched the library spectrum represented in this
                dictionary value. The sum of these query peak intensities will be returned as the 'ionCount' column
                value in the output file.
    dictionary 'ppmDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function library_file_to_dict().
        value - list of floats represents the ppm difference of all peak matches that were within the given ppm tolerance. This is
                most directly used for determining the ppm offset and standard deviation of an uncorrected
                csodiaq output to be used in generating a corrected csodiaq output. For a corrected csodiaq output,
                it is primarily used for figure generation to compare with uncorrected csodiaq output.
'''
def peak_comparison(libSpectrum, libKeys, expMz, expIntensity, ppmTol, ppmYOffset):
    # final dictionary returned is initialized. See function description for details on contents.
    cosDict = {k: [0.0, 0.0, 0.0] for k in libKeys}
    countDict = {k: 0 for k in libKeys}
    ionDict = {k: set() for k in libKeys}
    ppmDict = {k: [] for k in libKeys}

    # By tracking the indices of the current library/query peaks we reduce the time complexity of the algorithm
    #   from O(n*m) to O(n+m) where n = (# of library peaks) and m = (# of query peaks).
    i, j = 0, 0
    expPeakMz = expMz[j] - ppm_offset(expMz[j], ppmYOffset)
    while i < len(libSpectrum) and j < len(expMz):

        # If the m/z of the peaks are not within the given ppm tolerance, the indices of the smaller of the two is incremented
        #   and we are returned to the top of the while loop.
        if not approx(libSpectrum[i][0],expPeakMz, ppmTol):
            if libSpectrum[i][0] > expPeakMz:
                j += 1
                if j < len(expMz): expPeakMz = expMz[j] - ppm_offset(expMz[j], ppmYOffset)
                continue
            if libSpectrum[i][0] < expPeakMz: i += 1; continue

        # To account for the matching of one query peak to multiple library peaks, library peaks are looped over
        #   after the initial match. Every matching peak contributes to the various variables in the final returned dictionary.
        p = i + 0
        while (p < len(libSpectrum)):
            ppm = approx(libSpectrum[p][0], expPeakMz, ppmTol)
            if p==len(libSpectrum) or not ppm: break

            # Data required to calculate the cosine score
            cosDict[libSpectrum[p][2]][0] += libSpectrum[p][1]*expIntensity[j]
            cosDict[libSpectrum[p][2]][1] += expIntensity[j]**2
            cosDict[libSpectrum[p][2]][2] += libSpectrum[p][1]**2

            # Number of matching peaks
            countDict[libSpectrum[p][2]] += 1

            # Ion Count
            ionDict[libSpectrum[p][2]].add(j)

            # ppm differences in matched peaks
            ppmDict[libSpectrum[p][2]].append(ppm)
            p += 1

        # Note that the possibility of one library peak matching to multiple query peaks is automatically accounted for
        #   by the fact that the query peak is the next default increment after all match calculations have been made.
        j += 1
        if j < len(expMz): expPeakMz = expMz[j] - ppm_offset(expMz[j], ppmYOffset)
    return cosDict, countDict, ionDict, ppmDict


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
Function: query_spectra_analysis()
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
def query_spectra_analysis( expSpectraFile, outFile, ppmFile, lib, ppmTol, ppmYOffset):
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
        'totalWindowWidth' # width of m/z that was captured in the query spectrum. Corresponds to MzEXP.
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

        # Time taken to analyze every 100 spectra is recorded and printed to the screen in a csv format (can be copied and pasted).
        #   Printing was chosen over direct writing for easy tracking of program progress.
        time = timer()
        prevtime = time
        print("Index,TimeElapsed,NumPeaks,ExpPrecursorMz,NumLibrarySpectra")

        # Beginning to loop over query spectra.
        with mzxml.read(expSpectraFile) as spectra:
            for spec in spectra:

                # TEMP - Normalize intensities by finding their square root
                spec['intensity array'] = [x**0.5 for x in spec['intensity array']]

                # Printing time taken to analyze every 100 spectra in csv format.
                count += 1
                if count % 100 == 0:
                    time = timer()
                    print(str(count)+','+str(time-prevtime)+','+str(len(spec['m/z array']))+','+str(spec['precursorMz'][0]['precursorMz'])+','+str(len(libKeys)))
                    prevtime = time

                # Only library spectra with a precursor mass that falls within the target window of the experimental spectrum are included. See lib_mz_match_query_window() function description for more details.
                libKeys = lib_mz_match_query_window( spec, allLibKeys )

                # Cosine score is generated for each library spectrum and returned in a list. See peak_comparison_for_cosine() function description for formatting.
                if len(libKeys) != 0:
                    cosDict, countDict, ionDict, ppmDict = peak_comparison(pool_lib_spectra(lib, libKeys), libKeys, spec['m/z array'], spec['intensity array'], ppmTol, ppmYOffset)
                else: continue

                # Cosine scores and supplementary data are written to the output file.
                # Note that this loop essentially goes over every library spectrum identified by the lib_mz_match_query_window() function.
                for key in libKeys:
                    cosine = cosine_similarity(cosDict[key])
                    ionCount = sum([ spec['intensity array'][j]+ppm_offset(spec['intensity array'][j],ppmYOffset) for j in ionDict[key] ])
                    # Library spectra that had too few matching peaks are excluded. numPeakMatch variable determines the threshold.
                    if countDict[key] > 3:
                        temp = [
                            expSpectraFile, #fileName
                            spec['num'], #scan
                            spec['precursorMz'][0]['precursorMz'], #MzEXP
                            spec['precursorMz'][0]['precursorCharge'], #zEXP
                            key[1], #peptide
                            lib[key]['ProteinName'], #protein
                            key[0], #MzLIB
                            lib[key]['PrecursorCharge'], #zLIB
                            cosine, #cosine
                            lib[key]['transition_group_id'], #name
                            len(spec['m/z array']), #Peaks(Query)
                            len(lib[key]['Peaks']), #Peaks(Library)
                            countDict[key], #shared
                            ionCount, #ionCount
                            spec['compensationVoltage'], #compensationVoltage
                            spec['precursorMz'][0]['windowWideness'] #totalWindowWidth
                        ]
                        writer.writerow(temp)
                        ppmWriter.writerow([spec['num'], key[1], lib[key]['ProteinName']] + sorted(ppmDict[key]))

        # Prints the final number of experimental spectra analyzed.
        print('#'+str(count))

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
Function: find_best_matchNum_fdr()
Purpose: Given the output of the query_spectra_analysis() function as contained in a csodiaq output file and an
            FDR cutoff point, this function determines the optimal minimum allowed peak number of the output.
            It also provides the number of rows available after filtering for peak number and FDR cutoff,
            allowing for easy filtering of the pandas dataframe later.
Parameters:
    'df' - Pandas Data Frame representing the contents of a csodiaq output file.
Returns:
    'bestMatchNum' - int representing the optimal minimum allowed peak number.
    'bestFDR' - int representing the number of rows that are kept after filtering for both peak matches equal to
        and above bestMatchNum and values above the FDR cutoff.
'''
def find_best_matchNum_fdr(df, FDRCutoff):
    # list is created representing every unique number of matches represented in the dataframe
    matches = sorted(list(set(df['shared'])))
    bestMatchNum = 0
    bestFDR = 0

    # loops over every unique number of matches represented in the dataframe.
    for m in matches:

        # Temporary dataframe created, representing the overall data set filtered for the number of matches in question
        tempDf = df[df['shared'] >= m].reset_index(drop=True)

        # The number of rows kept after filtering for the FDR cutoff is calculated.
        fdrLength, decoys = fdr_calculation(tempDf, FDRCutoff=FDRCutoff)

        # The number of matches with the highest number of kept rows is considered the optimal number of matches allowed.
        if fdrLength > bestFDR:
            bestMatchNum = m
            bestFDR = fdrLength
    return bestMatchNum, bestFDR

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
def write_ppm_spread(cosFile, ppmFile, outFile):
    # Data is read in.
    df = pd.read_csv(cosFile)
    df.fillna('',inplace=True)

    ppmDict = read_ppm_file_to_dict(ppmFile)

    # list of keys corresponding to ppmDict are generated from the csodiaq data frame.
    listOfKeys = [(df['scan'].loc[i],df['peptide'].loc[i],df['protein'].loc[i]) for i in range(len(df))]

    # all values from the ppmFile corresponding to those keys are saved into a single list.
    ppmList = []
    for key in listOfKeys: ppmList += ppmDict[key]

    # CSV file is written with the ppm spread. It's essentially one long row of ppm difference values.
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(ppmList)


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
def find_offset_tol(data, histFile):

    # offset is calculated as the mean value of the provided data
    offset = sum(data)/len(data)

    # tolerance is calculated as the second standard deviation of the provided data
    tolerance = statistics.pstdev(data)*2

    # if a histogram file is provided, it is created, with offset (black) and tolerance (red) lines drawn for reference
    if histFile:
        hist, bins = np.histogram(data, bins=200)
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.clf()
        plt.bar(center, hist, align='center', width=width)
        plt.title("Gaussian Histogram of Matched Peak PPM Difference Spread")
        plt.xlabel("Difference between Matched Peaks (PPM)")
        plt.ylabel("Frequency")
        plt.axvline(x=offset, color='black', linestyle = 'dashed')
        plt.axvline(x=offset-tolerance, color='red', linestyle = 'dashed')
        plt.axvline(x=offset+tolerance, color='red', linestyle = 'dashed')
        plt.savefig(histFile)
    return offset, tolerance


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
    overallDf = pd.read_csv(inFile).sort_values('cosine', ascending=False).reset_index(drop=True)

    # spectral FDR is calculated and written to dataframe 'spectralDf'
    spectralDf = add_fdr_to_csodiaq_output(overallDf)

    # peptide FDR is calculated and written to dataframe 'peptideDf'
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    # proteins/protein groups considered present are determined using the IDPicker algorithm.
    #   To eliminate invalid peptides from consideration, only proteins found in the peptide FDR dataframe 'peptideDf' are considered.

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
    tempProtDf = add_leading_protein_column(peptideDf, verifiedProteinDict)

    # Protein FDR is calculated using the highest-scoring peptide for each protein group.
    proteinDf = add_fdr_to_csodiaq_output(tempProtDf, filterType='leadingProtein', DEBUG=True)

    # Because the optimal minimum number of allowed peak matches can differ between overall peptide and protein FDR calculations,
    #   the peptide FDR must be recalculated using the optimal protein minimum number for inclusion with the protein FDR output
    proteinBestNum = min(set(proteinDf['shared']))
    tempPeptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide', bestMatchNum=proteinBestNum)

    # The re-calculated peptide FDR is used as the template for the final protein FDR output.
    proteinDict = proteinDf.set_index('leadingProtein').T.to_dict()
    proteinDf = tempPeptideDf.copy()
    proteinDf = add_leading_protein_column(proteinDf, verifiedProteinDict)

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

    # Data from all of the above dataframes are written to their respective files.
    spectralDf.to_csv(specFile, index=False)
    peptideDf.to_csv(pepFile, index=False)
    proteinDf.to_csv(protFile, index=False)

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

    # For proteins that were part of an IDPicker-determined protein group the same row was added for every protein in the group. This is dropping those duplicate rows.
    finalDf = finalDf.drop_duplicates(keep='first').reset_index(drop=True)
    return finalDf


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

    # if no bestMatchNum value provided, the optimal minimum number of matches allowed is calculated and used.
    #   NOTE: this next part is essentially the find_best_matchNum_fdr() function, but also allows for filtering. I may want
    #   to add a new parameter to that function to reduce code duplication.
    if not bestMatchNum:
        matches = sorted(list(set(df['shared'])))
        matchNumHits = []
        for m in matches:
            finalDf = df[df['shared'] >= m].reset_index(drop=True)

            # if a filterType was provided, the dataframe is filtered by the corresponding filterType
            if filterType != 'spectral': finalDf = finalDf.drop_duplicates(subset=filterType, keep='first').reset_index(drop=True)
            fdrRows, dec = fdr_calculation(finalDf)
            matchNumHits.append(fdrRows)
        bestMatchNum = matches[matchNumHits.index(max(matchNumHits))]

    # final dataframe is initialized as a filtered version of the paraemter input 'df' - filters are minimum number of
    #   allowed matches and highest-scoring unique value of a column type if given
    finalDf = df[df['shared'] >= bestMatchNum].reset_index(drop=True)
    if filterType != 'spectral': finalDf = finalDf.drop_duplicates(subset=filterType, keep='first').reset_index(drop=True)

    # FDR is calculated
    fdrList, decoyNum = fdr_calculation(finalDf, fdrList=True)

    # All rows below the FDR cutoff are removed
    finalDf = finalDf.truncate(after=len(fdrList)-1)

    # FDR column is added to the dataframe
    finalDf[filterType + 'FDR'] = fdrList
    finalDf = finalDf.reset_index(drop=True)
    return finalDf
