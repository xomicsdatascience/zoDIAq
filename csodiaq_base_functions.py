import pandas as pd
from pyteomics import mzxml
from bisect import bisect
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import mgf
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
Function: traml_library_upload_csv()
Purpose: Given a traml file in .tsv or .csv format representing the library spectra, this function generates
            a dictionary to be used for future analysis. key:value pairings are described in the "Returns"
            section of this comment. Note that the tuple key format was chosen because library spectra will
            need to be sorted and chosen by precursor mass, so identifiers already included in the library
            spectra were explicitly not used as keys.
Parameters:
    'fileName' - string representing the path to the library spectra file (required traml .tsv or .csv format).
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
Function: lib_mz_match_query_window()
Purpose: This function determines which library spectra fit within the baseline precursor mass window of the
            experimental spectrum. This significantly reduces the time complexity of the overall analysis by eliminating
            all irrelevant library spectra from consideration, the vast majority of them.
Parameters:
    'spec' - dictionary corresponding to pyteomics.mzxml.read() values. This variable contains all data corresponding
        to the query spectrum in question.
    'sortedLibKeys' - a list of (float, string) tuples. List represents all keys of the dictionary 'lib' from
        traml_library_upload_csv(). See dictionary 'lib' key explanation in function traml_library_upload_csv().
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
    'lib' - dictionary as returned by the traml_library_upload_csv() function.
    'libKeys' - list of keys corresponding to the 'lib' parameter as returned by the lib_mz_match_query_window()
        function.
Returns:
    A list of (float, float, tuple) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        traml_library_upload_csv().
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
        traml_library_upload_csv().
    'libKeys' - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
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
        key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
        value - list of floats
            float - sum of products of intensities between matched library and query peaks. In context of the cosine
                similarity score algorithm, this is the "A*B" sum.
            float - sum of squared query peak intensities when matched to a library peaks. In context of the cosine
                similarity score algorithm, this is the "B^2" sum.
            float - sum of squared library peak intensities when matched to a query peaks. In context of the cosine
                similarity score algorithm, this is the "A^2" sum.
    dictionary 'countDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
        value - int
            Representing the number of matched peaks. This value will be the 'shared' column value in the output file.
    dictionary 'ionDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
        value - set representing the indices of the query peaks that matched the library spectrum represented in this
                dictionary value. The sum of these query peak intensities will be returned as the 'ionCount' column
                value in the output file.
    dictionary 'ppmDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
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
    'lib' - dictionary as returned by the traml_library_upload_csv() function.
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
def fdr_calculation(df, FDRCutoff, fdrList=False):
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
        if curFDR > FDRCutoff:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if len(fdrValues) < 1/FDRCutoff:
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
    'FDRCutoff' - float representing the minumum allowed FDR of the csodiaq output.
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

def find_offset_tol(data, histFile):

    hist, bins = np.histogram(data, bins=200)

    # Threshold frequency
    freq = sum(hist)/len(hist)

    hist2 = [x for x in hist]
    # Zero out low values
    hist[np.where(hist <= freq)] = 0

    # Plot
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    index_max = max(range(len(hist)), key=hist.__getitem__)

    min_i = 0
    max_i = 0
    for i in range(len(hist)):
        if hist[i] != 0: min_i = i; break
    for i in range(min_i, len(hist)):
        if hist[i] == 0: max_i = i-1; break

    offset2 = center[index_max]
    tolerance2 = (center[max_i] - center[min_i])/2

    offset = sum(data)/len(data)
    tolerance = statistics.pstdev(data)*2
    if histFile:
        plt.clf()
        plt.bar(center, hist2, align='center', width=width)
        plt.title("Gaussian Histogram")
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.axvline(x=offset, color='black', linestyle = 'dashed')
        plt.axvline(x=offset-tolerance, color='red', linestyle = 'dashed')
        plt.axvline(x=offset+tolerance, color='red', linestyle = 'dashed')
        plt.savefig(histFile)
    return offset, tolerance

def write_ppm_spread_decoy(cosFile, ppmFile, outFile):
    # Data is read in.
    print(1)
    df = pd.read_csv(cosFile)
    print(2)
    ppmDict = read_ppm_file_to_dict(ppmFile)
    print(3)
    # list of keys corresponding to ppmDict are generated from the csodiaq data frame.
    listOfKeys = [(df['scan'].loc[i],df['peptide'].loc[i],df['protein'].loc[i]) for i in range(len(df))]

    # all values from the ppmFile corresponding to those keys are saved into a single list.
    ppmList = []
    decoyList = []
    for key in listOfKeys:
        if 'DECOY' not in key[2]:
            ppmList += ppmDict[key]
#            ppmList.append(sum([float(x) for x in ppmDict[key]])/len(ppmDict[key]))

    print(4)
    for key in listOfKeys:
        if 'DECOY' in key[2]:
            decoyList += ppmDict[key]
#            decoyList.append(sum([float(x) for x in ppmDict[key]])/len(ppmDict[key]))
    print(5)

    # CSV file is written with the ppm spread. It's essentially one long row of ppm difference values.
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(ppmList)
        writer.writerow(decoyList)

def mgf_library_upload(fileName, numLibPeaks):

    libMGF = mgf.read(fileName)
    print('#Enter library dictionary upload: ')
    print('#'+str(timedelta(seconds=timer())))
    libDict = {}
    for spec in libMGF:
        key = (spec['params']['pepmass'][0], spec['params']['seq'])
        charge = spec['params']['charge'][0]
        name = spec['params']['title']
        if 'protein' in spec['params']: protein = spec['params']['protein']
        else: protein = ''
        mz = spec['m/z array']
        intensity = spec['intensity array']
        intensity = [x**0.5 for x in intensity]
        keyList = [key for x in mz]
        peaks = list(tuple(zip(mz,intensity,keyList)))
        peaks.sort(key=lambda x:x[1],reverse=True)
        peaks = peaks[:numLibPeaks]
        peaks.sort(key=lambda x:x[0])
        tempDict = {
            'PrecursorCharge':charge,
            'transition_group_id':name,
            'ProteinName':protein,
            'Peaks':peaks
        }
        libDict[key] = tempDict
    return libDict

def library_file_to_dict(inFile, numLibPeaks):
    fileType = inFile.split('.')[-1]
    if fileType == 'mgf':
        lib = mgf_library_upload(inFile, numLibPeaks)
    else:
        lib = traml_library_upload_csv(inFile, numLibPeaks)
    return lib

def write_csodiaq_fdr_outputs(inFile, specFile, pepFile, protFile):
    overallDf = pd.read_csv(inFile).sort_values('cosine', ascending=False).reset_index(drop=True)

    count = 0


    print(1)
    spectralDf = add_fdr_to_csodiaq_output(overallDf)

    print(2)
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    print(3)
    # getting valid proteins
    peptideProteinConnections = []
    for i in range(len(peptideDf)):
        peptide = peptideDf['peptide'].loc[i]
        proteinGroup = peptideDf['protein'].loc[i]
        proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
        for pro in proteins:
            protein = pro[0] + pro[1]
            peptideProteinConnections.append((peptide,protein))
    verifiedProteinList = idp.find_valid_proteins(peptideProteinConnections)

    print(len(verifiedProteinList))


    print(4)
    # for each protein in the verified list, add it to a new dataframe with the protein in the "leading protein" column
    tempProtDf = add_leading_protein_column(peptideDf, verifiedProteinList)
    tempProtDf.to_csv(protFile)

    print(5)
    proteinDf = add_fdr_to_csodiaq_output(tempProtDf, filterType='leadingProtein', DEBUG=True)

    print(6)
    proteinBestNum = min(set(proteinDf['shared']))
    tempPeptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide', bestMatchNum=proteinBestNum)

    proteinDict = proteinDf.set_index('leadingProtein').T.to_dict()
    proteinDf = tempPeptideDf.copy()
    proteinDf = add_leading_protein_column(proteinDf, verifiedProteinList)
    print(7)
    proteinDf.to_csv(protFile)

    proteinCosine = []
    proteinFDR = []
    removables = []
    for i in range(len(proteinDf)):
        protein = proteinDf['leadingProtein'].loc[i]
        print('--------------')
        print('protein: ' + str(protein))
        print('index: '+str(i))

        if protein in proteinDict:
            proteinCosine.append(proteinDict[protein]['cosine'])
            proteinFDR.append(proteinDict[protein]['leadingProteinFDR'])
            print('cosine: '+str(proteinDict[protein]['cosine']))
            print('leadingProteinFDR: '+str(proteinDict[protein]['leadingProteinFDR']))
        else:
            removables.append(i)

    proteinDf = proteinDf.drop(proteinDf.index[removables]).reset_index(drop=True)
    proteinDf['proteinCosine'] = proteinCosine
    proteinDf['leadingProteinFDR'] = proteinFDR
    print(8)

    proteinDf = proteinDf.sort_values(['proteinCosine','protein', 'cosine'], ascending=[False,False,False]).reset_index(drop=True)
    tempDf = proteinDf.drop_duplicates(subset='protein', keep='first').reset_index(drop=True)

    spectralDf.to_csv(specFile)
    peptideDf.to_csv(pepFile)
    proteinDf.to_csv(protFile)

def add_leading_protein_column(df, verifiedProteins):
    tempDf = pd.DataFrame(columns = df.columns)
    leadingProteins = []
    for i in range(len(df)):
        proteinGroup = df['protein'].loc[i]
        proteins = re.findall('(DECOY_0_)?(sp\|\w{6}\|)', proteinGroup)
        for pro in proteins:
            protein = pro[0] + pro[1]
            if protein in verifiedProteins:
                tempDf = tempDf.append(df.loc[i])
                leadingProteins.append(protein)
    tempDf['leadingProtein'] = leadingProteins
    tempDf = tempDf.reset_index(drop=True)
    return tempDf

def add_fdr_to_csodiaq_output(df, filterType='spectral', bestMatchNum=0, DEBUG=True):

    if not bestMatchNum:
        matches = sorted(list(set(df['shared'])))
        matchNumHits = []
        for m in matches:
            tempDf = df[df['shared'] >= m].reset_index(drop=True)
            if filterType != 'spectral': tempDf = tempDf.drop_duplicates(subset=filterType, keep='first').reset_index(drop=True)
            fdrRows, dec = fdr_calculation(tempDf, 0.01)
            matchNumHits.append(fdrRows)
        bestMatchNum = matches[matchNumHits.index(max(matchNumHits))]

    tempDf = df[df['shared'] >= bestMatchNum].reset_index(drop=True)
    if filterType != 'spectral': tempDf = tempDf.drop_duplicates(subset=filterType, keep='first').reset_index(drop=True)

    fdrList, decoyNum = fdr_calculation(tempDf, 0.01, fdrList=True)
    tempDf = tempDf.truncate(after=len(fdrList)-1)
    tempDf[filterType + 'FDR'] = fdrList
    tempDf = tempDf.reset_index(drop=True)

    return tempDf


'''
    for x in verifiedProteinList:
        if 'DECOY' not in x:
            count += 1
            if count < 10: print(x)
        if 'DECOY' in x:
            decoyCount += 1
            if decoyCount < 10: print(x)
        if decoyCount > 10 and count > 10: break
'''
