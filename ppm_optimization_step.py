import pandas as pd
from pyteomics import mzxml
from bisect import bisect
from timeit import default_timer as timer
from datetime import timedelta
import csv
import statistics
import matplotlib.pyplot as plt

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
    fileName - string representing the path to the library spectra file (required traml .tsv or .csv format).
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

    # TEMP - Normalize intensities by finding their square root
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
        lib[key]['Peaks'] = list(tuple(zip(mz,intensity,keyList)))
    return lib


'''
Function: lib_mz_match_query_window()
Purpose: This function determines which library spectra fit within the baseline precursor mass window of the
            experimental spectrum. This significantly reduces the time complexity of the overall analysis by eliminating
            all irrelevant library spectra from consideration, the vast majority of them.
Parameters:
    spec - dictionary corresponding to pyteomics.mzxml.read() values. This variable contains all data corresponding
        to the query spectrum in question.
    sortedLibKeys - a list of (float, string) tuples. List represents all keys of the dictionary 'lib' from
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
    lib - dictionary as returned by the traml_library_upload_csv() function.
    libKeys - list of keys corresponding to the 'lib' parameter as returned by the lib_mz_match_query_window()
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
            If the difference is outside the provided tolerance, 0 is returned instead.
Parameters:
    x - first float value, corresponding to the theoretical mass (library mz).
    y - second float value, corresponding to the experimental mass (query mz).
    ppm - ppm used to determine tolerance.
Returns:
    boolean - True indicates the two values are approximately equal, False otherwise.
'''
def approx(x, y, ppmTol):
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)


'''
Function: peak_comparison_for_cosine()
Purpose: This function determines peaks that are within a sufficiently close tolerance to one another as to be included
            in the cosine similarity score calculation. Specifically, this function iterates over peaks of the
            pooled library spectrum and the experimental spectrum in an O(n+m) manner. The cosine score is not explicitly
            calculated in this function - that occurs in the cosine_similarity() function. Rather, this function
            provides the variables necessary for the cosine similarity score to be calculated.
Parameters:
    libSpectrum - list of (float, float, tuple) tuples. Represents spectrum peaks - see 'Peaks' key explanation in function
        traml_library_upload_csv().
    libKeys - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
    expMz - list of floats. Corresponds to m/z values from the experimental spectrum being analyzed. Matched by index
        with the 'expIntensity' parameter.
    expIntensity -list of floats. Corresponds to intensity values from the experimental spectrum being analyzed. Matched
        by index with the 'expMz' parameter.
    ppm - see 'ppm' parameter description for function approx().
Returns:
    dictionary 'libDict'
        key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
        value - list of various values.
            float - sum of products of intensities between matched library and query peaks. In context of the cosine
                similarity score algorithm, this is the "A*B" sum.
            float - sum of squared query peak intensities when matched to a library peaks. In context of the cosine
                similarity score algorithm, this is the "B^2" sum.
            float - sum of squared library peak intensities when matched to a query peaks. In context of the cosine
                similarity score algorithm, this is the "A^2" sum.
            int - representing the number of matched peaks. This value will be the 'shared' column value in the output file.
            set - Representing the indices of the query peaks that matched the library spectrum represented in this
                dictionary value. The sum of these query peak intensities will be returned as the 'ionCount' column
                value in the output file.
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
Purpose: Provided the output of the peak_comparison_for_cosine() function, this function explicitly calculates the
            cosine similarity score and returns it. To save on time, this function also takes in output from the
            peak_comparison_for_cosine() corresponding to the number of matched peaks and the ion count and returns
            them with the cosine score.
Parameters:
    key - (float, string) tuple. See dictionary 'lib' key explanation in function traml_library_upload_csv().
    row - list of various values. See dictionary 'libDict' value explanation in function peak_comparison_for_cosine().
    expIntensity - list of floats. See 'expIntensity' parameter explanation in function peak_comparison_for_cosine().
Returns:
    a list of various values.
        (float, string) tuple - see dictionary 'lib' key explanation in function traml_library_upload_csv().
        float - corresponds to the final cosine similarity score between library spectrum indicated by the 'key'
            paramter/first return value and the experimental spectrum being analyzed.
        int - representing the number of matched peaks. This value will be the 'shared' column value in the output file.
        float - representing the sum of the matched query peak intensities. This value will be returned as the
        'ionCount' column value in the output file.
'''
def cosine_similarity(row):
    magnitude = (row[1]**0.5) * (row[2]**0.5) # (sqrt(sum(A^2))*sqrt(sum(B^2)))
    return (row[0] / magnitude if magnitude else 0)

'''
Function: query_spectra_analysis()
Purpose: This function loops through all query spectra and calculates the cosine similarity score between
            it and every library spectra with a precursor m/z value within its designated window.
Parameters:
    expSpectraFile - string representing the path to the query spectra file (required .mzXML format).
    outFile - string representing the path to the output/results file.
    lib - dictionary as returned by the traml_library_upload_csv() function.
Returns:
    No Return Value. Results are written directly to the output file. The description of specific columns of output
        file are provided in the function comments.
'''
def query_spectra_analysis( expSpectraFile, outFile, ppmFile, lib, numPeakMatch, ppmTol, ppmYOffset):
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
        'cosine', # cosine score comparing the library spectrum corresponding to this row with the query spectrum.
        'name', # Title - corresponds to the column "transition_group_id," a library spectrum identifier.
        'Peak(Query)', # The number of peaks in the query spectrum.
        'Peaks(Library)', # The number of peaks in the library spectrum.
        'shared', # The number of peaks that matched between query spectrum/library spectrum.
        'ionCount', # Sum of query spectrum intensities, excluding possible duplicates - currently uncalculated, set to 0.
        'CompensationVoltage', #the compensation voltage of the query spectrum.
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
                    if countDict[key] > numPeakMatch:
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
                        ppmWriter.writerow([spec['num'], key[1], lib[key]['ProteinName']] + ppmDict[key])

        # Prints the final number of experimental spectra analyzed.
        print('#'+str(count))

def fdr_calculation(df, FDRCutoff=0.01):
    fdr = []
    numDecoys = 0
    for i in range(len(df)):
        if 'DECOY' in df.loc[i]['protein']:
            numDecoys += 1
        newFDR = numDecoys/(i+1)
        if newFDR > FDRCutoff:
            if len(fdr) < 1/FDRCutoff:
                return [], 0
            return fdr, numDecoys-1
        fdr.append(newFDR)
    return fdr, numDecoys-1

def find_best_matchNum_fdr(df, FDRCutoff):
    matches = sorted(list(set(df['shared'])))
    bestMatchNum = 0
    bestFDR = 0
    for m in matches:
        tempDf = df[df['shared'] >= m].reset_index(drop=True)
        fdr, decoys = fdr_calculation(tempDf, FDRCutoff=FDRCutoff)
        if len(fdr) > bestFDR:
            bestMatchNum = m
            bestFDR = len(fdr)
    return bestMatchNum, bestFDR

def read_ppm_file_to_dict(ppmFile):
    ppmDict = {}
    with open(ppmFile, 'r') as r:
        csvReader = csv.reader(r)
        for row in csvReader:
            ppmDict[(int(row[0]),row[1],row[2])] = row[3:]
    return ppmDict

def write_ppm_spread(cosFile, ppmFile, outFile):
    df = pd.read_csv(cosFile)
    ppmDict = read_ppm_file_to_dict(ppmFile)
    listOfKeys = [(df['scan'].loc[i],df['peptide'].loc[i],df['protein'].loc[i]) for i in range(len(df))]
    ppmList = []
    for key in listOfKeys: ppmList += ppmDict[key]
    with open(outFile, 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(ppmList)

def ppm_offset(mz, ppm):
    return (mz*ppm)/1000000
