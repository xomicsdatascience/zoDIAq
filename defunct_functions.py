

'''
DEFUNCT: This only pools library spectra, it does NOT pool query spectra
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
def pooled_library_query_spectra_analysis( expSpectraFile, outFile, ppmFile, lib, ppmTol, ppmYOffset):
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
        'csoDIAq_Score'
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

        # Beginning to loop over query spectra.
        with mzxml.read(expSpectraFile) as spectra:

            # Time taken to analyze every 100 spectra is recorded and printed to the screen in a csv format (can be copied and pasted).
            #   Printing was chosen over direct writing for easy tracking of program progress.
            time = timer()
            prevtime = time

            for spec in spectra:
#                if int(spec['num']) < 3601:

                # TEMP - Normalize intensities by finding their square root
                spec['intensity array'] = [x**0.5 for x in spec['intensity array']]

                # Only library spectra with a precursor mass that falls within the target window of the experimental spectrum are included. See lib_mz_match_query_window() function description for more details.
                top_mz = spec['precursorMz'][0]['precursorMz'] + spec['precursorMz'][0]['windowWideness'] / 2
                bottom_mz = spec['precursorMz'][0]['precursorMz'] - spec['precursorMz'][0]['windowWideness'] / 2
                libKeys = lib_mz_match_query_window( top_mz, bottom_mz, allLibKeys )

                # Printing time taken to analyze every 100 spectra.
                count += 1
                if count % 100 == 0:
                    time = timer()
                    print(str(count)+','+str(time-prevtime)+','+str(len(spec['m/z array']))+','+str(spec['precursorMz'][0]['precursorMz'])+','+str(len(libKeys))+','+outFile)
                    prevtime = time

                # Cosine score is generated for each library spectrum and returned in a list. See peak_comparison_for_cosine() function description for formatting.
                if len(libKeys) != 0:
                    cosDict, countDict, ionDict, ppmDict = pooled_library_peak_comparison(pool_lib_spectra(lib, libKeys), libKeys, spec['m/z array'], spec['intensity array'], ppmTol, ppmYOffset)
                else: continue

                # Cosine scores and supplementary data are written to the output file.
                # Note that this loop essentially goes over every library spectrum identified by the lib_mz_match_query_window() function.
                for key in libKeys:
                    # Library spectra that had too few matching peaks are excluded. numPeakMatch variable determines the threshold.
                    if countDict[key] > 2:
                        cosine = cosine_similarity(cosDict[key])
                        ionCount = sum([ spec['intensity array'][j]+ppm_offset(spec['intensity array'][j],ppmYOffset) for j in ionDict[key] ])
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
                            spec['precursorMz'][0]['windowWideness'], #totalWindowWidth
                            (countDict[key]**(1/5))*cosine
                        ]
                        writer.writerow(temp)
                        ppmWriter.writerow([spec['num'], key[1], lib[key]['ProteinName']] + sorted(ppmDict[key]))

        # Prints the final number of experimental spectra analyzed.
        print('#'+str(count))

'''
DEFUNCT: See pooled_library_query_spectra_analysis()
Function: pooled_library_peak_comparison()
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
def pooled_library_peak_comparison(libSpectrum, libKeys, expMz, expIntensity, ppmTol, ppmYOffset):
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
