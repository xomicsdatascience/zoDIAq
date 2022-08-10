import pandas as pd
from pyteomics import mzxml, mgf, mass
from bisect import bisect, bisect_left
from timeit import default_timer as timer
import csv
import numpy as np
from . import idpicker as idp
import re
from collections import defaultdict
from numba import njit
from . import IdentificationSpectraMatcher
from . import spectra_matcher_functions as smf


def library_file_to_dict(inFile):
    fileType = inFile.split('.')[-1]
    if fileType == 'mgf':
        lib = mgf_library_upload(inFile)
    else:
        lib = traml_library_upload(inFile)
    return lib


def perform_spectra_pooling_and_analysis(querySpectraFile, outFile, lib, tolerance, maxQuerySpectraToPool, corrected, histFile):

    smf.print_milestone('Begin Grouping Scans by m/z Windows:')
    queWindowDict, queScanValuesDict = pool_scans_by_mz_windows(
        querySpectraFile)

    print('Number of Unpooled MS/MS Query Spectra: ' +
          str(len(queScanValuesDict)))
    print('Number of Pooled MS/MS Query Spectra/Mz Windows: ' +
          str(len(queWindowDict)), flush=True)

    # To enhance the print experience, status prints will be given at intervals tailored to the number of identified windows.
    #  example: if there are 1-99 pooled query spectra, print statements are made after every pooled query spectra analysis is complete.
    #           if there are 100-999, print after every 10 pooled spectra. And so on.
    printFriendlyCounter = 100
    while printFriendlyCounter < len(queWindowDict):
        printFriendlyCounter *= 10
    printFriendlyCounter /= 100

    allLibKeys, libIdToKeyDict, libIdToDecoyDict = gather_library_metadata(lib)
    allSpectraMatches = IdentificationSpectraMatcher.IdentificationSpectraMatcher()
    numWindowsAnalyzed = 0

    prevtime = timer()
    smf.print_milestone('Begin Pooled Spectra Analysis:')
    with mzxml.read(querySpectraFile, use_index=True) as spectra:

        for precMz_win, scans in queWindowDict.items():
            top_mz = precMz_win[0] + precMz_win[1] / 2
            bottom_mz = precMz_win[0] - precMz_win[1] / 2
            libKeys = identify_lib_spectra_in_window(
                top_mz, bottom_mz, allLibKeys)
            if len(libKeys) == 0:
                continue
            pooledLibSpectra = pool_lib_spectra(lib, libKeys)
            pooledQueSpectra = []

            for i in range(len(scans)):
                scanNumber = scans[i]
                queSpectrum = spectra.get_by_id(scanNumber)
                pooledQueSpectra += smf.format_spectra_for_pooling(
                    queSpectrum, scanNumber)

                if (i % maxQuerySpectraToPool == 0 and i != 0) or i == len(scans)-1:
                    pooledQueSpectra.sort()
                    windowSpectraMatches = IdentificationSpectraMatcher.IdentificationSpectraMatcher()
                    windowSpectraMatches.compare_spectra(
                        pooledLibSpectra, pooledQueSpectra, tolerance, libIdToDecoyDict)
                    allSpectraMatches.extend_all_spectra(windowSpectraMatches)
                    pooledQueSpectra.clear()

            numWindowsAnalyzed += 1
            if numWindowsAnalyzed % printFriendlyCounter == 0:
                time = timer()
                print('\nNumber of Pooled Experimental Spectra Analyzed: ' +
                      str(numWindowsAnalyzed))
                print('Number of Spectra in Current Pooled Spectra: ' +
                      str(len(scans)))
                print('Time Since Last Checkpoint: ' +
                      str(round(time-prevtime, 2)) + ' Seconds', flush=True)
                prevtime = time

    smf.print_milestone('Begin FDR Analysis:')
    maccCutoff = allSpectraMatches.find_score_fdr_cutoff()

    if corrected != -1:
        smf.print_milestone('Begin Correction Process:')
        allSpectraMatches.filter_by_corrected_ppm_window(
            corrected, maccCutoff, histFile)

        smf.print_milestone('Begin Corrected FDR Analysis:')
        maccCutoff = allSpectraMatches.find_score_fdr_cutoff()

    smf.print_milestone('\nBegin Writing to File: ')
    allSpectraMatches.write_output(
        outFile, querySpectraFile, maccCutoff, queScanValuesDict, libIdToKeyDict, lib)


def write_fdr_outputs(inFile, specFile, pepFile, protFile):

    smf.print_milestone('Generating FDR Analysis Files:')
    overallDf = pd.read_csv(inFile).sort_values(
        'MaCC_Score', ascending=False).reset_index(drop=True)
    spectralDf = add_fdr_to_csodiaq_output(overallDf)
    peptideDf = add_fdr_to_csodiaq_output(overallDf, filterType='peptide')

    spectralDf.to_csv(specFile, index=False)
    peptideDf.to_csv(pepFile, index=False)

    if protFile:
        peptideProteinConnections = format_peptide_protein_connections(
            peptideDf)
        verifiedProteinDict = idp.find_valid_proteins(
            peptideProteinConnections)
        proteinDf = add_leading_protein_column(peptideDf, verifiedProteinDict)
        tempProtDf = add_fdr_to_csodiaq_output(
            proteinDf, filterType='leadingProtein')
        proteinMetaInfoDict = tempProtDf.set_index(
            'leadingProtein').T.to_dict()
        proteinDf = remove_invalid_peptides_and_add_metadata(
            proteinDf, proteinMetaInfoDict)
        proteinDf = mark_peptides_unique_to_proteins(proteinDf)
        proteinDf.to_csv(protFile, index=False)


def filter_fdr_output_for_targeted_reanalysis(fdrFile, proteins, heavy):
    smf.print_milestone('Generate DISPA Targeted Reanalysis Files:')
    fdrDf = pd.read_csv(fdrFile)
    fdrDf = fdrDf[~fdrDf['protein'].str.contains(
        'DECOY', na=False)].reset_index(drop=True)
    if heavy:
        fdrDf = fdrDf[fdrDf['peptide'].str.endswith(
            'R') | fdrDf['peptide'].str.endswith('K')].reset_index(drop=True)
    if proteins:
        fdrDf = fdrDf[fdrDf['uniquePeptide'] == 1].sort_values(
            'ionCount', ascending=False).reset_index(drop=True)
        fdrDf = fdrDf.groupby(['leadingProtein']).head(proteins).reset_index()
    return fdrDf


def write_targeted_reanalysis_outputs(header, fdrDf, heavy):
    bins = False
    CVs = gather_all_possible_cv_values(fdrDf)
    allDfs = []
    for CV in CVs:
        if CV:
            cvDf = fdrDf[fdrDf['CompensationVoltage']
                         == CV].reset_index(drop=True)
        else:
            cvDf = fdrDf
        if bins:
            cvDf, data = compile_reanalysis_files_with_bins(cvDf, heavy)
        else:
            cvDf, data = compile_reanalysis_files_without_bins(cvDf, heavy)
        allDfs.append(cvDf)
        finalDf = pd.DataFrame(
            data, columns=['Compound', 'Formula', 'Adduct', 'm.z', 'z', 'MSXID'])
        outFile = header
        if CV:
            outFile += '_' + str(CV)
        outFile += '.txt'
        finalDf.to_csv(outFile, sep='\t', index=False)

    compDf = pd.concat(allDfs)
    if bins:
        binMark = '_withBins_'
    else:
        binMark = '_withoutBins_'
    compDf.to_csv(header+binMark+'allCVs.csv', index=False)


def traml_library_upload(fileName):
    if fileName.endswith('.tsv'):
        lib_df = pd.read_csv(fileName, sep='\t')
    else:
        lib_df = pd.read_csv(fileName)
    smf.print_milestone('Enter library dictionary upload: ')

    # Pan human and spectraST libraries have different column names. This normalizes the columns.
    headings = traml_column_headings(lib_df.columns)
    lib_df = lib_df.loc[:, lib_df.columns.intersection([headings['PrecursorMz'], headings['FullUniModPeptideName'], headings['PrecursorCharge'],
                                                       headings['ProductMz'], headings['LibraryIntensity'], headings['transition_group_id'], headings['ProteinName']])]
    lib_df = lib_df[[headings['PrecursorMz'], headings['FullUniModPeptideName'], headings['PrecursorCharge'],
                     headings['ProductMz'], headings['LibraryIntensity'], headings['transition_group_id'], headings['ProteinName']]]
    lib_df.columns = ['PrecursorMz', 'FullUniModPeptideName', 'PrecursorCharge',
                      'ProductMz', 'LibraryIntensity', 'transition_group_id', 'ProteinName']

    lib_df['LibraryIntensity'] = [
        x**0.5 for x in list(lib_df['LibraryIntensity'])]
    lib_df['ID'] = list(zip(lib_df['PrecursorMz'].tolist(),
                        lib_df['FullUniModPeptideName'].tolist()))

    mz_dict = lib_df.groupby("ID")['ProductMz'].apply(list).to_dict()
    intensity_dict = lib_df.groupby(
        "ID")['LibraryIntensity'].apply(list).to_dict()
    lib_df.drop_duplicates(subset="ID", inplace=True)
    lib_df = lib_df.loc[:, lib_df.columns.intersection(
        ['ID', 'PrecursorCharge', 'transition_group_id', 'ProteinName'])]
    lib_df.set_index("ID", drop=True, inplace=True)
    lib = lib_df.to_dict(orient="index")

    # pan human library formats are different, including how the peptides are matched to proteins (esp. decoys). This section of code adjusts for this discrepancy.
    if headings['type'] == 'PanHuman':
        for key, value in lib.items():
            proteins = lib[key]['ProteinName'].split('/')
            num = proteins.pop(0)
            newProteins = [x for x in proteins if 'DECOY' not in x]
            proteinStr = str(len(newProteins))
            for x in newProteins:
                if 'DECOY' in num:
                    proteinStr += ('/DECOY_'+x)
                else:
                    proteinStr += ('/'+x)
            lib[key]['ProteinName'] = proteinStr

    id = 0
    for key in lib:
        id += 1
        mz, intensity = (list(t) for t in zip(
            *sorted(zip(mz_dict[key], intensity_dict[key]))))
        keyList = [id for i in range(len(mz))]
        peaks = list(tuple(zip(mz, intensity, keyList)))

        peaks.sort(key=lambda x: x[1], reverse=True)
        if len(peaks) > 10:
            peaks = peaks[:10]

        peaks.sort(key=lambda x: x[0])
        lib[key]['Peaks'] = peaks
        lib[key]['ID'] = id
        if 'DECOY' in lib[key]['ProteinName']:
            lib[key]['Decoy'] = 1
        else:
            lib[key]['Decoy'] = 0
    return lib


def traml_column_headings(columns):
    if 'FullUniModPeptideName' in columns:  # SpectraST
        return {
            'type': 'SpectraST',
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'FullUniModPeptideName',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'transition_group_id',
            'ProteinName': 'ProteinName',
        }
    else:  # pan human
        return {
            'type': 'PanHuman',
            'PrecursorMz': 'PrecursorMz',
            'FullUniModPeptideName': 'ModifiedPeptideSequence',
            'PrecursorCharge': 'PrecursorCharge',
            'ProductMz': 'ProductMz',
            'LibraryIntensity': 'LibraryIntensity',
            'transition_group_id': 'TransitionGroupId',
            'ProteinName': 'ProteinId',
        }


def mgf_library_upload(fileName):
    libMGF = mgf.read(fileName)
    smf.print_milestone('Enter library dictionary upload: ')
    lib = {}
    id = 0
    for spec in libMGF:
        id += 1
        key = (spec['params']['pepmass'][0], spec['params']['seq'])
        charge = int(re.sub('[+-]', '', str(spec['params']['charge'][0])))
        name = spec['params']['title']
        if 'protein' in spec['params']:
            protein = spec['params']['protein']
        else:
            protein = ''
        if 'DECOY' in name:
            decoy = 1
        else:
            decoy = 0
        mz = spec['m/z array']
        intensity = spec['intensity array']
        intensity = [x**0.5 for x in intensity]
        keyList = [id for x in mz]
        peaks = list(tuple(zip(mz, intensity, keyList)))
        peaks.sort(key=lambda x: x[1], reverse=True)
        if len(peaks) > 10:
            peaks = peaks[:10]
        peaks.sort(key=lambda x: x[0])
        tempDict = {
            'PrecursorCharge': charge,
            'transition_group_id': name,
            'ProteinName': protein,
            'Peaks': peaks,
            'ID': id,
            'Decoy': decoy,
        }
        lib[key] = tempDict
    return lib


def identify_lib_spectra_in_window(top_mz, bottom_mz, sortedLibKeys):
    temp = sortedLibKeys[:]
    top_key = (top_mz, "z")
    bottom_key = (bottom_mz, "")

    i1 = bisect(temp, bottom_key)
    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
    if i2-i1 == 1:
        return []
    temp.insert(i2, top_key)

    return temp[i1+1:i2]


def pool_lib_spectra(lib, libKeys):
    finalList = []
    for key in libKeys:
        finalList += lib[key]['Peaks']
    return sorted(finalList)


def pool_scans_by_mz_windows(querySpectraFile):
    queWindowDict = defaultdict(list)
    queScanValuesDict = defaultdict(dict)

    with mzxml.read(querySpectraFile, use_index=True) as spectra:
        for spec in spectra:

            if 'precursorMz' not in spec:
                continue
            scan = spec['num']
            precMz = spec['precursorMz'][0]['precursorMz']
            windowWidth = spec['precursorMz'][0]['windowWideness']
            queWindowDict[precMz, windowWidth].append(scan)

            queScanValuesDict[scan]['precursorMz'] = precMz
            queScanValuesDict[scan]['windowWideness'] = windowWidth
            queScanValuesDict[scan]['peaksCount'] = spec['peaksCount']
            if 'compensationVoltage' in spec:
                CV = spec['compensationVoltage']
            else:
                CV = ''
            queScanValuesDict[scan]['CV'] = CV

    return queWindowDict, queScanValuesDict


def gather_library_metadata(lib):
    allLibKeys = sorted(lib.keys())
    libIdToKeyDict = {}
    libIdToDecoyDict = {}
    for key in allLibKeys:
        libIdToKeyDict[lib[key]['ID']] = key
        libIdToDecoyDict[lib[key]['ID']] = lib[key]['Decoy']
    return allLibKeys, libIdToKeyDict, libIdToDecoyDict


def format_peptide_protein_connections(peptideDf):
    peptideProteinConnections = []

    for i in range(len(peptideDf)):
        peptide = peptideDf['peptide'].loc[i]
        proteinGroup = peptideDf['protein'].loc[i]

        for protein in proteinGroup.split('/')[1:]:
            peptideProteinConnections.append((peptide, protein))
    return peptideProteinConnections


def remove_invalid_peptides_and_add_metadata(proteinDf, proteinMetaInfoDict):
    removables, proteinCosine, proteinFDR = [], [], []

    for i in range(len(proteinDf)):
        protein = proteinDf['leadingProtein'].loc[i]
        if protein in proteinMetaInfoDict:  # dict only contains proteins with FDR < 0.01
            proteinCosine.append(proteinMetaInfoDict[protein]['cosine'])
            proteinFDR.append(
                proteinMetaInfoDict[protein]['leadingProteinFDR'])
        else:
            removables.append(i)
    proteinDf = proteinDf.drop(
        proteinDf.index[removables]).reset_index(drop=True)
    proteinDf['proteinCosine'] = proteinCosine
    proteinDf['leadingProteinFDR'] = proteinFDR
    return proteinDf


def mark_peptides_unique_to_proteins(proteinDf):
    proteinDf = proteinDf.sort_values(['proteinCosine', 'leadingProtein', 'cosine'], ascending=[
                                      False, False, False]).reset_index(drop=True)
    uniquePepsDict = defaultdict(set)
    for i in range(len(proteinDf)):
        uniquePepsDict[proteinDf.loc[i]['peptide']].add(
            proteinDf.loc[i]['leadingProtein'])

    uniquePeps = []
    for i in range(len(proteinDf)):
        p = proteinDf.loc[i]['peptide']
        if len(uniquePepsDict[proteinDf.loc[i]['peptide']]) == 1:
            uniquePeps.append(1)
        else:
            uniquePeps.append(0)

    proteinDf['uniquePeptide'] = uniquePeps
    return proteinDf


def add_fdr_to_csodiaq_output(df, filterType='spectral', bestMatchNum=0):
    finalDf = df.copy()
    if filterType != 'spectral':
        finalDf = finalDf.drop_duplicates(
            subset=filterType, keep='first').reset_index(drop=True)
    fdrList, decoyNum = fdr_calculation(finalDf)
    finalDf = finalDf.truncate(after=len(fdrList)-1)
    finalDf[filterType + 'FDR'] = fdrList
    finalDf = finalDf.reset_index(drop=True)
    return finalDf


def add_leading_protein_column(df, verifiedProteinDict):
    finalDf = pd.DataFrame(columns=df.columns)
    leadingProteins = []
    for i in range(len(df)):
        proteinGroup = df['protein'].loc[i]
        proteins = proteinGroup.split('/')[1:]
        for protein in proteins:
            if protein in verifiedProteinDict:
                finalDf = finalDf.append(df.loc[i])
                leadingProteins.append(verifiedProteinDict[protein])
    finalDf['leadingProtein'] = leadingProteins
    finalDf = finalDf.reset_index(drop=True)
    finalDf = finalDf.drop_duplicates(keep='first').reset_index(drop=True)
    return finalDf


def fdr_calculation(df):  # ***NOTE***make return type
    fdrValues = []
    indices = []
    numDecoys = 0
    df.fillna("nan", inplace=True)
    count = 0
    for index, row in df.iterrows():
        # current criteria for 'decoys' is to have 'decoy' in the protein name. This may change in the future.
        if 'DECOY' in row['protein']:
            numDecoys += 1

        # calculates the FDR up to this point in the data frame.
        curFDR = numDecoys/(count+1)

        # conditional statement comparing the current FDR to the FDR Cutoff. If larger, function values are returned.
        if curFDR > 0.01:

            # if the number of rows has not yet reached the minimum number that allows for the FDR cutoff, 0 is returned instead.
            if len(fdrValues) < 1/0.01:
                return [], 0
            return fdrValues, numDecoys-1
        fdrValues.append(curFDR)
        indices.append(index)
        count += 1

    return fdrValues, numDecoys-1


def gather_all_possible_cv_values(fdrDf):
    CVs = set(fdrDf['CompensationVoltage'])
    # get list of all CVs, or CV==0 if there are no CVs
    def notNan(x): return ~np.isnan(x)
    CVs = set(filter(notNan, CVs))
    if len(CVs) == 0:
        CVs.add('')
    return CVs


def calculate_all_light_heavy_mzs(cvDf):
    lightMzs = []
    heavyMzs = []
    peptides = []
    for index, row in cvDf.iterrows():
        peptide = row['peptide']
        lightMz = float(row['MzLIB'])
        charge = row['zLIB']
        heavyMz = smf.calculate_heavy_mz(peptide, lightMz, charge)
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

    binDict = defaultdict(list)
    for i in range(len(cvDf)):
        if heavy:
            binDict[scanLightMzs[i], scanHeavyMzs[i]].append(peptides[i])
        else:
            binDict[scanLightMzs[i], 0].append(peptides[i])

    data = []
    binKeys = sorted(binDict)
    for i in range(len(binKeys)):
        data.append(['/'.join(binDict[binKeys[i]]), '',
                    '(no adduct)', binKeys[i][0]-(binWidth/2), 2, i+1])
        if heavy:
            data.append(['/'.join(binDict[binKeys[i]]), '',
                        '(no adduct)', binKeys[i][1]-(binWidth/2), 2, i+1])
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
        heavyMz = smf.calculate_heavy_mz(compound, lightMz, charge)
        MSXID = i+1
        scanLightMzs.append(round(lightMz, ndigits=2))
        scanHeavyMzs.append(round(heavyMz, ndigits=2))
        data.append([compound, formula, adduct,
                    scanLightMzs[-1], charge, MSXID])
        if heavy:
            data.append([compound, formula, adduct,
                        scanHeavyMzs[-1], charge, MSXID])

    cvDf['scanLightMzs'] = scanLightMzs
    cvDf['scanHeavyMzs'] = scanHeavyMzs
    return cvDf, data


def bin_assignment(mzValues, binWidths):
    maxi = max(mzValues)+binWidths
    mini = min(mzValues)
    bins = np.arange(mini, maxi, binWidths)
    bins = [round(x, 2) for x in bins]
    values = np.digitize(mzValues, bins)
    return values, bins
