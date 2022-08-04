import pandas as pd
from pyteomics import mzxml, mgf, mass
import numpy as np
import re
from collections import defaultdict
from numba import njit
from . import QuantificationSpectraMatcher
from . import spectra_matcher_functions as smf


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
            i = smf.approx_list(lightfragmz, mz)
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
            peaks.append((smf.calculate_heavy_mz(fragList[i][2], fragMz, 1), fragInt, (1, i, seq)))

        peaks.sort(key=lambda x:x[0])

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
            heavyMz = smf.calculate_heavy_mz(fragSeq, fragMz, lib[key]['FragmentCharge'])
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
        seq, mz, z, CV = fragDf.loc[i]['peptide'], fragDf.loc[i]['MzLIB'], fragDf.loc[i]['zLIB'], fragDf.loc[i]['CompensationVoltage']
        lightMz = round(mz, 2)
        key = (round(fragDf.loc[i]['scanLightMzs'],2), round(fragDf.loc[i]['scanHeavyMzs'],2), CV)
        if key in mzxmlToCsodiaqDict:
            scan = mzxmlToCsodiaqDict[key]
            libScanDict[lightMz, seq] = scan
            fragVarDict[scan].append({'seq':seq, 'mz':mz, 'z':z, 'CV':CV})
    return fragVarDict, libScanDict

def connect_mzxml_to_csodiaq_and_library(idFile, libFile, mzxmlFiles, maxPeaks):
    smf.print_milestone('Preparing Quantification Dictionaries:')
    metadataToScanDict = create_mzxml_to_csodiaq_dict(mzxmlFiles[0])
    fileType = idFile.split('.')[-1]
    if fileType == 'csv': fragDf = pd.read_csv(idFile)
    else: fragDf = pd.read_csv(idFile, sep='\t')

    scanToCsodiaqDict, libMetadataToScanDict = connect_csodiaq_data_to_scans(idFile, metadataToScanDict, fragDf)
    scanToLibPeaksDict = pool_library_spectra_by_scan(libFile,libMetadataToScanDict, fragDf, maxPeaks)
    return scanToCsodiaqDict, scanToLibPeaksDict


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

    for f in mzxmlFiles:
        ppmDiffs = []
        allSpectraMatch = QuantificationSpectraMatcher.QuantificationSpectraMatcher()
        scanToNoiseIntensityCutoffDict = dict()
        with mzxml.read(f, use_index =True) as file:
            for scan in sorted(libDict.keys()):

                spec = file.get_by_id(scan)

                scanToNoiseIntensityCutoffDict[int(scan)] = np.mean(sorted(spec['intensity array'])[:10])/2

                expSpectrum = smf.format_spectra_for_pooling(spec, scan, sqrt=False)
                expSpectrum.sort()

                libSpectra = sorted(libDict[scan])

                quantSpectraMatch = QuantificationSpectraMatcher.QuantificationSpectraMatcher()
                quantSpectraMatch.compare_spectra(libSpectra, expSpectrum, massTol, minMatch)
                allSpectraMatch.extend_all_spectra(quantSpectraMatch)

        if correction != -1: allSpectraMatch.filter_by_corrected_ppm_window(correction, hist, minMatch)
        ratioDict = defaultdict(initialize_ratio_dict_values)
        if len(allSpectraMatch.libraryIntensities) != 0: ratioDict = allSpectraMatch.determine_ratios(ratioDict, scanToNoiseIntensityCutoffDict, ratioType, minMatch)


        finalDf[f] = [ratioDict[(int(row['scan']),row['peptide'])] for index, row in finalDf.iterrows()]
    smf.print_milestone('Finish SILAC Quantification')
    return finalDf
