import pandas as pd
import os
import numpy as np
from psims.mzml.writer import MzMLWriter

class Spectrum(object):
    def __init__(
        self,
        id,
        mz_array,
        intensity_array,
        precursor_mz,
        precursor_charge,
        window_width,

    ):
        self.id = id
        self.mz_array = mz_array
        self.intensity_array = intensity_array
        self.precursor_mz = precursor_mz
        self.precursor_charge = precursor_charge
        self.window_width = window_width

def get_parent_dir():
    return os.path.dirname(os.path.abspath(__file__))

def make_expected_output_and_query_spectra_files(libraryPath):
    spectraBreakdown = make_library_spectra_cosine_breakdown(libraryPath)
    outputDf = make_output_df_from_spectra_breakdown(spectraBreakdown)
    write_mzml_query_file(spectraBreakdown)

def make_library_spectra_cosine_breakdown(libraryPath):
    librarySpectra = make_target_decoy_library_dict(libraryPath)
    return make_spectra_breakdown_from_library_spectra(librarySpectra)

def make_target_decoy_library_dict(libraryPath):
    libraryDf = pd.read_csv(libraryPath)
    targets = []
    decoys = []
    for vars, df in libraryDf.groupby([
        "precursorMz",
        "peptide",
        "precursorCharge",
        "id",
        "protein",
    ]):
        if 'DECOY' in vars[-1]:
            peptideList = decoys
        else:
            peptideList = targets
        tempDict = {}
        tempDict["precursorMz"] = vars[0]
        tempDict["peptide"] = vars[1]
        tempDict["precursorCharge"] = vars[2]
        tempDict["id"] = vars[3]
        tempDict["protein"] = vars[4]
        tempDict["mzs"] = np.array(df["peakMz"])
        tempDict["intensities"] = np.array(df["peakIntensity"])
        peptideList.append(tempDict)
    return {'targets':targets, 'decoys':decoys}

def make_spectra_breakdown_from_library_spectra(librarySpectra):
    spectraBreakdown = [
        {
            'title': 'highCosineTargets',
            'scan': 1,
            'cosine': 0.99,
            'libSpectra': librarySpectra['targets'][:200]
        },
        {
            'title': 'highCosineDecoys',
            'scan': 2,
            'cosine': 0.98,
            'libSpectra': librarySpectra['decoys'][:2]
        },
        {
            'title': 'midCosineTargets',
            'scan': 3,
            'cosine': 0.97,
            'libSpectra': librarySpectra['targets'][200:400]
        },
        {
            'title': 'midCosineDecoys',
            'scan': 4,
            'cosine': 0.96,
            'libSpectra': librarySpectra['decoys'][2:5]
        },
        {
            'title': 'lowCosineTargets',
            'scan': 5,
            'cosine': 0.95,
            'libSpectra': librarySpectra['targets'][400:]
        },
        {
            'title': 'lowCosineDecoys',
            'scan': 6,
            'cosine': 0.94,
            'libSpectra': librarySpectra['decoys'][5:]
        },
    ]
    for cosineGroup in spectraBreakdown:
        cosineGroup['queryScanMz'], cosineGroup['windowWidth'] = find_mz_window_of_query_scan_from_spectra(cosineGroup['libSpectra'])
        cosineGroup['querySpectra'] = [create_query_spectrum_from_library_spectrum_and_cosine_score(spectrum, cosineGroup['cosine']) for spectrum in cosineGroup['libSpectra']]
    return spectraBreakdown

def make_output_df_from_spectra_breakdown(spectraBreakdown):
    data = []
    for cosineGroup in spectraBreakdown:
        cosine = cosineGroup['cosine']
        numQueryPeaks = find_number_of_total_query_peaks_in_scan(cosineGroup['querySpectra'])
        for i in range(len(cosineGroup["libSpectra"])):
            libSpectra = cosineGroup["libSpectra"][i]
            querySpectra = cosineGroup["querySpectra"][i]
            queryMzLargerThanPrecursorIdx = np.argwhere(np.array(querySpectra['mzs']) > libSpectra['precursorMz'])
            queryMzLargerThanPrecursorIntensities = np.array(querySpectra['intensities'])[queryMzLargerThanPrecursorIdx]
            if 'DECOY' in libSpectra["protein"]: decoy = 1
            else: decoy = 0
            data.append([
                cosineGroup['scan'],
                cosineGroup['queryScanMz'],
                libSpectra["peptide"],
                libSpectra["protein"],
                decoy,
                libSpectra["precursorMz"],
                libSpectra["precursorCharge"],
                cosine,
                libSpectra["id"],
                numQueryPeaks,
                len(libSpectra["mzs"]),
                len(libSpectra["mzs"]),
                sum(queryMzLargerThanPrecursorIntensities)[0],
                f'-{cosineGroup["scan"]}',
                cosineGroup['windowWidth'],
                len(libSpectra["mzs"]) - len(queryMzLargerThanPrecursorIntensities)
            ])
    return pd.DataFrame(data, columns = [
        "scan",
        "MzEXP",
        "peptide",
        "protein",
        "isDecoy",
        "MzLIB",
        "zLIB",
        "cosine",
        "name",
        "Peaks(Query)",
        "Peaks(Library)",
        "shared",
        "ionCount",
        "CompensationVoltage",
        "totalWindowWidth",
        "exclude_num",
    ])

def find_scan_from_title(title):
    if title == 'highCosineTargets': return 1
    elif title == 'highCosineDecoys': return 2
    elif title == 'midCosineTargets': return 3
    elif title == 'midCosineDecoys': return 4
    elif title == 'lowCosineTargets': return 5
    elif title == 'lowCosineDecoys': return 6
    else:
        raise ValueError('woops! wrong title name!')

def find_mz_window_of_query_scan_from_spectra(spectra):
    allPrecursorMzs = [spectrum['precursorMz'] for spectrum in spectra]
    minMz = min(allPrecursorMzs)
    maxMz = max(allPrecursorMzs)
    centralMz = (minMz + maxMz) / 2
    return (centralMz, (centralMz-minMz)*2+.0001)

def find_number_of_total_query_peaks_in_scan(spectra):
    lengths = [len(spectrum['mzs']) for spectrum in spectra]
    return sum(lengths)

def make_input_query_spectra_from_spectra_breakdown(spectraBreakdown):

    pass

def create_vector_that_can_be_used_to_create_cosine_score(vectorA, cosineScore):
    unitVector = vectorA / np.linalg.norm(vectorA)
    positiveVector = np.array(np.full((len(vectorA)), 1000.0))
    perpendicularVector = positiveVector - positiveVector.dot(unitVector)*unitVector
    perpendicularVector = perpendicularVector / np.linalg.norm(perpendicularVector)
    vectorB = cosineScore*unitVector + np.sqrt(1 - cosineScore**2)*perpendicularVector
    return vectorB

def create_query_spectrum_from_library_spectrum_and_cosine_score(spectrum, cosineScore):
    return {
        'mzs': spectrum['mzs'],
        'intensities': create_vector_that_can_be_used_to_create_cosine_score(spectrum['intensities'], cosineScore)
    }

def write_mzml_query_file(spectraBreakdown):
    spectra = organize_scan_data(spectraBreakdown)
    filePath = os.path.join(get_parent_dir(), 'inputs', 'query.mzML')
    write_spectra_to_file(spectra, filePath)


def organize_scan_data(spectraBreakdown):
    spectra = []
    for cosineGroup in spectraBreakdown:
        id = f'scan={cosineGroup["scan"]}'
        mzs = [mz for spectrum in cosineGroup['querySpectra'] for mz in spectrum['mzs']]
        intensities = [intensity for spectrum in cosineGroup['querySpectra'] for intensity in spectrum['intensities']]
        spectra.append(Spectrum(id, mzs, intensities, cosineGroup['queryScanMz'], 0, cosineGroup['windowWidth']))
    return spectra

def write_spectra_to_file(scans, filePath):
    with MzMLWriter(open(filePath, "wb")) as out:
        out.controlled_vocabularies()

        out.file_description(
            [  # the list of file contents terms
                "MS1 spectrum",
                "MSn spectrum",
                "centroid spectrum",
            ]
        )

        out.software_list(
            [
                {
                    "id": "psims-writer",
                    "version": "0.1.2",
                    "params": [
                        "python-psims",
                    ],
                }
            ]
        )

        source = out.Source(1, ["electrospray ionization", "electrospray inlet"])
        analyzer = out.Analyzer(
            2, ["fourier transform ion cyclotron resonance mass spectrometer"]
        )
        detector = out.Detector(3, ["inductive detector"])
        config = out.InstrumentConfiguration(
            id="IC1", component_list=[source, analyzer, detector], params=["LTQ-FT"]
        )
        out.instrument_configuration_list([config])

        methods = []

        methods.append(
            out.ProcessingMethod(
                order=1,
                software_reference="psims-writer",
                params=[
                    "Gaussian smoothing",
                    "median baseline reduction",
                    "MS:1000035",  # peak picking
                    "Conversion to mzML",
                ],
            )
        )
        processing = out.DataProcessing(methods, id="DP1")
        out.data_processing_list([processing])

        # Open the run and spectrum list sections
        with out.run(id="my_analysis"):
            spectrum_count = len(scans)
            with out.spectrum_list(count=spectrum_count):
                for prod in scans:

                    out.write_spectrum(
                        prod.mz_array,
                        prod.intensity_array,
                        id=prod.id,
                        params=[
                            "MSn Spectrum",
                            {"ms level": 2},
                        ],
                        # Include precursor information
                        precursor_information={
                            "mz": prod.precursor_mz,
                            "charge": prod.precursor_charge,
                            "activation": [
                                "beam-type collisional dissociation",
                                {"collision energy": 25},
                            ],
                            "isolation_window": [
                                prod.window_width,
                            ],
                        },
                    )



if __name__ == "__main__":
    libraryPath = os.path.join(get_parent_dir(), 'libraries', 'template_test_library.csv')
    make_expected_output_and_query_spectra_files(libraryPath)

