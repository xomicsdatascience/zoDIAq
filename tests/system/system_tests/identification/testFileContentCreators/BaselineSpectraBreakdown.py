import os
import numpy as np
import pandas as pd
from psims.mzml.writer import MzMLWriter
import pyopenms


class BaselineSpectraBreakdown:
    def __init__(self, templateLibraryDf):
        targetDecoyLibrarySpectra = _separate_target_and_decoy_library_spectra(
            templateLibraryDf
        )
        librarySpectraBreakdown = self.make_library_spectra_breakdown_summary(
            targetDecoyLibrarySpectra
        )
        mzWindowSpectraBreakdown = (
            self.add_query_file_components_to_library_spectra_summary(
                librarySpectraBreakdown
            )
        )
        self.expectedOutputDf = self.make_expected_output_df(mzWindowSpectraBreakdown)
        self.inputSpectra = self.make_input_query_spectra(mzWindowSpectraBreakdown)

    def make_library_spectra_breakdown_summary(self, targetDecoyLibrarySpectra):
        return [
            {
                "title": "highCosineTargets",
                "scan": 1,
                "cosine": 0.99,
                "libSpectra": targetDecoyLibrarySpectra["targets"][:200],
            },
            {
                "title": "highCosineDecoys",
                "scan": 2,
                "cosine": 0.98,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][:2],
            },
            {
                "title": "midCosineTargets",
                "scan": 3,
                "cosine": 0.97,
                "libSpectra": targetDecoyLibrarySpectra["targets"][200:400],
            },
            {
                "title": "midCosineDecoys",
                "scan": 4,
                "cosine": 0.96,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][2:5],
            },
            {
                "title": "lowCosineTargets",
                "scan": 5,
                "cosine": 0.95,
                "libSpectra": targetDecoyLibrarySpectra["targets"][400:],
            },
            {
                "title": "lowCosineDecoys",
                "scan": 6,
                "cosine": 0.94,
                "libSpectra": targetDecoyLibrarySpectra["decoys"][5:],
            },
        ]

    def add_query_file_components_to_library_spectra_summary(
        self, librarySpectraBreakdown
    ):
        for mzWindowGroup in librarySpectraBreakdown:
            (
                mzWindowGroup["queryScanMz"],
                mzWindowGroup["windowWidth"],
            ) = find_mz_window_of_query_scan_from_spectra(mzWindowGroup["libSpectra"])
            mzWindowGroup["querySpectra"] = [
                self.create_query_spectrum_from_library_spectrum_and_cosine_score(
                    spectrum, mzWindowGroup
                )
                for spectrum in mzWindowGroup["libSpectra"]
            ]
        return librarySpectraBreakdown

    def create_query_spectrum_from_library_spectrum_and_cosine_score(
        self, spectrum, mzWindowGroup
    ):
        queryMz = self.format_query_mz_values(spectrum["mzs"], mzWindowGroup)
        return {
            "mzs": queryMz,
            "intensities": create_vector_that_can_be_used_to_create_cosine_score(
                spectrum["intensities"], mzWindowGroup["cosine"]
            ),
        }

    def format_query_mz_values(self, libraryMzs, mzWindowGroup):
        return libraryMzs

    def make_expected_output_df(self, mzWindowSpectraBreakdown):
        data = []
        for cosineGroup in mzWindowSpectraBreakdown:
            cosine = cosineGroup["cosine"]
            numQueryPeaks = find_number_of_total_query_peaks_in_scan(
                cosineGroup["querySpectra"]
            )
            for i in range(len(cosineGroup["libSpectra"])):
                libSpectra = cosineGroup["libSpectra"][i]
                querySpectra = cosineGroup["querySpectra"][i]
                queryMzLargerThanPrecursorIdx = np.argwhere(
                    np.array(querySpectra["mzs"]) > cosineGroup["queryScanMz"]
                )
                queryMzLargerThanPrecursorIntensities = np.array(
                    querySpectra["intensities"]
                )[queryMzLargerThanPrecursorIdx]
                ionCount = 0
                if len(queryMzLargerThanPrecursorIntensities):
                    ionCount = sum(queryMzLargerThanPrecursorIntensities)[0]
                if "DECOY" in libSpectra["protein"]:
                    decoy = 1
                else:
                    decoy = 0
                data.append(
                    [
                        cosineGroup["scan"],
                        cosineGroup["queryScanMz"],
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
                        ionCount,
                        -cosineGroup["scan"],
                        cosineGroup["windowWidth"],
                        len(libSpectra["mzs"])
                        - len(queryMzLargerThanPrecursorIntensities),
                    ]
                )
        return pd.DataFrame(
            data,
            columns=[
                "scan",
                "MzEXP",
                "peptide",
                "protein",
                "isDecoy",
                "MzLIB",
                "zLIB",
                "cosine",
                "name",
                "Peak(Query)",
                "Peaks(Library)",
                "shared",
                "ionCount",
                "CompensationVoltage",
                "totalWindowWidth",
                "exclude_num",
            ],
        )

    def make_input_query_spectra(self, mzWindowSpectraBreakdown):
        spectra = []
        for cosineGroup in mzWindowSpectraBreakdown:
            id = f'scan={cosineGroup["scan"]}'
            mzs = [
                mz for spectrum in cosineGroup["querySpectra"] for mz in spectrum["mzs"]
            ]
            intensities = [
                intensity
                for spectrum in cosineGroup["querySpectra"]
                for intensity in spectrum["intensities"]
            ]
            spectra.append(
                _Spectrum(
                    id,
                    mzs,
                    intensities,
                    cosineGroup["queryScanMz"],
                    0,
                    cosineGroup["windowWidth"],
                    f'-{cosineGroup["scan"]}',
                )
            )
        return spectra

    def write_query_scan_data_input_files(self, inputDir, testHeader):
        mzmlFile = os.path.join(inputDir, f"{testHeader}.mzML")
        _write_query_scan_data_to_mzml_file(self.inputSpectra, mzmlFile)
        _convert_mzml_file_to_mzxml_file(mzmlFile)


class _Spectrum(object):
    def __init__(
        self,
        id,
        mzArray,
        intensityArray,
        precursorMz,
        precursorCharge,
        windowWidth,
        compensationVoltage=None,
    ):
        self.id = id
        self.mzArray = mzArray
        self.intensityArray = intensityArray
        self.precursorMz = precursorMz
        self.precursorCharge = precursorCharge
        self.windowWidth = windowWidth
        self.compensationVoltage = compensationVoltage


def _separate_target_and_decoy_library_spectra(libraryDf):
    targets = []
    decoys = []
    for vars, df in libraryDf.groupby(
        [
            "precursorMz",
            "peptide",
            "precursorCharge",
            "id",
            "protein",
        ]
    ):
        if "DECOY" in vars[-1]:
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
    return {"targets": targets, "decoys": decoys}


def find_mz_window_of_query_scan_from_spectra(spectra):
    allPrecursorMzs = [spectrum["precursorMz"] for spectrum in spectra]
    minMz = min(allPrecursorMzs)
    maxMz = max(allPrecursorMzs)
    centralMz = (minMz + maxMz) / 2
    return (centralMz, (centralMz - minMz) * 2 + 0.0001)


def find_number_of_total_query_peaks_in_scan(spectra):
    lengths = [len(spectrum["mzs"]) for spectrum in spectra]
    return sum(lengths)


def create_vector_that_can_be_used_to_create_cosine_score(vectorA, cosineScore):
    vectorA = np.sqrt(vectorA)
    unitVector = vectorA / np.linalg.norm(vectorA)
    positiveVector = np.array(np.full((len(vectorA)), 1000.0))
    perpendicularVector = positiveVector - positiveVector.dot(unitVector) * unitVector
    perpendicularVector = perpendicularVector / np.linalg.norm(perpendicularVector)
    vectorB = (
        cosineScore * unitVector + np.sqrt(1 - cosineScore**2) * perpendicularVector
    )
    return vectorB**2


def calculate_mz_value_from_ppm(mz, ppm):
    return mz * (1e6 - ppm) / 1e6


def _write_query_scan_data_to_mzml_file(queryScanData, filePath):
    with MzMLWriter(open(filePath, "wb")) as outWriter:
        _initialize_mzml_file_parameters(outWriter)
        with outWriter.run(id="zodiaq_test_analysis"):
            with outWriter.spectrum_list(count=len(queryScanData)):
                _write_spectra_to_file(queryScanData, outWriter)


def _initialize_mzml_file_parameters(outWriter):
    outWriter.controlled_vocabularies()
    outWriter.file_description(
        [
            "MS1 spectrum",
            "MSn spectrum",
            "centroid spectrum",
        ]
    )

    outWriter.software_list(
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

    source = outWriter.Source(1, ["electrospray ionization", "electrospray inlet"])
    analyzer = outWriter.Analyzer(
        2, ["fourier transform ion cyclotron resonance mass spectrometer"]
    )
    detector = outWriter.Detector(3, ["inductive detector"])
    config = outWriter.InstrumentConfiguration(
        id="IC1", component_list=[source, analyzer, detector], params=["LTQ-FT"]
    )
    outWriter.instrument_configuration_list([config])

    methods = []

    methods.append(
        outWriter.ProcessingMethod(
            order=1,
            software_reference="psims-writer",
            params=[
                "Gaussian smoothing",
                "median baseline reduction",
                "MS:1000035",
                "Conversion to mzML",
            ],
        )
    )
    processing = outWriter.DataProcessing(methods, id="DP1")
    outWriter.data_processing_list([processing])


def _write_spectra_to_file(queryScanData, outWriter):
    for spectrum in queryScanData:
        outWriter.write_spectrum(
            spectrum.mzArray,
            spectrum.intensityArray,
            id=spectrum.id,
            params=[
                "MSn Spectrum",
                {"ms level": 2},
                {"compensationVoltage": spectrum.compensationVoltage},
            ],
            precursor_information={
                "mz": spectrum.precursorMz,
                "charge": spectrum.precursorCharge,
                "activation": [
                    "beam-type collisional dissociation",
                    {"collision energy": 25},
                ],
                "isolation_window": [
                    spectrum.windowWidth,
                ],
            },
        )


def _convert_mzml_file_to_mzxml_file(mzmlFilePath):
    msExperiment = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(mzmlFilePath, msExperiment)
    pyopenms.MzXMLFile().store(
        f"{os.path.splitext(mzmlFilePath)[0]}.mzXML", msExperiment
    )
