from psims.mzml.writer import MzMLWriter

class Spectrum(object):
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


def organize_query_scan_data_for_writing(mzWindowSpectraBreakdown):
    spectra = []
    for cosineGroup in spectraBreakdown:
        id = f'scan={cosineGroup["scan"]}'
        mzs = [mz for spectrum in cosineGroup['querySpectra'] for mz in spectrum['mzs']]
        intensities = [intensity for spectrum in cosineGroup['querySpectra'] for intensity in spectrum['intensities']]
        spectra.append(Spectrum(id, mzs, intensities, cosineGroup['queryScanMz'], 0, cosineGroup['windowWidth'], f'-{cosineGroup["scan"]}'))
    return spectra

def write_query_scan_data_to_mzml_file(queryScanData, filePath):
    with MzMLWriter(open(filePath, "wb")) as outWriter:
        initialize_mzml_file_parameters(outWriter)
        with outWriter.run(id="csodiaq_test_analysis"):
            with outWriter.spectrum_list(count=len(queryScanData)):
                write_spectra_to_file(queryScanData, outWriter)

def initialize_mzml_file_parameters(outWriter):
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

def write_spectra_to_file(queryScanData, outWriter):
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