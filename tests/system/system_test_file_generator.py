class Spectrum(object):
    def __init__(
        self,
        id,
        mz_array,
        intensity_array,
        precursor_mz=None,
        precursor_charge=None,
        precursor_intensity=None,
        precursor_scan_id=None,
    ):
        self.id = id
        self.mz_array = mz_array
        self.intensity_array = intensity_array
        self.precursor_mz = precursor_mz
        self.precursor_intensity = precursor_intensity
        self.precursor_charge = precursor_charge
        self.precursor_scan_id = precursor_scan_id


mz_array = [
    255.22935009,
    283.26141863,
    284.26105318,
    301.23572871,
    304.908247,
    329.26327093,
    755.25267878,
    910.6317435,
    960.96747136,
    971.31396162,
    972.649568,
    991.66036894,
    1017.87649113,
    1060.29899182,
    1112.67519902,
    1113.11545762,
    1113.86200673,
    1114.42377982,
    1152.34596544,
    1177.73119994,
    1188.36935517,
    1214.70161813,
    1265.9795606,
    1266.16111855,
    1293.14606767,
    1294.68263447,
    1367.01605133,
    1565.95282753,
    1700.23290184,
    969.65408783,
    1110.37794027,
    1170.89893785,
    1175.34669421,
    1183.40737076,
    861.61958381,
    1292.94207114,
    1295.4429046,
    876.96085225,
    1335.39755355,
    1357.92354342,
    1365.96972386,
    925.64504217,
    958.65480011,
    1438.49014079,
    1452.48402739,
    967.98859816,
    986.63557879,
    1480.46326372,
    1001.97032019,
    1007.67089513,
    1525.51932956,
    1016.67747057,
    1080.3583722,
    1090.03199733,
    1133.01762219,
    1186.72735997,
    960.79647487,
    1274.44354121,
    1918.2693496,
    961.85503396,
]

intensity_array = [
    1.90348869e03,
    1.92160377e03,
    3.26032338e02,
    1.05527732e03,
    9.50991606e02,
    1.52403574e03,
    8.63154395e02,
    4.11169655e02,
    2.33462730e03,
    2.62603673e02,
    2.73669694e02,
    8.62436899e02,
    4.22323174e02,
    2.54371429e02,
    1.02364420e03,
    5.44244205e02,
    4.93101348e02,
    2.64984906e02,
    9.36500725e02,
    4.79373626e02,
    9.26742857e02,
    4.52209221e02,
    3.02178809e03,
    4.94385979e02,
    1.67240655e03,
    9.41320838e02,
    7.25744090e02,
    1.27260012e03,
    8.24236545e02,
    4.93518583e02,
    7.33281806e04,
    9.60817582e02,
    2.64893187e03,
    7.95614286e03,
    4.94586514e03,
    2.35346153e04,
    1.50301526e03,
    5.36435167e03,
    3.78042332e02,
    1.09926345e03,
    3.10133857e03,
    7.41566590e02,
    2.77340229e05,
    1.00796887e05,
    4.69519356e03,
    1.55343822e04,
    5.45621612e03,
    5.53939031e03,
    9.49732490e03,
    8.05000735e03,
    2.65457068e03,
    1.36766228e04,
    2.69348480e03,
    6.71802368e03,
    4.46828571e02,
    1.39065143e04,
    4.29267365e03,
    2.73782365e03,
    1.35373492e03,
    1.17601397e03,
]

scan_data = [
    [
        Spectrum("scan=1", mz_array, intensity_array),
        (Spectrum("scan=2", mz_array, intensity_array, 404.4, 2, 1, "scan=1"),),
    ]
]


def get_scan_data():
    return scan_data


def test_example():
    from psims.mzml.writer import MzMLWriter

    scans = get_scan_data()
    fileName = "/Users/cranneyc/Desktop/systemTestExperimentation/write.mzML"

    with MzMLWriter(open(fileName, "wb")) as out:
        # Add default controlled vocabularies
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
            spectrum_count = len(scans) + sum([len(products) for _, products in scans])
            with out.spectrum_list(count=spectrum_count):
                for scan, products in scans:
                    # Write Precursor scan
                    out.write_spectrum(
                        scan.mz_array,
                        scan.intensity_array,
                        id=scan.id,
                        params=[
                            "MS1 Spectrum",
                            {"ms level": 1},
                            {"total ion current": sum(scan.intensity_array)},
                        ],
                    )
                    # Write MSn scans
                    for prod in products:
                        out.write_spectrum(
                            prod.mz_array,
                            prod.intensity_array,
                            id=prod.id,
                            params=[
                                "MSn Spectrum",
                                {"ms level": 2},
                                {"total ion current": sum(prod.intensity_array)},
                            ],
                            # Include precursor information
                            precursor_information={
                                "mz": prod.precursor_mz,
                                "intensity": prod.precursor_intensity,
                                "charge": prod.precursor_charge,
                                "scan_id": prod.precursor_scan_id,
                                "activation": [
                                    "beam-type collisional dissociation",
                                    {"collision energy": 25},
                                ],
                                "isolation_window": [
                                    prod.precursor_mz - 1,
                                    prod.precursor_mz,
                                    prod.precursor_mz + 1,
                                ],
                            },
                        )


if __name__ == "__main__":
    test_example()
