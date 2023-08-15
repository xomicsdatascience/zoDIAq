import pandas as pd
import numpy as np
import pytest
from csodiaq.targetedReanalysis.targetedReanalysisFunctions import (
    calculate_mz_of_heavy_version_of_peptide,
    filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue,
    filter_to_only_keep_top_peptides_unique_to_protein,
    calculate_mz_of_heavy_isotope_of_each_peptide,
    make_bin_assignments_for_mz_values,
    calculate_binning_information_by_compensation_voltage,
    create_targeted_reanalysis_dataframe,
    organize_for_targeted_reanalysis_of_identified_peptides,
    filter_out_peptides_based_on_user_settings,
    create_targeted_reanalysis_dataframes_by_compensation_voltage,
    create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides,
)


def test__output_formatting_functions__calculate_mz_of_heavy_version_of_peptide():
    testMz = 0.0
    charge = 1
    numLys = 1
    numArg = 1
    peptide = numLys * "K" + numArg * "R"
    expectedOneChargeOneLysOneArg = 18.022469
    oneChargeOneLysOneArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedOneChargeOneLysOneArg == oneChargeOneLysOneArg

    charge = 2
    expectedTwoChargeOneLysOneArg = expectedOneChargeOneLysOneArg / 2
    twoChargeOneLysOneArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedTwoChargeOneLysOneArg == twoChargeOneLysOneArg

    numLys = 2
    charge = 1
    peptide = numLys * "K" + numArg * "R"
    expectedOneChargeTwoLysOneArg = 26.036668
    oneChargeTwoLysOneArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedOneChargeTwoLysOneArg == oneChargeTwoLysOneArg

    charge = 2
    expectedTwoChargeTwoLysOneArg = expectedOneChargeTwoLysOneArg / 2
    twoChargeTwoLysOneArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedTwoChargeTwoLysOneArg == twoChargeTwoLysOneArg

    numLys = 1
    numArg = 2
    charge = 1
    peptide = numLys * "K" + numArg * "R"
    expectedOneChargeOneLysTwoArg = 28.030738999999997
    oneChargeOneLysTwoArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedOneChargeOneLysTwoArg == oneChargeOneLysTwoArg

    charge = 2
    expectedTwoChargeOneLysTwoArg = expectedOneChargeOneLysTwoArg / 2
    twoChargeOneLysTwoArg = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedTwoChargeOneLysTwoArg == twoChargeOneLysTwoArg

    testMz = 100.0
    expectedTwoChargeOneLysTwoArgHundredMz = expectedTwoChargeOneLysTwoArg + testMz
    twoChargeOneLysTwoArgHundredMz = calculate_mz_of_heavy_version_of_peptide(
        peptide, testMz, z=charge
    )
    assert expectedTwoChargeOneLysTwoArgHundredMz == twoChargeOneLysTwoArgHundredMz


def test__output_formatting_functions__filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue():
    data = [
        ["A"],
        ["B"],
        ["C"],
        ["K"],
        ["K"],
        ["R"],
        ["R"],
        ["R"],
    ]
    df = pd.DataFrame(data, columns=["peptide"])
    expectedData = data[3:]
    expectedOutput = pd.DataFrame(expectedData, columns=["peptide"])
    output = filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(
        df
    )
    assert expectedOutput.equals(output)
    pass


def test__output_formatting_functions__filter_to_only_keep_top_peptides_unique_to_protein():
    inputData = [
        ["protein1", 100.0, 1],
        ["protein1", 200.0, 1],
        ["protein1", 300.0, 1],
        ["protein1", 400.0, 1],
        ["protein2", 100.0, 1],
        ["protein2", 200.0, 1],
        ["protein2", 300.0, 1],
        ["protein3", 100.0, 1],
        ["protein3", 200.0, 1],
        ["protein4", 100.0, 1],
        ["protein5", 100.0, 0],
    ]
    inputDf = pd.DataFrame(
        inputData, columns=["leadingProtein", "ionCount", "uniquePeptide"]
    )
    topProteinsToKeep = 2
    expectedOutputData = [
        ["protein1", 400.0, 1],
        ["protein1", 300.0, 1],
        ["protein2", 300.0, 1],
        ["protein2", 200.0, 1],
        ["protein3", 200.0, 1],
        ["protein3", 100.0, 1],
        ["protein4", 100.0, 1],
    ]
    expectedOutputDf = pd.DataFrame(
        expectedOutputData, columns=["leadingProtein", "ionCount", "uniquePeptide"]
    )
    outputDf = filter_to_only_keep_top_peptides_unique_to_protein(
        inputDf, topProteinsToKeep
    )
    assert expectedOutputDf.equals(outputDf)


@pytest.fixture
def inputFilteringDf():
    inputData = [
        ["peptide01R", "protein1", 100.0, 1],
        ["peptide02K", "protein1", 200.0, 1],
        ["peptide03R", "protein1", 300.0, 1],
        ["peptide04K", "protein1", 400.0, 1],
        ["peptide05R", "protein2", 100.0, 1],
        ["peptide06K", "protein2", 200.0, 1],
        ["peptide07R", "protein2", 300.0, 1],
        ["peptide08K", "protein3", 100.0, 1],
        ["peptide09R", "protein3", 200.0, 1],
        ["peptide10K", "protein4", 100.0, 1],
        ["peptide11R", "protein5", 100.0, 0],
        ["peptide12K", "protein6", 100.0, 0],
        ["peptide13", "protein7", 100.0, 1],
    ]
    return pd.DataFrame(
        inputData, columns=["peptide", "leadingProtein", "ionCount", "uniquePeptide"]
    )


def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__no_heavy_no_proteins(
    inputFilteringDf,
):
    isIncludeHeavyIsotopes = False
    maximumPeptidesPerProtein = 0
    expectedOutputDf = inputFilteringDf.copy()
    outputDf = filter_out_peptides_based_on_user_settings(
        inputFilteringDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__is_heavy_no_proteins(
    inputFilteringDf,
):
    isIncludeHeavyIsotopes = True
    maximumPeptidesPerProtein = 0

    expectedOutputData = [
        ["peptide01R", "protein1", 100.0, 1],
        ["peptide02K", "protein1", 200.0, 1],
        ["peptide03R", "protein1", 300.0, 1],
        ["peptide04K", "protein1", 400.0, 1],
        ["peptide05R", "protein2", 100.0, 1],
        ["peptide06K", "protein2", 200.0, 1],
        ["peptide07R", "protein2", 300.0, 1],
        ["peptide08K", "protein3", 100.0, 1],
        ["peptide09R", "protein3", 200.0, 1],
        ["peptide10K", "protein4", 100.0, 1],
        ["peptide11R", "protein5", 100.0, 0],
        ["peptide12K", "protein6", 100.0, 0],
    ]
    expectedOutputDf = pd.DataFrame(
        expectedOutputData, columns=inputFilteringDf.columns
    )
    outputDf = filter_out_peptides_based_on_user_settings(
        inputFilteringDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__no_heavy_has_proteins(
    inputFilteringDf,
):
    isIncludeHeavyIsotopes = False
    maximumPeptidesPerProtein = 2
    expectedOutputData = [
        ["peptide04K", "protein1", 400.0, 1],
        ["peptide03R", "protein1", 300.0, 1],
        ["peptide07R", "protein2", 300.0, 1],
        ["peptide06K", "protein2", 200.0, 1],
        ["peptide09R", "protein3", 200.0, 1],
        ["peptide08K", "protein3", 100.0, 1],
        ["peptide10K", "protein4", 100.0, 1],
        ["peptide13", "protein7", 100.0, 1],
    ]
    expectedOutputDf = pd.DataFrame(
        expectedOutputData, columns=inputFilteringDf.columns
    )
    outputDf = filter_out_peptides_based_on_user_settings(
        inputFilteringDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__is_heavy_has_proteins(
    inputFilteringDf,
):
    isIncludeHeavyIsotopes = True
    maximumPeptidesPerProtein = 2
    expectedOutputData = [
        ["peptide04K", "protein1", 400.0, 1],
        ["peptide03R", "protein1", 300.0, 1],
        ["peptide07R", "protein2", 300.0, 1],
        ["peptide06K", "protein2", 200.0, 1],
        ["peptide09R", "protein3", 200.0, 1],
        ["peptide08K", "protein3", 100.0, 1],
        ["peptide10K", "protein4", 100.0, 1],
    ]
    expectedOutputDf = pd.DataFrame(
        expectedOutputData, columns=inputFilteringDf.columns
    )
    outputDf = filter_out_peptides_based_on_user_settings(
        inputFilteringDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__calculate_mz_of_heavy_isotope_of_each_peptide():
    data = [
        ["KR", 100.0, 1],
        ["KKR", 100.0, 1],
        ["KRR", 100.0, 1],
    ]
    inputDf = pd.DataFrame(data, columns=["peptide", "MzLIB", "zLIB"])
    expectedOutput = [118.022469, 126.036668, 128.030738999999997]
    output = calculate_mz_of_heavy_isotope_of_each_peptide(inputDf)
    np.testing.assert_array_almost_equal(np.array(expectedOutput), np.array(output))


@pytest.fixture
def inputBinningDf():
    inputData = [
        ["R", 100.0, 1],
        ["R", 100.2, 1],
        ["K", 100.49, 2],
        ["K", 100.5, 2],
        ["R", 102.0, 1],
        ["K", 115.0, 1],
        ["RK", 115.1, 1],
    ]
    return pd.DataFrame(inputData, columns=["peptide", "MzLIB", "zLIB"])


@pytest.fixture
def expectedLightMzBins():
    return np.array(
        [
            99.75,
            99.75,
            99.75,
            101.25,
            102.75,
            114.75,
            114.75,
        ]
    )


def test__output_formatting_functions__make_bin_assignments_for_mz_values(
    inputBinningDf, expectedLightMzBins
):
    binValueProximity = 0.75
    bins = make_bin_assignments_for_mz_values(
        inputBinningDf["MzLIB"], maxDistanceFromBinValue=binValueProximity
    )
    np.testing.assert_array_equal(expectedLightMzBins, bins)


def test__output_formatting_functions__make_bin_assignments_for_mz_values__custom_bin_value_works(
    inputBinningDf,
):
    expectedLightMzBinsForBinValues1 = [
        100.0,
        100.0,
        100.0,
        100.0,
        102.0,
        116.0,
        116.0,
    ]
    binValueProximity = 1.0
    bins = make_bin_assignments_for_mz_values(
        inputBinningDf["MzLIB"], maxDistanceFromBinValue=binValueProximity
    )
    np.testing.assert_array_equal(expectedLightMzBinsForBinValues1, bins)


def test__output_formatting_functions__organize_for_targeted_reanalysis_of_identified_peptides__no_heavy(
    inputBinningDf, expectedLightMzBins
):
    isIncludeHeavyIsotopes = False
    binValueProximity = 0.75
    expectedOutputDf = inputBinningDf.copy()
    expectedOutputDf["lightMzBin"] = expectedLightMzBins
    outputDf = organize_for_targeted_reanalysis_of_identified_peptides(
        inputBinningDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        binValueProximity=binValueProximity,
    )
    assert expectedOutputDf.equals(outputDf)


@pytest.fixture
def expectedHeavyMzColumn(inputBinningDf):
    lightAndHeavyLysKMassDiff = 8.014199
    lightAndHeavyArgRMassDiff = 10.00827
    mzValues = inputBinningDf["MzLIB"]
    chargeValues = inputBinningDf["zLIB"]
    return [
        mzValues[0] + lightAndHeavyArgRMassDiff / chargeValues[0],
        mzValues[1] + lightAndHeavyArgRMassDiff / chargeValues[1],
        mzValues[2] + lightAndHeavyLysKMassDiff / chargeValues[2],
        mzValues[3] + lightAndHeavyLysKMassDiff / chargeValues[3],
        mzValues[4] + lightAndHeavyArgRMassDiff / chargeValues[4],
        mzValues[5] + lightAndHeavyLysKMassDiff / chargeValues[5],
        mzValues[6]
        + lightAndHeavyArgRMassDiff / chargeValues[6]
        + lightAndHeavyLysKMassDiff / chargeValues[6],
    ]


@pytest.fixture
def expectedHeavyMzBinColumn():
    return np.array(
        [
            109.75,
            109.75,
            103.75,
            105.25,
            112.75,
            123.25,
            133.75,
        ]
    )


def test__output_formatting_functions__organize_for_targeted_reanalysis_of_identified_peptides__has_heavy(
    inputBinningDf, expectedLightMzBins, expectedHeavyMzColumn, expectedHeavyMzBinColumn
):
    isIncludeHeavyIsotopes = True
    binValueProximity = 0.75
    expectedOutputDf = inputBinningDf.copy()
    expectedOutputDf["lightMzBin"] = expectedLightMzBins
    expectedOutputDf["heavyMz"] = expectedHeavyMzColumn
    expectedOutputDf["heavyMzBin"] = expectedHeavyMzBinColumn
    outputDf = organize_for_targeted_reanalysis_of_identified_peptides(
        inputBinningDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        binValueProximity=binValueProximity,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__calculate_binning_information_by_compensation_voltage__no_heavy(
    inputBinningDf, expectedLightMzBins
):
    isIncludeHeavyIsotopes = False
    binValueProximity = 0.75
    inputBinningDfCV30 = inputBinningDf.copy()
    inputBinningDfCV30["CompensationVoltage"] = [-30] * len(inputBinningDfCV30.index)
    inputBinningDfCV40 = inputBinningDf.copy()
    inputBinningDfCV40["CompensationVoltage"] = [-40] * len(inputBinningDfCV40.index)
    inputDf = pd.concat([inputBinningDfCV40, inputBinningDfCV30])
    expectedOutputDf = inputDf.copy()
    expectedOutputDf["lightMzBin"] = np.append(expectedLightMzBins, expectedLightMzBins)
    outputDf = calculate_binning_information_by_compensation_voltage(
        inputDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        binValueProximity=binValueProximity,
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__calculate_binning_information_by_compensation_voltage__with_heavy(
    inputBinningDf, expectedLightMzBins, expectedHeavyMzColumn, expectedHeavyMzBinColumn
):
    binValueProximity = 0.75
    isIncludeHeavyIsotopes = True
    inputBinningDfCV30 = inputBinningDf.copy()
    inputBinningDfCV30["CompensationVoltage"] = [-30] * len(inputBinningDfCV30.index)
    inputBinningDfCV40 = inputBinningDf.copy()
    inputBinningDfCV40["CompensationVoltage"] = [-40] * len(inputBinningDfCV40.index)
    inputDf = pd.concat([inputBinningDfCV40, inputBinningDfCV30])
    expectedOutputDf = inputDf.copy()
    expectedOutputDf["lightMzBin"] = np.append(expectedLightMzBins, expectedLightMzBins)
    expectedOutputDf["heavyMz"] = expectedHeavyMzColumn + expectedHeavyMzColumn
    expectedOutputDf["heavyMzBin"] = np.append(
        expectedHeavyMzBinColumn, expectedHeavyMzBinColumn
    )
    outputDf = calculate_binning_information_by_compensation_voltage(
        inputDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        binValueProximity=binValueProximity,
    )
    assert expectedOutputDf.equals(outputDf)


@pytest.fixture
def inputFormattedDf():
    inputData = [
        ["peptide1", 20.0, 30.0],
        ["peptide2", 20.0, 30.0],
        ["peptide3", 10.0, 20.0],
    ]
    return pd.DataFrame(inputData, columns=["peptide", "lightMzBin", "heavyMzBin"])


@pytest.fixture
def targetedReanalysisNoHeavyDf():
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    expectedOutputData = [
        ["1/peptide3", formula, adduct, 10.0, charge, 1],
        ["2/peptide1/peptide2", formula, adduct, 20.0, charge, 2],
    ]
    return pd.DataFrame(
        expectedOutputData,
        columns=["Compound", "Formula", "Adduct", "m.z", "z", "MSXID"],
    )


def test__output_formatting_functions__create_targeted_reanalysis_dataframe__no_heavy(
    inputFormattedDf, targetedReanalysisNoHeavyDf
):
    isIncludeHeavyIsotopes = False
    outputDf = create_targeted_reanalysis_dataframe(
        inputFormattedDf, isIncludeHeavyIsotopes=isIncludeHeavyIsotopes
    )
    assert targetedReanalysisNoHeavyDf.equals(outputDf)


@pytest.fixture
def targetedReanalysisWithHeavyDf():
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    expectedOutputData = [
        ["1/peptide3", formula, adduct, 10.0, charge, 1],
        ["1/peptide3", formula, adduct, 20.0, charge, 1],
        ["2/peptide1/peptide2", formula, adduct, 20.0, charge, 2],
        ["2/peptide1/peptide2", formula, adduct, 30.0, charge, 2],
    ]
    return pd.DataFrame(
        expectedOutputData,
        columns=["Compound", "Formula", "Adduct", "m.z", "z", "MSXID"],
    )


def test__output_formatting_functions__create_targeted_reanalysis_dataframe__with_heavy(
    inputFormattedDf, targetedReanalysisWithHeavyDf
):
    isIncludeHeavyIsotopes = True
    outputDf = create_targeted_reanalysis_dataframe(
        inputFormattedDf, isIncludeHeavyIsotopes=isIncludeHeavyIsotopes
    )
    assert targetedReanalysisWithHeavyDf.equals(outputDf)


def test__output_formatting_functions__create_targeted_reanalysis_dataframes_by_compensation_voltage__no_heavy(
    inputFormattedDf, targetedReanalysisNoHeavyDf
):
    isIncludeHeavyIsotopes = False
    inputFormattedDfCV30 = inputFormattedDf.copy()
    inputFormattedDfCV30["CompensationVoltage"] = [-30] * len(
        inputFormattedDfCV30.index
    )
    inputFormattedDfCV40 = inputFormattedDf.copy()
    inputFormattedDfCV40["CompensationVoltage"] = [-40] * len(
        inputFormattedDfCV40.index
    )
    inputDf = pd.concat([inputFormattedDfCV30, inputFormattedDfCV40])
    expectedOutput = {
        "CV_30": targetedReanalysisNoHeavyDf,
        "CV_40": targetedReanalysisNoHeavyDf,
    }
    output = create_targeted_reanalysis_dataframes_by_compensation_voltage(
        inputDf, isIncludeHeavyIsotopes
    )
    for cv, expectedTargetedReanalysisDf in expectedOutput.items():
        assert cv in output
        assert expectedTargetedReanalysisDf.equals(output[cv])


def test__output_formatting_functions__create_targeted_reanalysis_dataframes_by_compensation_voltage__with_heavy(
    inputFormattedDf, targetedReanalysisWithHeavyDf
):
    isIncludeHeavyIsotopes = True
    inputFormattedDfCV30 = inputFormattedDf.copy()
    inputFormattedDfCV30["CompensationVoltage"] = [-30] * len(
        inputFormattedDfCV30.index
    )
    inputFormattedDfCV40 = inputFormattedDf.copy()
    inputFormattedDfCV40["CompensationVoltage"] = [-40] * len(
        inputFormattedDfCV40.index
    )
    inputDf = pd.concat([inputFormattedDfCV30, inputFormattedDfCV40])
    expectedOutput = {
        "CV_30": targetedReanalysisWithHeavyDf,
        "CV_40": targetedReanalysisWithHeavyDf,
    }
    output = create_targeted_reanalysis_dataframes_by_compensation_voltage(
        inputDf, isIncludeHeavyIsotopes
    )
    for cv, expectedTargetedReanalysisDf in expectedOutput.items():
        assert cv in output
        assert expectedTargetedReanalysisDf.equals(output[cv])


@pytest.fixture
def inputProteinFdrDf():
    genericMz = 100.0
    inputData = [
        ["peptide01R", "protein1", 100.0, 1, genericMz, 1, ""],
        ["peptide02R", "protein1", 200.0, 1, genericMz * 2, 2, ""],
        ["peptide03K", "protein1", 300.0, 1, genericMz * 3, 1, ""],
        ["peptide04K", "protein1", 400.0, 1, genericMz * 4, 1, ""],
        ["peptide05R", "protein2", 100.0, 1, genericMz * 5, 1, ""],
        ["peptide06R", "protein2", 200.0, 1, genericMz, 2, ""],
        ["peptide07K", "protein2", 300.0, 1, genericMz * 2, 1, ""],
        ["peptide08KK", "protein3", 100.0, 1, genericMz * 3, 2, ""],
        ["peptide09R", "protein3", 200.0, 1, genericMz * 4, 1, ""],
        ["peptide10RR", "protein4", 100.0, 1, genericMz * 5, 1, ""],
        ["peptide11K", "protein5", 100.0, 0, genericMz, 1, ""],
        ["peptide12K", "protein6", 100.0, 0, genericMz * 2, 2, ""],
        ["peptide13", "protein7", 100.0, 1, genericMz * 3, 1, ""],
    ]
    return pd.DataFrame(
        inputData,
        columns=[
            "peptide",
            "leadingProtein",
            "ionCount",
            "uniquePeptide",
            "MzLIB",
            "zLIB",
            "CompensationVoltage",
        ],
    )


def test__output_formatting_functions__create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides__no_heavy_no_protein_no_cv(
    inputProteinFdrDf,
):
    isIncludeHeavyIsotopes = False
    maximumPeptidesPerProtein = 0
    binValueProximity = 0.75
    fullDf = inputProteinFdrDf.copy()
    fullDf["lightMzBin"] = [
        99.75,
        200.25,
        300.75,
        399.75,
        500.25,
        99.75,
        200.25,
        300.75,
        399.75,
        500.25,
        99.75,
        200.25,
        300.75,
    ]
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    targetedReanalysisData = [
        ["3/peptide01R/peptide06R/peptide11K", formula, adduct, 99.75, charge, 1],
        ["3/peptide02R/peptide07K/peptide12K", formula, adduct, 200.25, charge, 2],
        ["3/peptide03K/peptide08KK/peptide13", formula, adduct, 300.75, charge, 3],
        ["2/peptide04K/peptide09R", formula, adduct, 399.75, charge, 4],
        ["2/peptide05R/peptide10RR", formula, adduct, 500.25, charge, 5],
    ]
    targetedReanalysisDf = pd.DataFrame(
        targetedReanalysisData,
        columns=["Compound", "Formula", "Adduct", "m.z", "z", "MSXID"],
    )
    expectedOutputDict = {
        "noCV": targetedReanalysisDf,
        "fullDf": fullDf,
    }
    outputDict = create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides(
        inputProteinFdrDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
        binValueProximity=binValueProximity,
    )
    for type, df in expectedOutputDict.items():
        assert type in outputDict
        assert df.equals(outputDict[type])


def test__output_formatting_functions__create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides__with_heavy_with_protein_with_cv(
    inputProteinFdrDf,
):
    isIncludeHeavyIsotopes = True
    maximumPeptidesPerProtein = 2
    binValueProximity = 0.75
    filteredSortedIdxs = [3, 2, 6, 5, 8, 7, 9]
    filteredDf = inputProteinFdrDf.loc[filteredSortedIdxs].reset_index(drop=True)
    inputDfCV30 = inputProteinFdrDf.copy()
    inputDfCV30["CompensationVoltage"] = [-30] * len(inputDfCV30.index)
    inputDfCV40 = inputProteinFdrDf.copy()
    inputDfCV40["CompensationVoltage"] = [-40] * len(inputDfCV40.index)
    inputDfCV40["leadingProtein"] = [
        "protein08",
        "protein08",
        "protein08",
        "protein08",
        "protein09",
        "protein09",
        "protein09",
        "protein10",
        "protein10",
        "protein11",
        "protein12",
        "protein13",
        "protein14",
    ]
    inputDf = pd.concat([inputDfCV40, inputDfCV30])
    filteredDf["lightMzBin"] = [
        399.75,
        300.75,
        200.25,
        99.75,
        399.75,
        300.75,
        500.25,
    ]
    filteredDf["heavyMz"] = [
        408.014199,
        308.014199,
        208.014199,
        105.004135,
        410.008270,
        308.014199,
        520.016540,
    ]
    filteredDf["heavyMzBin"] = [
        407.75,
        308.75,
        208.25,
        104.75,
        410.75,
        308.75,
        520.25,
    ]
    filteredDf["CompensationVoltage"] = [-30] * len(filteredDf.index)
    secondCvFilteredDf = filteredDf.copy()
    secondCvFilteredDf["leadingProtein"] = list(
        inputDfCV40.loc[filteredSortedIdxs]["leadingProtein"]
    )
    secondCvFilteredDf["CompensationVoltage"] = [-40] * len(secondCvFilteredDf.index)

    filteredDf = pd.concat([secondCvFilteredDf, filteredDf])
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    targetedReanalysisData = [
        ["1/peptide06R", formula, adduct, 99.75, charge, 1],
        ["1/peptide06R", formula, adduct, 104.75, charge, 1],
        ["1/peptide07K", formula, adduct, 200.25, charge, 2],
        ["1/peptide07K", formula, adduct, 208.25, charge, 2],
        ["2/peptide03K/peptide08KK", formula, adduct, 300.75, charge, 3],
        ["2/peptide03K/peptide08KK", formula, adduct, 308.75, charge, 3],
        ["1/peptide04K", formula, adduct, 399.75, charge, 4],
        ["1/peptide04K", formula, adduct, 407.75, charge, 4],
        ["1/peptide09R", formula, adduct, 399.75, charge, 5],
        ["1/peptide09R", formula, adduct, 410.75, charge, 5],
        ["1/peptide10RR", formula, adduct, 500.25, charge, 6],
        ["1/peptide10RR", formula, adduct, 520.25, charge, 6],
    ]
    targetedReanalysisDf = pd.DataFrame(
        targetedReanalysisData,
        columns=["Compound", "Formula", "Adduct", "m.z", "z", "MSXID"],
    )
    expectedOutputDict = {
        "CV_30": targetedReanalysisDf,
        "CV_40": targetedReanalysisDf,
        "fullDf": filteredDf,
    }
    outputDict = create_mass_spec_input_dataframes_for_targeted_reanalysis_of_identified_peptides(
        inputDf,
        isIncludeHeavyIsotopes=isIncludeHeavyIsotopes,
        maximumPeptidesPerProtein=maximumPeptidesPerProtein,
        binValueProximity=binValueProximity,
    )
    expectedFullDf = (
        expectedOutputDict["fullDf"]
        .sort_values(["peptide", "leadingProtein", "MzLIB", "CompensationVoltage"])
        .reset_index(drop=True)
    )
    fullDf = (
        outputDict["fullDf"]
        .sort_values(["peptide", "leadingProtein", "MzLIB", "CompensationVoltage"])
        .reset_index(drop=True)
    )
    assert expectedFullDf.equals(fullDf)
    for type, df in expectedOutputDict.items():
        if type == "fullDf":
            continue
        assert type in outputDict
        assert df.equals(outputDict[type])
