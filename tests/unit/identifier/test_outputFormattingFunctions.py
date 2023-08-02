from csodiaq.identifier.outputFormattingFunctions import (
    format_output_line,
    extract_metadata_from_match_and_score_dataframes,
    format_output_as_pandas_dataframe,
    drop_duplicate_values_from_df_in_given_column,
    create_spectral_fdr_output_from_full_output,
    create_peptide_fdr_output_from_full_output,
    identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff,
    organize_peptide_df_by_leading_proteins,
    determine_if_peptides_are_unique_to_leading_protein,
    calculate_mz_of_heavy_version_of_peptide,
    filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue,
    filter_to_only_keep_top_peptides_unique_to_protein,
    calculate_mz_of_heavy_isotope_of_each_peptide,
    make_bin_assignments_for_mz_values,
    calculate_binning_information_by_compensation_voltage,
    create_targeted_reanalysis_dataframe,
    organize_for_targeted_reanalysis_of_identified_peptides,
    filter_out_peptides_based_on_user_settings,
    create_targeted_reanalysis_dataframe_by_compensation_voltage__no_heavy,
)
import pandas as pd
import pytest
import numpy as np

@pytest.fixture
def identifierOutputData():
    return [
        "3",
        4,
        "testPeptide",
        "testProtein",
        0,
        0,
        1,
        8,
        "testIdentifiers",
        5,
        2,
        9,
        10,
        6,
        7,
        11,
        12,
    ]


def test__output_formatting_functions__format_output_line(identifierOutputData):
    libDict = {
        "peptide": "testPeptide",
        "proteinName": "testProtein",
        "isDecoy": 0,
        "precursorMz": 0,
        "precursorCharge": 1,
        "identifier": "testIdentifiers",
        "peaks": list(range(2)),
    }
    queryDict = {
        "scan": "3",
        "precursorMz": 4,
        "peaksCount": 5,
        "CV": 6,
        "windowWidth": 7,
    }
    matchDict = {
        "cosineSimilarityScore": 8,
        "shared": 9,
        "ionCount": 10,
        "maccScore": 11,
        "exclude_num": 12,
    }
    output = format_output_line(libDict, queryDict, matchDict)
    assert output == identifierOutputData


def test__output_formatting_functions__extract_metadata_from_match_and_score_dataframes():
    lib1Idx = 0
    lib2Idx = 1
    queryIdx = 0
    genericIntensity = 100.0
    genericPpmDifference = 10.0
    lib1PeakCount = 10
    lib1PeakCountAbovePrecursorMz = 3
    lib2PeakCount = 4
    precursorMz = 100.0
    abovePrecursorMz = precursorMz + 10
    belowPrecursorMz = precursorMz - 10
    lib1CosineScore = 1.0
    lib2CosineScore = 0.9
    lib1Score = 0.8
    lib2Score = 0.7
    lib1MatchAbovePrecursorMz = [
        [
            lib1Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            abovePrecursorMz,
            genericPpmDifference,
        ]
    ] * lib1PeakCountAbovePrecursorMz
    lib1MatchBelowPrecursorMz = [
        [
            lib1Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            belowPrecursorMz,
            genericPpmDifference,
        ]
    ] * (lib1PeakCount - lib1PeakCountAbovePrecursorMz)
    lib2Match = [
        [
            lib2Idx,
            genericIntensity,
            queryIdx,
            genericIntensity,
            abovePrecursorMz,
            genericPpmDifference,
        ]
    ] * lib2PeakCount
    columns = [
        "libraryIdx",
        "libraryIntensity",
        "queryIdx",
        "queryIntensity",
        "queryMz",
        "ppmDifference",
    ]
    matchDf = pd.DataFrame(
        lib1MatchAbovePrecursorMz + lib1MatchBelowPrecursorMz + lib2Match,
        columns=columns,
    )

    scoreData = [
        [lib1Idx, queryIdx, lib1CosineScore, lib1Score],
        [lib2Idx, queryIdx, lib2CosineScore, lib2Score],
    ]
    scoreDf = pd.DataFrame(
        scoreData, columns=["libraryIdx", "queryIdx", "cosineScore", "maccScore"]
    )
    expectedOutput = {
        (lib1Idx, queryIdx): {
            "cosineSimilarityScore": lib1CosineScore,
            "shared": lib1PeakCount,
            "ionCount": genericIntensity * lib1PeakCountAbovePrecursorMz,
            "maccScore": lib1Score,
            "exclude_num": lib1PeakCount - lib1PeakCountAbovePrecursorMz,
        },
        (lib2Idx, queryIdx): {
            "cosineSimilarityScore": lib2CosineScore,
            "shared": lib2PeakCount,
            "ionCount": lib2PeakCount * genericIntensity,
            "maccScore": lib2Score,
            "exclude_num": 0,
        },
    }
    queryDict = {
        str(queryIdx): {
            "precursorMz": precursorMz,
        },
    }
    output = extract_metadata_from_match_and_score_dataframes(
        matchDf, scoreDf, queryDict
    )
    for idKey, metadataDict in expectedOutput.items():
        assert idKey in output
        for metadataType, metadataValue in metadataDict.items():
            assert metadataType in output[idKey]
            assert output[idKey][metadataType] == metadataValue


def test__output_formatting_functions__format_output_as_pandas_dataframe(
    identifierOutputData,
):
    expectedColumns = [
        "fileName",
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
        "MaCC_Score",
        "exclude_num",
    ]
    inputFileName = "dummyFile"
    expectedOutputDf = pd.DataFrame(
        [[inputFileName] + identifierOutputData], columns=expectedColumns
    )
    outputDf = format_output_as_pandas_dataframe(inputFileName, [identifierOutputData])
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__drop_duplicate_values_from_df_in_given_column():
    columnName = "test"
    data = [
        [0],
        [0],
        [1],
        [1],
        [2],
    ]
    df = pd.DataFrame(data, columns=[columnName])
    expectedOutputData = [
        [0],
        [1],
        [2],
    ]
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=[columnName])
    outputDf = drop_duplicate_values_from_df_in_given_column(df, columnName)
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__create_spectral_fdr_output_from_full_output():
    numNonDecoys = 100
    numDecoys = 2
    isDecoyColumn = ([0] * numNonDecoys) + ([1] * numDecoys)
    scoreColumn = list(range(len(isDecoyColumn) - 1, -1, -1))
    inputDf = pd.DataFrame()
    inputDf["isDecoy"] = isDecoyColumn

    expectedOutputDf = pd.DataFrame()
    expectedOutputDf["isDecoy"] = isDecoyColumn[:-1]
    expectedOutputDf["spectralFDR"] = [0] * numNonDecoys + [1 / (numNonDecoys + 1)]

    outputDf = create_spectral_fdr_output_from_full_output(inputDf)
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__create_peptide_fdr_output_from_full_output():
    numDuplicatePeptides = 2
    numNonDuplicatePeptides = 99
    numNonDecoys = numDuplicatePeptides + numNonDuplicatePeptides
    numDecoys = 2
    peptideColumn = (
        [0] * numDuplicatePeptides
        + list(range(1, numNonDuplicatePeptides + 1))
        + [f"decoy{i}" for i in range(numDecoys)]
    )
    isDecoyColumn = ([0] * numNonDecoys) + ([1] * numDecoys)
    spectralFDRColumn = (
        [0] * numNonDecoys + [1 / (numNonDecoys + 1)] + [2 / (numNonDecoys + 2)]
    )
    inputDf = pd.DataFrame()
    inputDf["peptide"] = peptideColumn
    inputDf["isDecoy"] = isDecoyColumn

    expectedOutputDf = pd.DataFrame()
    expectedOutputDf["peptide"] = peptideColumn[1:-1]
    expectedOutputDf["isDecoy"] = isDecoyColumn[1:-1]
    expectedOutputDf["peptideFDR"] = [0] * (numNonDecoys - 1) + [1 / numNonDecoys]
    outputDf = create_peptide_fdr_output_from_full_output(inputDf)

    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__organize_peptide_df_by_leading_proteins():
    peptideProteinData = [
        ["peptide01", "1/protein7"],
        ["peptide02", "3/protein4/protein6/protein9"],
        ["peptide03", "1/protein1"],
        ["peptide04", "2/protein1/protein5"],
        ["peptide05", "1/protein7"],
        ["peptide06", "2/protein3/protein6"],
        ["peptide07", "1/protein1"],
        ["peptide08", "4/protein1/protein2/protein5/protein8"],
        ["peptide09", "1/protein1"],
        ["peptide10", "2/protein4/protein9"],
    ]
    peptideProteinDf = pd.DataFrame(peptideProteinData, columns=["peptide", "protein"])
    leadingProteins = set(
        [
            ("protein1",),
            ("protein4", "protein9"),
            ("protein6",),
            ("protein7",),
        ]
    )
    expectedOutputData = [
        ["peptide01", "1/protein7", "1/protein7"],
        ["peptide02", "3/protein4/protein6/protein9", "2/protein4/protein9"],
        ["peptide02", "3/protein4/protein6/protein9", "1/protein6"],
        ["peptide03", "1/protein1", "1/protein1"],
        ["peptide04", "2/protein1/protein5", "1/protein1"],
        ["peptide05", "1/protein7", "1/protein7"],
        ["peptide06", "2/protein3/protein6", "1/protein6"],
        ["peptide07", "1/protein1", "1/protein1"],
        ["peptide08", "4/protein1/protein2/protein5/protein8", "1/protein1"],
        ["peptide09", "1/protein1", "1/protein1"],
        ["peptide10", "2/protein4/protein9", "2/protein4/protein9"],
    ]
    expectedOutputDf = pd.DataFrame(
        expectedOutputData, columns=["peptide", "protein", "leadingProtein"]
    )
    outputDf = organize_peptide_df_by_leading_proteins(
        peptideProteinDf, leadingProteins
    )
    assert expectedOutputDf.equals(outputDf)


def test__output_formatting_functions__identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff():
    numLeadingProteins = 100
    duplicateLeadingProtein = 0
    decoyLeadingProteins = ["decoy1", "decoy2"]
    leadingProteinColumn = (
        list(range(numLeadingProteins))
        + [duplicateLeadingProtein]
        + decoyLeadingProteins
    )
    isDecoyColumn = [
        0 for i in range(len(leadingProteinColumn) - len(decoyLeadingProteins))
    ]
    isDecoyColumn += [1 for i in range(len(decoyLeadingProteins))]
    df = pd.DataFrame(
        {
            "leadingProtein": leadingProteinColumn,
            "isDecoy": isDecoyColumn,
        }
    )
    expectedOutput = {i: 0.0 for i in range(numLeadingProteins)}
    expectedOutput["decoy1"] = 1 / (numLeadingProteins + 1)
    output = identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff(
        df
    )
    assert expectedOutput == output


def test__output_formatting_functions__determine_if_peptides_are_unique_to_leading_protein():
    inputData = [
        ["peptide1", "protein1"],
        ["peptide2", "protein1"],
        ["peptide3", "protein2"],
    ]
    inputDf = pd.DataFrame(inputData, columns=["peptide", "leadingProtein"])
    expectedOutput = [
        0,
        0,
        1,
    ]
    output = determine_if_peptides_are_unique_to_leading_protein(inputDf)
    assert expectedOutput == output

def test__output_formatting_functions__calculate_mz_of_heavy_version_of_peptide():
    testMz = 0.0

    charge = 1
    numLys = 1
    numArg = 1
    peptide = numLys * 'K' + numArg * 'R'
    expectedOneChargeOneLysOneArg = 18.022469
    oneChargeOneLysOneArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedOneChargeOneLysOneArg == oneChargeOneLysOneArg

    charge = 2
    expectedTwoChargeOneLysOneArg = expectedOneChargeOneLysOneArg / 2
    twoChargeOneLysOneArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedTwoChargeOneLysOneArg == twoChargeOneLysOneArg

    numLys = 2
    charge = 1
    peptide = numLys * 'K' + numArg * 'R'
    expectedOneChargeTwoLysOneArg = 26.036668
    oneChargeTwoLysOneArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedOneChargeTwoLysOneArg == oneChargeTwoLysOneArg

    charge = 2
    expectedTwoChargeTwoLysOneArg = expectedOneChargeTwoLysOneArg / 2
    twoChargeTwoLysOneArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedTwoChargeTwoLysOneArg == twoChargeTwoLysOneArg

    numLys = 1
    numArg = 2
    charge = 1
    peptide = numLys * 'K' + numArg * 'R'
    expectedOneChargeOneLysTwoArg = 28.030738999999997
    oneChargeOneLysTwoArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedOneChargeOneLysTwoArg == oneChargeOneLysTwoArg

    charge = 2
    expectedTwoChargeOneLysTwoArg = expectedOneChargeOneLysTwoArg / 2
    twoChargeOneLysTwoArg = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
    assert expectedTwoChargeOneLysTwoArg == twoChargeOneLysTwoArg

    testMz = 100.0
    expectedTwoChargeOneLysTwoArgHundredMz = expectedTwoChargeOneLysTwoArg + testMz
    twoChargeOneLysTwoArgHundredMz = calculate_mz_of_heavy_version_of_peptide(peptide, testMz, z=charge)
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
    output = filter_to_only_keep_peptides_with_possibly_heavy_K_or_R_terminal_residue(df)
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
    inputDf = pd.DataFrame(inputData, columns=["leadingProtein","ionCount", "uniquePeptide"])
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
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=["leadingProtein","ionCount", "uniquePeptide"])
    outputDf = filter_to_only_keep_top_peptides_unique_to_protein(inputDf, topProteinsToKeep)
    assert expectedOutputDf.equals(outputDf)

@pytest.fixture
def filteringInputDf():
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
    return pd.DataFrame(inputData, columns=["peptide","leadingProtein","ionCount","uniquePeptide"])

def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__no_heavy_no_proteins(filteringInputDf):
    isIncludeHeavy = False
    expectedOutputDf = filteringInputDf.copy()
    outputDf = filter_out_peptides_based_on_user_settings(filteringInputDf, isIncludeHeavy=isIncludeHeavy)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__is_heavy_no_proteins(filteringInputDf):
    isIncludeHeavy = True
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
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=filteringInputDf.columns)
    outputDf = filter_out_peptides_based_on_user_settings(filteringInputDf, isIncludeHeavy=isIncludeHeavy)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__no_heavy_has_proteins(filteringInputDf):
    isIncludeHeavy = False
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
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=filteringInputDf.columns)
    outputDf = filter_out_peptides_based_on_user_settings(filteringInputDf, isIncludeHeavy=isIncludeHeavy, maximumPeptidesPerProtein=maximumPeptidesPerProtein)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__filter_out_peptides_based_on_user_settings__is_heavy_has_proteins(filteringInputDf):
    isIncludeHeavy = True
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
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=filteringInputDf.columns)
    outputDf = filter_out_peptides_based_on_user_settings(filteringInputDf, isIncludeHeavy=isIncludeHeavy, maximumPeptidesPerProtein=maximumPeptidesPerProtein)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__calculate_mz_of_heavy_isotope_of_each_peptide():
    data = [
        ["KR",100.0, 1],
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
    return np.array([
        99.75,
        99.75,
        99.75,
        101.25,
        102.75,
        114.75,
        114.75,
    ])

def test__output_formatting_functions__make_bin_assignments_for_mz_values(inputBinningDf, expectedLightMzBins):
    bins = make_bin_assignments_for_mz_values(inputBinningDf["MzLIB"])
    np.testing.assert_array_equal(expectedLightMzBins,bins)

def test__output_formatting_functions__organize_for_targeted_reanalysis_of_identified_peptides__no_heavy(inputBinningDf, expectedLightMzBins):
    isIncludeHeavy = False
    expectedOutputDf = inputBinningDf.copy()
    expectedOutputDf["lightMzBin"] = expectedLightMzBins
    outputDf = organize_for_targeted_reanalysis_of_identified_peptides(inputBinningDf, isIncludeHeavy=isIncludeHeavy)
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
        mzValues[6] + lightAndHeavyArgRMassDiff / chargeValues[6] + lightAndHeavyLysKMassDiff / chargeValues[6],
    ]

@pytest.fixture
def expectedHeavyMzBinColumn():
    return np.array([
        109.75,
        109.75,
        103.75,
        105.25,
        112.75,
        123.25,
        133.75,
    ])

def test__output_formatting_functions__organize_for_targeted_reanalysis_of_identified_peptides__has_heavy(inputBinningDf, expectedLightMzBins, expectedHeavyMzColumn, expectedHeavyMzBinColumn):
    isIncludeHeavy = True
    expectedOutputDf = inputBinningDf.copy()
    expectedOutputDf["lightMzBin"] = expectedLightMzBins
    expectedOutputDf["heavyMz"] = expectedHeavyMzColumn
    expectedOutputDf["heavyMzBin"] = expectedHeavyMzBinColumn
    outputDf = organize_for_targeted_reanalysis_of_identified_peptides(inputBinningDf, isIncludeHeavy=isIncludeHeavy)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__calculate_binning_information_by_compensation_voltage__no_heavy(inputBinningDf, expectedLightMzBins):
    isIncludeHeavy = False
    inputBinningDfCV30 = inputBinningDf.copy()
    inputBinningDfCV30["CompensationVoltage"] = [-30] * len(inputBinningDfCV30.index)
    inputBinningDfCV40 = inputBinningDf.copy()
    inputBinningDfCV40["CompensationVoltage"] = [-40] * len(inputBinningDfCV40.index)
    inputDf = pd.concat([inputBinningDfCV40, inputBinningDfCV30])
    expectedOutputDf = inputDf.copy()
    expectedOutputDf["lightMzBin"] = np.append(expectedLightMzBins, expectedLightMzBins)
    outputDf = calculate_binning_information_by_compensation_voltage(inputDf, isIncludeHeavy=isIncludeHeavy)
    assert expectedOutputDf.equals(outputDf)

def test__output_formatting_functions__calculate_binning_information_by_compensation_voltage__with_heavy(inputBinningDf, expectedLightMzBins, expectedHeavyMzColumn, expectedHeavyMzBinColumn):
    isIncludeHeavy = True
    inputBinningDfCV30 = inputBinningDf.copy()
    inputBinningDfCV30["CompensationVoltage"] = [-30] * len(inputBinningDfCV30.index)
    inputBinningDfCV40 = inputBinningDf.copy()
    inputBinningDfCV40["CompensationVoltage"] = [-40] * len(inputBinningDfCV40.index)
    inputDf = pd.concat([inputBinningDfCV40, inputBinningDfCV30])
    expectedOutputDf = inputDf.copy()
    expectedOutputDf["lightMzBin"] = np.append(expectedLightMzBins, expectedLightMzBins)
    expectedOutputDf["heavyMz"] = expectedHeavyMzColumn + expectedHeavyMzColumn
    expectedOutputDf["heavyMzBin"] = np.append(expectedHeavyMzBinColumn, expectedHeavyMzBinColumn)
    outputDf = calculate_binning_information_by_compensation_voltage(inputDf, isIncludeHeavy=isIncludeHeavy)
    assert expectedOutputDf.equals(outputDf)


@pytest.fixture
def inputFormattedDf():
    inputData = [
        ["peptide1",20.0, 30.0],
        ["peptide2",20.0, 30.0],
        ["peptide3",10.0, 20.0],
    ]
    return pd.DataFrame(inputData, columns=["peptide","lightMzBin", "heavyMzBin"])

@pytest.fixture
def targetedReanalysisNoHeavyDf():
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    expectedOutputData = [
        ["1/peptide3", formula, adduct, 10.0, charge, 0],
        ["2/peptide1/peptide2", formula, adduct, 20.0, charge, 1],
    ]
    return pd.DataFrame(expectedOutputData, columns=["Compound","Formula","Adduct","m.z","z","MSXID"])

def test__output_formatting_functions__create_targeted_reanalysis_dataframe__no_heavy(inputFormattedDf, targetedReanalysisNoHeavyDf):
    isIncludeHeavy = False
    outputDf = create_targeted_reanalysis_dataframe(inputFormattedDf, isIncludeHeavy=isIncludeHeavy)
    assert targetedReanalysisNoHeavyDf.equals(outputDf)

@pytest.fixture
def targetedReanalysisWithHeavyDf():
    formula = ""
    adduct = "(no adduct)"
    charge = 2
    expectedOutputData = [
        ["1/peptide3", formula, adduct, 10.0, charge, 0],
        ["1/peptide3", formula, adduct, 20.0, charge, 0],
        ["2/peptide1/peptide2", formula, adduct, 20.0, charge, 1],
        ["2/peptide1/peptide2", formula, adduct, 30.0, charge, 1],
    ]
    return pd.DataFrame(expectedOutputData, columns=["Compound","Formula","Adduct","m.z","z","MSXID"])

def test__output_formatting_functions__create_targeted_reanalysis_dataframe__with_heavy(inputFormattedDf, targetedReanalysisWithHeavyDf):
    isIncludeHeavy = True
    outputDf = create_targeted_reanalysis_dataframe(inputFormattedDf, isIncludeHeavy=isIncludeHeavy)
    assert targetedReanalysisWithHeavyDf.equals(outputDf)

def test__output_formatting_functions__create_targeted_reanalysis_dataframe_by_compensation_voltage__no_heavy(inputFormattedDf, targetedReanalysisNoHeavyDf):
    isIncludeHeavy = False
    inputFormattedDfCV30 = inputFormattedDf.copy()
    inputFormattedDfCV30["CompensationVoltage"] = [-30] * len(inputFormattedDfCV30.index)
    inputFormattedDfCV40 = inputFormattedDf.copy()
    inputFormattedDfCV40["CompensationVoltage"] = [-40] * len(inputFormattedDfCV40.index)
    inputDf = pd.concat([inputFormattedDfCV30, inputFormattedDfCV40])
    expectedOutput = {
        "-30": targetedReanalysisNoHeavyDf,
        "-40": targetedReanalysisNoHeavyDf,
    }
    output = create_targeted_reanalysis_dataframe_by_compensation_voltage__no_heavy(inputDf, isIncludeHeavy)
    for cv, expectedTargetedReanalysisDf in expectedOutput.items():
        assert cv in output
        assert expectedTargetedReanalysisDf.equals(output[cv])

def test__output_formatting_functions__create_targeted_reanalysis_dataframe_by_compensation_voltage__with_heavy(inputFormattedDf, targetedReanalysisWithHeavyDf):
    isIncludeHeavy = True
    inputFormattedDfCV30 = inputFormattedDf.copy()
    inputFormattedDfCV30["CompensationVoltage"] = [-30] * len(inputFormattedDfCV30.index)
    inputFormattedDfCV40 = inputFormattedDf.copy()
    inputFormattedDfCV40["CompensationVoltage"] = [-40] * len(inputFormattedDfCV40.index)
    inputDf = pd.concat([inputFormattedDfCV30, inputFormattedDfCV40])
    expectedOutput = {
        "-30": targetedReanalysisWithHeavyDf,
        "-40": targetedReanalysisWithHeavyDf,
    }
    output = create_targeted_reanalysis_dataframe_by_compensation_voltage__no_heavy(inputDf, isIncludeHeavy)
    for cv, expectedTargetedReanalysisDf in expectedOutput.items():
        assert cv in output
        assert expectedTargetedReanalysisDf.equals(output[cv])
