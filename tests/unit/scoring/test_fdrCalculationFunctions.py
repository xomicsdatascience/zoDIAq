import pytest
import pandas as pd

from csodiaq.scoring.fdrCalculationFunctions import (
    drop_duplicate_values_from_df_in_given_column,
    create_spectral_fdr_output_from_full_output_sorted_by_desired_score,
    create_peptide_fdr_output_from_full_output_sorted_by_desired_score,
    identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff,
    organize_peptide_df_by_leading_proteins,
    determine_if_peptides_are_unique_to_leading_protein,
)


def test__fdr_calculation_functions__drop_duplicate_values_from_df_in_given_column():
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


def test__fdr_calculation_functions__create_spectral_fdr_output_from_full_output():
    numNonDecoys = 100
    numDecoys = 2
    isDecoyColumn = ([0] * numNonDecoys) + ([1] * numDecoys)
    scoreColumn = list(range(len(isDecoyColumn) - 1, -1, -1))
    inputDf = pd.DataFrame()
    inputDf["isDecoy"] = isDecoyColumn

    expectedOutputDf = pd.DataFrame()
    expectedOutputDf["isDecoy"] = isDecoyColumn[:-1]
    expectedOutputDf["spectralFDR"] = [0] * numNonDecoys + [1 / (numNonDecoys + 1)]

    outputDf = create_spectral_fdr_output_from_full_output_sorted_by_desired_score(
        inputDf
    )
    assert expectedOutputDf.equals(outputDf)


def test__fdr_calculation_functions__create_peptide_fdr_output_from_full_output():
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
    outputDf = create_peptide_fdr_output_from_full_output_sorted_by_desired_score(
        inputDf
    )

    assert expectedOutputDf.equals(outputDf)


def test__fdr_calculation_functions__organize_peptide_df_by_leading_proteins():
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


def test__fdr_calculation_functions__identify_leading_protein_to_fdr_dictionary_for_leading_proteins_below_fdr_cutoff():
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


def test__fdr_calculation_functions__determine_if_peptides_are_unique_to_leading_protein():
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
