import pytest
import pandas as pd

from csodiaq.scoring.quantificationFunctions import (
    calculate_ion_count_from_peptides_of_protein,
    calculate_ion_count_for_each_protein_in_protein_fdr_df,
    compile_ion_count_comparison_across_runs_df,
)


@pytest.fixture
def ionCountList():
    return [
        100.0,
        200.0,
        300.0,
    ]


def test__fdr_calculation_functions__calculate_ion_count_from_peptides_of_protein(
    ionCountList,
):
    expectedOutput = 200.0
    output = calculate_ion_count_from_peptides_of_protein(ionCountList)
    assert expectedOutput == output


def test__fdr_calculation_functions__calculate_ion_count_for_each_protein_in_protein_fdr_df():
    inputData = [
        ["1/protein1", 100.0],
        ["1/protein1", 200.0],
        ["1/protein1", 300.0],
        ["2/protein2/protein3", 400.0],
        ["2/protein2/protein3", 500.0],
    ]

    inputDf = pd.DataFrame(inputData, columns=["leadingProtein", "ionCount"])
    expectedOutputData = [
        ["protein1", 200.0],
        ["protein2", 450.0],
        ["protein3", 450.0],
    ]
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=["protein", "ionCount"])
    outputDf = calculate_ion_count_for_each_protein_in_protein_fdr_df(inputDf)
    assert expectedOutputDf.equals(outputDf)


def test__fdr_calculation_functions__compile_ion_count_comparison_across_runs_df():
    columnName = "value"
    inputDf1 = pd.DataFrame(
        [
            ["p1", 100.0],
            ["p2", 200.0],
            ["p3", 300.0],
        ],
        columns=[columnName, "ionCount"],
    )
    inputDf2 = pd.DataFrame(
        [
            ["p3", 400.0],
            ["p4", 500.0],
            ["p5", 600.0],
        ],
        columns=[columnName, "ionCount"],
    )
    inputDf3 = pd.DataFrame(
        [
            ["p2", 700.0],
            ["p3", 800.0],
            ["p4", 900.0],
            ["p6", 1000.0],
        ],
        columns=[columnName, "ionCount"],
    )
    inputDfs = {
        "df1": inputDf1,
        "df2": inputDf2,
        "df3": inputDf3,
    }
    expectedOutputDf = pd.DataFrame(
        [
            [100.0, 200.0, 300.0, 0, 0, 0],
            [0, 0, 400.0, 500.0, 600.0, 0],
            [0, 700.0, 800.0, 900.0, 0, 1000.0],
        ],
        columns=[f"p{i}" for i in range(1, 7)],
        index=["df1", "df2", "df3"],
    )
    outputDf = compile_ion_count_comparison_across_runs_df(inputDfs, columnName)
    assert expectedOutputDf.equals(outputDf)
