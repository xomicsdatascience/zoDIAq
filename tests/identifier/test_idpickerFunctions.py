import pandas as pd
from csodiaq.identifier.idpickerFunctions import initialize__format_peptide_protein_connections, collapse__group_identically_connected_peptides_and_proteins, separate__identify_and_label_independent_clusters, reduce__identify_minimum_number_of_most_connected_proteins, identify_high_confidence_proteins
import numpy as np

def test__idpicker_functions__initialize__format_peptide_protein_connections():
    peptideProteinData = [
        ['peptide01', '1/protein7'],
        ['peptide02', '3/protein4/protein6/protein9'],
        ['peptide03', '1/protein1'],
        ['peptide04', '2/protein1/protein5'],
        ['peptide05', '1/protein7'],
        ['peptide06', '2/protein3/protein6'],
        ['peptide07', '1/protein1'],
        ['peptide08', '4/protein1/protein2/protein5/protein8'],
        ['peptide09', '1/protein1'],
        ['peptide10', '2/protein4/protein9'],
    ]
    peptideProteinDf = pd.DataFrame(peptideProteinData, columns=["peptide","protein"])
    expectedOutputData = [
        ('peptide01', 'protein7'),
        ('peptide02', 'protein4'),
        ('peptide02', 'protein6'),
        ('peptide02', 'protein9'),
        ('peptide03', 'protein1'),
        ('peptide04', 'protein1'),
        ('peptide04', 'protein5'),
        ('peptide05', 'protein7'),
        ('peptide06', 'protein3'),
        ('peptide06', 'protein6'),
        ('peptide07', 'protein1'),
        ('peptide08', 'protein1'),
        ('peptide08', 'protein2'),
        ('peptide08', 'protein5'),
        ('peptide08', 'protein8'),
        ('peptide09', 'protein1'),
        ('peptide10', 'protein4'),
        ('peptide10', 'protein9'),
    ]
    expectedOutputDf = pd.DataFrame(expectedOutputData, columns=["peptide","protein"])
    output = initialize__format_peptide_protein_connections(peptideProteinDf)
    assert expectedOutputDf.equals(output)


def test__idpicker_functions__collapse__group_identically_connected_peptides_and_proteins():
    data = [
        (1,7),
        (2,4),
        (2,6),
        (2,9),
        (3,1),
        (4,1),
        (4,5),
        (5,7),
        (6,3),
        (6,6),
        (7,1),
        (8,1),
        (8,2),
        (8,5),
        (8,8),
        (9,1),
        (10,4),
        (10,9)
    ]
    df = pd.DataFrame(data, columns=["peptide","protein"])
    groupedData = [
        ((1, 5), (7,)),
        ((2,), (4, 9)),
        ((2,), (6,)),
        ((3, 7, 9), (1,)),
        ((4,), (1,)),
        ((4,), (5,)),
        ((6,), (3,)),
        ((6,), (6,)),
        ((8,), (1,)),
        ((8,), (2, 8)),
        ((8,), (5,)),
        ((10,), (4, 9))
    ]
    expectedOutput = pd.DataFrame(groupedData, columns=["peptide","protein"])
    output = collapse__group_identically_connected_peptides_and_proteins(df)
    assert expectedOutput.equals(output)

def test__idpicker_functions__separate__identify_and_label_independent_clusters():
    data = [
        ((1, 5), (7,)),
        ((2,), (4, 9)),
        ((2,), (6,)),
        ((3, 7, 9), (1,)),
        ((4,), (1,)),
        ((4,), (5,)),
        ((6,), (3,)),
        ((6,), (6,)),
        ((8,), (1,)),
        ((8,), (2, 8)),
        ((8,), (5,)),
        ((10,), (4, 9))
    ]
    df = pd.DataFrame(data, columns=["peptide","protein"])
    expectedClusters = np.array([
        0,
        1,
        1,
        2,
        2,
        2,
        1,
        1,
        2,
        2,
        2,
        1
    ])
    clusters = separate__identify_and_label_independent_clusters(df)
    np.testing.assert_array_equal(expectedClusters, clusters)

def test__idpicker_functions__reduce__identify_minimum_number_of_most_connected_proteins():
    data = [
        ((1, 5), (7,)),
        ((2,), (4, 9)),
        ((2,), (6,)),
        ((3, 7, 9), (1,)),
        ((4,), (1,)),
        ((4,), (5,)),
        ((6,), (3,)),
        ((6,), (6,)),
        ((8,), (1,)),
        ((8,), (2, 8)),
        ((8,), (5,)),
        ((10,), (4, 9))
    ]
    df = pd.DataFrame(data, columns=["peptide","protein"])
    clusters = np.array([
        0,
        1,
        1,
        2,
        2,
        2,
        1,
        1,
        2,
        2,
        2,
        1
    ])
    df["cluster"] = clusters
    expectedProteins = set([
        (1,),
        (4, 9),
        (6,),
        (7,),
    ])
    proteins = reduce__identify_minimum_number_of_most_connected_proteins(df)
    assert expectedProteins == proteins

def test__idpicker_functions__identify_high_confidence_proteins():
    peptideProteinData = [
        ['peptide01', '1/protein7'],
        ['peptide02', '3/protein4/protein6/protein9'],
        ['peptide03', '1/protein1'],
        ['peptide04', '2/protein1/protein5'],
        ['peptide05', '1/protein7'],
        ['peptide06', '2/protein3/protein6'],
        ['peptide07', '1/protein1'],
        ['peptide08', '4/protein1/protein2/protein5/protein8'],
        ['peptide09', '1/protein1'],
        ['peptide10', '2/protein4/protein9'],
    ]
    peptideProteinDf = pd.DataFrame(peptideProteinData, columns=["peptide","protein"])
    expectedProteins = set([
        ('protein1',),
        ('protein4', 'protein9'),
        ('protein6',),
        ('protein7',),
    ])
    proteins = identify_high_confidence_proteins(peptideProteinDf)
    assert expectedProteins == proteins

