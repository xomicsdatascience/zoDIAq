import pandas as pd
from csodiaq.identifier.idpickerFunctions import collapse__group_identically_connected_peptides_and_proteins, separate__identify_and_label_independent_clusters
import numpy as np

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



