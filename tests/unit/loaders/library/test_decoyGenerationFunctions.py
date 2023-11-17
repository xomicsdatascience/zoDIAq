from zodiaq.loaders.library.decoyGenerationFunctions import shuffle_peptide_sequence_with_preserved_cleavage_points, calculate_similarities_between_strings, calculate_ion_mz


def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points():
    testPeptide = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    kIdx, pIdx, rIdx = 10, 15, 17
    aIdx, lIdx, qIdx, zIdx = 0, 11, 16, 25

    mismatches = []
    numRepeats = 20
    for i in range(numRepeats):
        shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
        mismatches.append(testPeptide[aIdx]!=shuffledPeptide[aIdx])
        mismatches.append(testPeptide[lIdx]!=shuffledPeptide[lIdx])
        mismatches.append(testPeptide[qIdx]!=shuffledPeptide[qIdx])
        mismatches.append(testPeptide[zIdx]!=shuffledPeptide[zIdx])
        assert testPeptide[kIdx] == shuffledPeptide[kIdx]
        assert testPeptide[pIdx] == shuffledPeptide[pIdx]
        assert testPeptide[rIdx] == shuffledPeptide[rIdx]
    mismatches = [m for m in mismatches if m]
    assert len(mismatches) > numRepeats * 3.5

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points__completely_homogenous_sequence():
    testPeptide = 'AAAAAAAAKAAABR'
    kIdx, rIdx = -6, -1
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[kIdx] == testPeptide[kIdx]
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert shuffledPeptide == testPeptide

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points__partially_homogenous_sequence():
    testPeptide = 'AAAAAAAAAAAAAAKABCDR'
    kIdx, rIdx = -6, -1
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[kIdx] == testPeptide[kIdx]
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert shuffledPeptide != testPeptide


def test__decoy_generation_functions__calculate_differences_between_strings():
    numChars = 10
    originalString = 'A' * 10
    changedString = originalString[:]
    for i in range(numChars):
        changedString = changedString[:i] + 'B' + changedString[i+1:]
        similarityProportion = calculate_similarities_between_strings(originalString, changedString)
        assert similarityProportion == (numChars-(i+1))/numChars

def test__decoy_generation_functions__calculate_ion_mz():
    sequence = "FANYIDKVR"
    mzs = [
        175.118953,
        274.187367,
        333.155733,
        397.231972,
        402.282331,
        454.253437,
        489.771994,
        517.309275,
        630.393339,
        793.456668,
    ]
    fragments = [
        ('y', 1, 1),
        ('y', 2, 1),
        ('b', 3, 1),
        ('y', 6, 2),
        ('y', 3, 1),
        ('y', 7, 2),
        ('y', 8, 2),
        ('y', 4, 1),
        ('y', 5, 1),
        ('y', 6, 1),
    ]
    for i in range(len(mzs)):
        calculatedMz = calculate_ion_mz(sequence, *fragments[i])
        assert abs(mzs[i] - calculatedMz) < 0.5