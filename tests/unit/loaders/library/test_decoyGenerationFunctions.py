from zodiaq.loaders.library.decoyGenerationFunctions import shuffle_peptide_sequence_with_preserved_cleavage_points, calculate_similarities_between_strings


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
