import random
import re
from pyteomics import mass
from zodiaq.loaders.library.decoyGenerationFunctions import shuffle_peptide_sequence_with_preserved_cleavage_points, calculate_similarities_between_strings, calculate_ion_mz, add_decoys_to_zodiaq_library, determine_if_decoys_should_be_generated
from test_LibraryLoaderStrategyTable import assert_final_dict_output_matches_expected

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points():
    random.seed(0)
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
    random.seed(0)
    testPeptide = 'AAAAAAAAKAAABR'
    kIdx, rIdx = -6, -1
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[kIdx] == testPeptide[kIdx]
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert shuffledPeptide == testPeptide

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points__partially_homogenous_sequence():
    random.seed(0)
    testPeptide = 'AAAAAAAAAAAAAAKABCDR'
    kIdx, rIdx = -6, -1
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[kIdx] == testPeptide[kIdx]
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert shuffledPeptide != testPeptide

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points__modifications_preserved_with_preceding_amino_acid():
    random.seed(0)
    testPeptide = 'LLAK(UniMod:1)C(UniMod:4)HABDEFGHIJLR'
    rIdx = -1
    modifiedKIdx = 4
    modifiedCIdx = 15
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert (shuffledPeptide[modifiedKIdx] != testPeptide[modifiedKIdx]) or (shuffledPeptide[modifiedCIdx] != testPeptide[modifiedCIdx])
    assert shuffledPeptide != testPeptide
    matches = set(re.findall(r'([A-Z](?:[\[\(][^\)\]]+[\]\)]))', shuffledPeptide))
    expectedMatches = set(['K(UniMod:1)','C(UniMod:4)'])
    assert matches == expectedMatches

def test__decoy_generation_functions__shuffle_peptide_sequence_with_preserved_cleavage_points__modifications_at_start_preserved_with_following_amino_acid():
    random.seed(0)
    testPeptide = '(UniMod:1)MDFQ(UniMod:2)HRPGGK'
    rIdx = -1
    modifiedKIdx = 4
    modifiedCIdx = 15
    shuffledPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(testPeptide)
    assert shuffledPeptide[rIdx] == testPeptide[rIdx]
    assert (shuffledPeptide[modifiedKIdx] != testPeptide[modifiedKIdx]) or (shuffledPeptide[modifiedCIdx] != testPeptide[modifiedCIdx])
    assert shuffledPeptide != testPeptide
    assert '1' not in shuffledPeptide
    assert '(UniMod:2)' in shuffledPeptide
    assert len('MDFQ(UniMod:2)HRPGGK') == len(shuffledPeptide)

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
        assert abs(mzs[i] - calculatedMz) < 0.0001

def test__decoy_generation_functions__calculate_ion_mz__with_unimod_modification():
    originalPeptide = 'ABC(UniMod:4)R'
    charge = 1
    fragment1 = ('y', 1, charge)
    fragment2 = ('y', 2, charge)
    unimod4Mass = 57.021464
    fragment1ExpectedMass = mass.calculate_mass(sequence='R', ion_type='y', charge=charge)
    fragment1Mass = calculate_ion_mz(originalPeptide, *fragment1)
    assert fragment1Mass == fragment1ExpectedMass
    fragment2ExpectedMass = mass.calculate_mass(sequence='CR', ion_type='y', charge=charge) + unimod4Mass
    fragment2Mass = calculate_ion_mz(originalPeptide, *fragment2)
    assert fragment2Mass == fragment2ExpectedMass


def test__decoy_generation_functions__add_decoys_to_zodiaq_library():
    random.seed(0)
    zodiaqLibraryDict = {
        (375.873226, "FANYIDKVR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "1/P08670",
            "peaks": [
                (175.118953, 2926.18, 0),
                (274.187367, 1647.689, 0),
                (333.155733, 1071.177, 0),
                (397.231972, 2078.822, 0),
                (402.282331, 4932.288, 0),
                (454.253437, 1301.4617, 0),
                (489.771994, 1395.553, 0),
                (517.309275, 10000.0, 0),
                (630.393339, 8233.006, 0),
                (793.456668, 5096.472, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
            "fragmentTypes": [
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
            ],
        }
    }
    outputDict = add_decoys_to_zodiaq_library(zodiaqLibraryDict.copy())
    expectedOutputDict = {
        (375.873226, "INALDYKVR"): {
            "precursorCharge": 3,
            "identification": "INALDYKVR_375.873226_DECOY",
            "proteinName": "1/DECOY_P08670",
            "peaks": [
                (175.11895217407, 2926.18, 1),
                (274.18736608706, 1647.689, 1),
                (299.17138166975, 1071.177, 1),
                (397.23197055067004, 2078.822, 1),
                (402.28232910106, 4932.288, 1),
                (432.750527443025, 1301.4617, 1),
                (489.77199116359503, 1395.553, 1),
                (565.34565763361, 10000.0, 1),
                (680.37260065744, 8233.006, 1),
                (793.4566646345701, 5096.472, 1),
            ],
            "zodiaqKeyIdx": 1,
            "isDecoy": 1,
            "fragmentTypes": [
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
            ],
        }
    } | zodiaqLibraryDict.copy()
    assert_final_dict_output_matches_expected(outputDict, expectedOutputDict)

def test__decoy_generation_functions__add_decoys_to_zodiaq_library__ambiguous_peptides_ignored():
    zodiaqLibraryDict = {
        (375.873226, "AAAAAAAAKAAABR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "1/P08670",
            "peaks": [
                (175.11895217407, 2926.18, 1),
                (274.18736608706, 1647.689, 1),
                (299.17138166975, 1071.177, 1),
                (397.23197055067004, 2078.822, 1),
                (402.28232910106, 4932.288, 1),
                (432.750527443025, 1301.4617, 1),
                (489.77199116359503, 1395.553, 1),
                (565.34565763361, 10000.0, 1),
                (680.37260065744, 8233.006, 1),
                (793.4566646345701, 5096.472, 1),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
            "fragmentTypes": [
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
            ],
        }
    }
    random.seed(0)
    outputDict = add_decoys_to_zodiaq_library(zodiaqLibraryDict)
    assert_final_dict_output_matches_expected(outputDict, zodiaqLibraryDict)

def test__decoy_generation_functions__determine_if_decoys_should_be_generated():
    zodiaqLibraryDict_shouldGenerate = {
        (375.873226, "FANYIDKVR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "1/P08670",
            "peaks": [
                (175.118953, 2926.18, 0),
                (274.187367, 1647.689, 0),
                (333.155733, 1071.177, 0),
                (397.231972, 2078.822, 0),
                (402.282331, 4932.288, 0),
                (454.253437, 1301.4617, 0),
                (489.771994, 1395.553, 0),
                (517.309275, 10000.0, 0),
                (630.393339, 8233.006, 0),
                (793.456668, 5096.472, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
            "fragmentTypes": [
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
            ],
        }
    }
    zodiaqLibraryDict_shouldNotGenerate_hasDecoys = {
        (375.873226, "INALDYKVR"): {
            "precursorCharge": 3,
            "identification": "INALDYKVR_375.873226_DECOY",
            "proteinName": "1/DECOY_P08670",
            "peaks": [
                (175.201, 2926.18, 1),
                (274.332, 1647.689, 1),
                (299.33872, 1071.177, 1),
                (397.46125, 2078.822, 1),
                (402.5043, 4932.288, 1),
                (433.00019999999995, 1301.4617, 1),
                (490.05150000000003, 1395.553, 1),
                (565.6775, 10000.0, 1),
                (680.7649, 8233.006, 1),
                (793.9225, 5096.472, 1),
            ],
            "zodiaqKeyIdx": 1,
            "isDecoy": 1,
            "fragmentTypes": [
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
            ],
        }
    } | zodiaqLibraryDict_shouldGenerate
    zodiaqLibraryDict_shouldNotGenerate_hasNoDecoysButNoFragmentTypes = {
        (375.873226, "FANYIDKVR"): {
            "precursorCharge": 3,
            "identification": "FANYIDKVR",
            "proteinName": "1/P08670",
            "peaks": [
                (175.118953, 2926.18, 0),
                (274.187367, 1647.689, 0),
                (333.155733, 1071.177, 0),
                (397.231972, 2078.822, 0),
                (402.282331, 4932.288, 0),
                (454.253437, 1301.4617, 0),
                (489.771994, 1395.553, 0),
                (517.309275, 10000.0, 0),
                (630.393339, 8233.006, 0),
                (793.456668, 5096.472, 0),
            ],
            "zodiaqKeyIdx": 0,
            "isDecoy": 0,
        }
    }
    assert determine_if_decoys_should_be_generated(zodiaqLibraryDict_shouldGenerate)
    assert not determine_if_decoys_should_be_generated(zodiaqLibraryDict_shouldNotGenerate_hasDecoys)
    assert not determine_if_decoys_should_be_generated(zodiaqLibraryDict_shouldNotGenerate_hasNoDecoysButNoFragmentTypes)


