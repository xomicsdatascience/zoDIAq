import re
import random
from pyteomics import mass
from zodiaq.utils import format_protein_string_to_list, format_protein_list_to_string
from zodiaq.loaders.library.modificationMassDict import modificationMassDict
nonCleavageAminoAcids = [ "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W", "Y", "V"]
cleavageAminoAcids = set(['K', 'P', 'R'])
modificationRegexStr = r'([A-Z](?:[\[\(][^\)\]]+[\]\)])?)'

def shuffle_peptide_sequence_with_preserved_cleavage_points(peptide):
    peptide = re.sub('^(?:[\[\(][^\)\]]+[\]\)])', '', peptide)
    originalPeptide = re.findall(modificationRegexStr, peptide)
    cleavageAALocations = [(i-len(originalPeptide)+1,char) for (i, char) in enumerate(originalPeptide) if char in cleavageAminoAcids]
    otherAAs = [char for char in originalPeptide if char not in cleavageAminoAcids]
    for i in range(100):
        shuffledPeptide = shuffle_non_cleavage_amino_acids(otherAAs, i)
        shuffledPeptide = insert_cleavage_amino_acids_into_shuffled_peptide(shuffledPeptide, cleavageAALocations, originalPeptide)
        if calculate_similarities_between_strings(originalPeptide, shuffledPeptide) < 0.7:
            return ''.join(shuffledPeptide)
    if calculate_similarities_between_strings(originalPeptide, shuffledPeptide) >= 0.7:
        return peptide
    return ''.join(shuffledPeptide)

def shuffle_non_cleavage_amino_acids(otherAAs, i):
    shuffledPeptide = otherAAs[:]
    random.shuffle(shuffledPeptide)
    if i % 10 == 0:
        randomly_swap_single_amino_acid_in_shuffled_peptide_to_increase_variability(shuffledPeptide)
    return shuffledPeptide

def randomly_swap_single_amino_acid_in_shuffled_peptide_to_increase_variability(shuffledPeptide):
    randomIdx = random.randint(0, len(shuffledPeptide) - 1)
    shuffledPeptide[randomIdx] = random.choice(nonCleavageAminoAcids)
    return shuffledPeptide

def insert_cleavage_amino_acids_into_shuffled_peptide(shuffledPeptide, cleavageAALocations, peptide):
    for j in range(len(cleavageAALocations) - 1, -1, -1):
        cleavageAAIndex = cleavageAALocations[j][0]
        if not cleavageAAIndex:
            cleavageAAIndex = len(peptide)
        cleavageAA = cleavageAALocations[j][1]
        shuffledPeptide.insert(cleavageAAIndex, cleavageAA)
    return shuffledPeptide

def insert_value(string, index, char):
    return string[:index] + [char] + string[index+1:]

def calculate_similarities_between_strings(s1, s2):
    minLength = min(len(s1),len(s2))
    similarities = 0
    for i in range(minLength):
        if s1[i] == s2[i]: similarities += 1
    return similarities / minLength


def calculate_ion_mz(sequence, type, position, charge):
    terminalSequence = determine_terminal_end_of_sequence(sequence, type, position)
    mass = calculate_molecular_mass_of_sequence(terminalSequence, type, charge)
    return mass / charge


def determine_terminal_end_of_sequence(sequence, type, position):
    modifiedSequence = re.findall(modificationRegexStr, sequence)
    if type == 'y':
        return ''.join(modifiedSequence[-position:])
    else:
        return ''.join(modifiedSequence[:position])


def calculate_molecular_mass_of_sequence(sequence, type, charge):
    modifications = re.findall(r'([\[\(][^\)\]]+[\]\)])', sequence)
    unmodifiedSequence = re.sub(r'([\[\(][^\)\]]+[\]\)])', '', sequence)
    unmodifiedMass = mass.calculate_mass(sequence=unmodifiedSequence, ion_type=type, charge=charge)
    modificationMass = sum([modificationMassDict[m] for m in modifications])
    return unmodifiedMass + modificationMass

def add_decoys_to_zodiaq_library(zodiaqLibraryDict):
    decoyDict = {}
    for key, value in zodiaqLibraryDict.items():
        decoyKey, decoyValue = create_decoy_from_target_in_zodiaq_library(key, value.copy())
        if decoyKey:
            decoyDict[decoyKey] = decoyValue
    outputDict = recalculate_zodiaq_library_key_idx_to_account_for_inserted_decoys(decoyDict | zodiaqLibraryDict)
    return outputDict

def create_decoy_from_target_in_zodiaq_library(targetKey, targetValue):
    targetPeptide = targetKey[1]
    targetPeaks = targetValue['peaks']
    fragmentTypes = targetValue['fragmentTypes']
    decoyPeptide = shuffle_peptide_sequence_with_preserved_cleavage_points(targetPeptide)
    if decoyPeptide == targetPeptide:
        return None, None
    decoyMzs = [calculate_ion_mz(decoyPeptide, *fragmentType) for fragmentType in fragmentTypes]
    decoyPeaks = [(decoyMzs[i], targetPeaks[i][1], -1) for i in range(len(targetPeaks))]
    targetProteinString = targetValue['proteinName']
    decoyProteins = format_protein_string_to_list(targetProteinString)
    decoyProteins = [f'DECOY_{protein}' for protein in decoyProteins]
    decoyProteinString = format_protein_list_to_string(decoyProteins)
    decoyValue = targetValue
    decoyValue['identification'] = f'{decoyPeptide}_{targetKey[0]}_DECOY'
    decoyValue['proteinName'] = decoyProteinString
    decoyValue['peaks'] = decoyPeaks
    decoyValue['isDecoy'] = 1
    decoyValue['zodiaqKeyIdx'] = -1
    return (targetKey[0], decoyPeptide), decoyValue

def recalculate_zodiaq_library_key_idx_to_account_for_inserted_decoys(outputDict):
    sortedKeys = sorted(outputDict)
    for i in range(len(sortedKeys)):
        key = sortedKeys[i]
        value = outputDict[key]
        value['peaks'] = [(x[0],x[1],i) for x in value['peaks']]
        value['zodiaqKeyIdx'] = i
        outputDict[key] = value
    return outputDict

def determine_if_decoys_should_be_generated(zodiaqLibraryDict):
    if there_are_any_decoys_in_library(zodiaqLibraryDict): return False
    if any_library_entries_dont_have_fragment_type_data_requirement_for_decoy_generation(zodiaqLibraryDict): return False
    return True

def there_are_any_decoys_in_library(zodiaqLibraryDict):
    for key, value in zodiaqLibraryDict.items():
        if value['isDecoy']: return True
    return False

def any_library_entries_dont_have_fragment_type_data_requirement_for_decoy_generation(zodiaqLibraryDict):
    for key, value in zodiaqLibraryDict.items():
        if 'fragmentTypes' not in value: return True
    return False



