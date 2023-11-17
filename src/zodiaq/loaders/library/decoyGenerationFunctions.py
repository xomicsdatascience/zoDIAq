import random
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

nonCleavageAminoAcids = [ "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W", "Y", "V"]
cleavageAminoAcids = set(['K', 'P', 'R'])

def shuffle_peptide_sequence_with_preserved_cleavage_points(peptide):
    cleavageAALocations = [(i-len(peptide)+1,char) for (i, char) in enumerate(peptide) if char in cleavageAminoAcids]
    otherAAs = [char for char in peptide if char not in cleavageAminoAcids]
    for i in range(100):
        shuffledPeptide = shuffle_non_cleavage_amino_acids(otherAAs, i)
        shuffledPeptideString = insert_cleavage_amino_acids_into_shuffled_peptide(shuffledPeptide, cleavageAALocations, peptide)
        if calculate_similarities_between_strings(peptide, shuffledPeptideString) < 0.7:
            return shuffledPeptideString
    if calculate_similarities_between_strings(peptide, shuffledPeptideString) >= 0.7:
        return peptide
    return shuffledPeptideString

def shuffle_non_cleavage_amino_acids(otherAAs, i):
    shuffledPeptide = otherAAs[:]
    random.shuffle(shuffledPeptide)
    if i % 10 == 0:
        insert_random_amino_acid_to_shuffled_peptide_to_increase_variability(shuffledPeptide)
    return shuffledPeptide

def insert_random_amino_acid_to_shuffled_peptide_to_increase_variability(shuffledPeptide):
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
    return ''.join(shuffledPeptide)

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
    if type == 'y':
        return Seq(sequence[-position:])
    else:
        return Seq(sequence[:position])


def calculate_molecular_mass_of_sequence(sequence, type, charge):
    H2O = 18.01468
    mass = molecular_weight(sequence, seq_type='protein')
    if type == 'b':
        mass -= H2O
    return mass + charge