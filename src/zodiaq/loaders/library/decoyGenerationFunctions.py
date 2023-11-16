import random
nonKRAminoAcids = [ "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W", "Y", "V"]

def shuffle_peptide_sequence_with_preserved_cleavage_points(peptide):
    cleavageAAs = set(['K','P','R'])
    cleavageAALocations = [(i-len(peptide)+1,char) for (i, char) in enumerate(peptide) if char in cleavageAAs]
    otherAAs = [char for char in peptide if char not in cleavageAAs]
    shuffledPeptide = ''
    for i in range(100):
        shuffledPeptide = otherAAs[:]
        random.shuffle(shuffledPeptide)
        if i % 10 == 0:
            randomIdx = random.randint(0,len(shuffledPeptide)-1)
            shuffledPeptide[randomIdx] = random.choice(nonKRAminoAcids)
        shuffledPeptide[0]
        for j in range(len(cleavageAALocations)-1, -1, -1):
            cleavageAAIndex = cleavageAALocations[j][0]
            if not cleavageAAIndex:
                cleavageAAIndex = len(peptide)
            cleavageAA = cleavageAALocations[j][1]
            shuffledPeptide.insert(cleavageAAIndex, cleavageAA)
        shuffledPeptideString = ''.join(shuffledPeptide)
        if calculate_similarities_between_strings(peptide, shuffledPeptideString) < 0.7:
            return shuffledPeptideString
    if calculate_similarities_between_strings(peptide, shuffledPeptideString) >= 0.7:
        return peptide
    return shuffledPeptideString

def insert_value(string, index, char):
    return string[:index] + [char] + string[index+1:]

def calculate_similarities_between_strings(s1, s2):
    minLength = min(len(s1),len(s2))
    similarities = 0
    for i in range(minLength):
        if s1[i] == s2[i]: similarities += 1
    return similarities / minLength