from collections import defaultdict

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

peptides, proteins = map(list,zip(*data))

def group_nodes_with_same_edge(l1, l2):
    initialDict = {}
    for i in range(len(l1)):
        if l1[i] not in initialDict: initialDict[l1[i]] = [l2[i]]
        else: initialDict[l1[i]].append(l2[i])

    reverseDict = defaultdict(list)
    for key, value in initialDict.items():
        reverseDict[tuple(sorted(value))].append(key)

    groups = {}
    for value in dict.values():
        for x in value:
            groups[x] = value

    return groups


proteinGroups = {key:[key] for key in proteins}
initialDict = {}
for i in range(len(proteins)):
    if proteins[i] not in initialDict: initialDict[proteins[i]] = [peptides[i]]
    else: initialDict[proteins[i]].append(peptides[i])

dict = defaultdict(list)
for key, value in initialDict.items():
    dict[tuple(sorted(value))].append(key)
for value in dict.values():
    if len(value) > 1:
        for x in value:
            proteinGroups[x] = value


peptideGroups = {key:[key] for key in peptides}
for key in peptideGroups:
    print('key: '+str(key))
    print('value: '+str(dict[key]))

initialDict = {}
for i in range(len(peptides)):
    if peptides[i] not in initialDict: initialDict[peptides[i]] = [proteins[i]]
    else: initialDict[peptides[i]].append(proteins[i])

dict = defaultdict(list)
for key, value in initialDict.items():
    dict[tuple(sorted(value))].append(key)
for key, value in dict.items():
    if len(value) > 1:
        for x in value:
            peptideGroups[x] = value


print('peptides-------------')
for key in peptideGroups:
    print('key: '+str(key))
    print('value: '+str(peptideGroups[key]))

print('proteins-------------')
for key in proteinGroups:
    print('key: '+str(key))
    print('value: '+str(proteinGroups[key]))


# Step 1 - consolidating
