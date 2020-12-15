from collections import defaultdict, OrderedDict

# step 1 - initialize
testData = [
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

# step 2 - collapse
def group_nodes_with_same_edge(listOfTupleConnections, first=True):
    if first: l1, l2 = map(list,zip(*listOfTupleConnections))
    else: l2, l1 = map(list,zip(*listOfTupleConnections))
    initialDict = {}
    for i in range(len(l1)):
        if l1[i] not in initialDict: initialDict[l1[i]] = [l2[i]]
        else: initialDict[l1[i]].append(l2[i])

    reverseDict = defaultdict(list)
    for key, value in initialDict.items():
        reverseDict[tuple(sorted(value))].append(key)

    groups = {}
    for value in reverseDict.values():
        for x in value:
            groups[x] = value

    for i in range(len(l1)):
        l1[i] = tuple(groups[l1[i]])


    if first: return list(set(tuple(zip(l1, l2))))
    else: return list(set(tuple(zip(l2, l1))))


# Step 3 - separate
def create_node_dict(data, first=True):
    v1, v2 = 0, 1
    if not first: v1, v2 = 1, 0

    nodeDict = defaultdict(set)
    compList = []
    for x in data:
        nodeDict[x[v1]].update(set(x[v2]))
    return nodeDict


def separate_into_clusters(data):
    nodeDict = create_node_dict(data)
    usedKeys = set()
    clusters = []
    while len(usedKeys) != len(nodeDict):
        unusedKeys = set(nodeDict.keys()) - usedKeys
        root = list(unusedKeys)[0]
        proteins = nodeDict[root]
        cluster = set([root])
        unusedKeys -= cluster
        priorLength = -1
        while priorLength != len(proteins):
            priorLength = len(proteins)
            for key in unusedKeys:
                if len(nodeDict[key].intersection(proteins)) > 0:
                    cluster.add(key)
                    proteins.update(nodeDict[key])
            unusedKeys -= cluster
        clusters.append(cluster)
        usedKeys.update(cluster)
    finalClusters = [[] for i in clusters]
    for x in data:
        for i in range(len(clusters)):
            if x[0] in clusters[i]: finalClusters[i].append(x)
    return finalClusters


# step 4 - reduce
def reduce_cluster(data):
    protDict, lenDict = defaultdict(list), defaultdict(int)
    totalPepSet = set()
    for x in data: protDict[x[1]].append(x[0]); lenDict[x[1]] += 1; totalPepSet.add(x[0])
#    lenDict = OrderedDict(sorted(lenDict.items()))

    sortedKeys = sorted(lenDict, key=lenDict.get, reverse=True)
    pepSet = set()
    proteins = []
    limit = -1

    for key in sortedKeys:
#        if len(pepSet) == len(totalPepSet): break
        if len(pepSet) == len(totalPepSet): limit = lenDict[key]
        if lenDict[key] < limit or limit == 1: break
        oldPepSetLength = len(pepSet)
        pepSet.update(set(protDict[key]))
        if len(pepSet) != oldPepSetLength or lenDict[key] > 1: proteins += list(key)
    return proteins

# comprehensive
def find_valid_proteins(data):

    data = group_nodes_with_same_edge(data)
    data = group_nodes_with_same_edge(data, first=False)

    cl = separate_into_clusters(data)
    finalProteins = []
    for c in cl: finalProteins += reduce_cluster(c)
    return set(finalProteins)

#finalProteins = find_valid_proteins(testData)
#print('final proteins: ' + str(finalProteins))
