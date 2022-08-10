from collections import defaultdict, OrderedDict

'''
Function: group_nodes_with_same_edge()
Purpose: Provided a list of connections, this function groups left-type or right-type (depending on value of 'first' boolean
            variable) that have identical connections. For example, if 'A' only connects with '1' and '2', and 'B' also only
            connects with '1' and '2', 'A' and 'B' will be combined into one group 'A,B'.
        NOTE: In our case, peptides are the first value, and proteins are the second value.
Parameters:
    'data' - list of 2-value tuples. Each tuple represents a connection between a peptide (first value) and
        a protein (second value). The function can accept all types inside the tuples, but generally they will be strings
        (single peptide or protein) or tuples (peptide or protein group).
    'first' - boolean value that indicates if peptides (value=True, default) or proteins (value=False) are being grouped.
Returns:
    list of 2-value tuples. Represents the same values as the input, but identically connected peptides/proteins have been
        combined. All values of the peptide/protein, regardless of whether or not the peptide/proteins have been grouped with
        others, are now cast as tuples.
'''


def group_nodes_with_same_edge(data, first=True):

    # separates the list of tuple connections into two lists of matching indices.
    # Note: which group of indices that are being "grouped" depends on the value of the function parameter 'first'
    if first:
        l1, l2 = map(list, zip(*data))
    else:
        l2, l1 = map(list, zip(*data))

    # all connections relating to a single node are consolidated into a dictionary.
    initialDict = {}
    for i in range(len(l1)):
        if l1[i] not in initialDict:
            initialDict[l1[i]] = [l2[i]]
        else:
            initialDict[l1[i]].append(l2[i])

    # nodes with identical connections are grouped together.
    reverseDict = defaultdict(list)
    for key, value in initialDict.items():
        reverseDict[tuple(sorted(value))].append(key)

    # a new dictionary is made, with key:value being node:node group
    groups = {}
    for value in reverseDict.values():
        for x in value:
            groups[x] = value

    # all nodes are replaced with their node group in the original list
    for i in range(len(l1)):
        l1[i] = tuple(groups[l1[i]])

    # The two lists are zipped back up into tuples, placed in a set to remove duplicates, and then returned as a new list
    if first:
        return list(set(tuple(zip(l1, l2))))
    else:
        return list(set(tuple(zip(l2, l1))))


# step 3 - cluster
'''
Function: separate_into_clusters()
Purpose: Given a list of node connections, this function separates them into clusters based on interconnectibility. In essence,
            there are no connections between two or more clusters - each is a uniquely interconnected group.
Parameters:
    'data' - list of (tuple, tuple) tuples. Each tuple in data represents a connection between a peptide group (first tuple)
        and a protein group (second tuple). Each sub-tuple contains an undetermined number of strings, depending on the number
        of peptides/proteins in the group.
Returns:
    list - list - (tuple, tuple) tuples. In essence, the original list has been subdivided into multiple lists, each sublist
        representing a cluster.
'''


def separate_into_clusters(data):
    # nodeDict has a key:value relationship of peptide:set(connected proteins). Same effect if you flip peptides with proteins, I just picked peptides first.
    nodeDict = defaultdict(set)
    for x in data:
        nodeDict[x[0]].update([x[1]])

    # number of peptides that have already been included in a cluster and should therefore NOT be considered again.
    usedKeys = set()

    # clusters of peptides specifically. Currently list of list of strings/peptides, will be converted to list of list of tuples/connections at the end (finalClusters variable)
    clusters = []

    # each iteration of this while loop generates a new cluster. Breaks the loop when the number of 'used' peptides equals the amount of total peptides.
    while len(usedKeys) != len(nodeDict):

        # all keys/peptides that don't already belong to a cluster are put into one set
        unusedKeys = set(nodeDict.keys()) - usedKeys

        # one of the keys/peptides is chosen arbitrarily.
        root = list(unusedKeys)[0]

        # proteins that have been considered are initialized with the proteins that connect to the root peptide.
        proteins = nodeDict[root]

        # the cluster is initialized with the root peptide.
        cluster = set([root])

        # the root peptide (the sole entry of 'cluster') is taken out of the unusedKeys variable
        unusedKeys -= cluster

        # This while loop iterates over the remaining peptides/keys until no peptides/proteins are added to the cluster.
        #   This is determined by checking if the number of proteins being considered has changed from what it was previously.
        priorLength = -1
        while priorLength != len(proteins):

            # 'priorLength' variable keeps track of the number of proteins in the previous iteration.
            priorLength = len(proteins)

            # the algorithm then loops over all remaining keys/peptides. If that key/peptide connects to any proteins that have already been added, it is added to the cluster.
            for key in unusedKeys:
                if len(nodeDict[key].intersection(proteins)) > 0:
                    cluster.add(key)
                    proteins.update(nodeDict[key])

            # peptides that were added in the above for loop are then subtracted from the keys still under consideration, as they have been added to the cluster
            unusedKeys -= cluster

        # Once no new peptides/proteins are added to the cluster, the cluster is complete and appended to the variable 'clusters'.
        clusters.append(cluster)
        usedKeys.update(cluster)

    # converts the 'clusters' variable from a list of list of strings/peptides to a list of lists of tuples/connections so the proteins are included.
    finalClusters = [[] for i in clusters]
    for x in data:
        for i in range(len(clusters)):
            if x[0] in clusters[i]:
                finalClusters[i].append(x)
    return finalClusters


# step 4 - reduce
'''
Function: reduce_cluster()
Purpose: Provided a list of connections representing a cluster (one of the sublist outputs of the separate_into_clusters()
            function) this function determines which protein group candidates are considered valid in the cluster. This is
            determined by sorting the proteins by the number of peptide connections, with the highest number being considered
            first. Proteins are then added until all peptides have been accounted for, and all remaining proteins after that
            selection are considered invalid.
        NOTE: block comments describe how to add various scoring parameters, should it be decided to include those in
            the future. Ideas currently include a "uniqueness" score (average number of protein connections among the
            connected peptides) or using the average cosine similarity score of the connected peptides. Also possibly
            including all proteins that are considered "close-knit" with the top scoring protein of the time. "close-knit"
            group is a group of proteins with equal opportunity of selection, specifically with
                a) the same number of connections and
                b) collectively connect to the same peptides.
Parameters:
    'data' - see 'data' parameter description for separate_into_clusters() function.
Returns:
    'proteins' - list of strings. Proteins are formatted as (number of proteins in protein group)/(list of proteins separated
    by '/' characters). Thus, by using the str.split() function and removing the first value you can get a list of the
    individual proteins (list of strings).
'''


def reduce_cluster(data):
    # proteins will be the final return value of this function, containing a list of strings that represent one protein group each
    proteins = []

    # These two variables are used to determine whether or not all peptides in the cluster (totalPeptides) are accounted for by peptides connected to the chosen proteins (usedPeptides)
    usedPeptides = set()
    totalPeptides = set([x[0] for x in data])

    # The 'overallScoreDict' variable is specific to considering the original number of connections between a protein and it's peptides in the cluster.
    #   This is used as a tie-breaker for proteins that have the same number of "current" connections. While the numbers of "current" connections and
    #   "former" connections starts the same, the number of "current" connections decreases as proteins and their connected peptides are removed.
    overallScoreDict = defaultdict(int)
    for x in data:
        overallScoreDict[x[1]] += 1

    # A protein is added to the list for each iteration of this while loop. When the number of peptides connected to the chosen peptides equals the number of total peptides, leave the loop.
    while len(usedPeptides) != len(totalPeptides):
        '''
        if using unique peptide scoring mechanism:
            write uniqueness score for each peptide, write to variable "pepDict"
        '''

        # 'protDict' is a dictionary with key:values of protein:set(connected peptides) relation. Note that because various peptides are excluded with each iteration of the while loop, the number of connected peptides can shrink with each iteration.
        protDict = defaultdict(set)
        for x in data:
            if x[0] not in usedPeptides:
                protDict[x[1]].update(set([x[0]]))

        # presuming the chosen scoring mechanisms can result in two or more top proteins with the same score(s), the proteins are ordered alphabetically as a final tier of consideration.
        #   It is recommended one make the scoring mechanism stringent enough to never use alphabetical choices as a criteria (as that has nothing to do with the actual probability the
        #   protein is found in a sample), but in case the scoring ultimately fails to identify a top candidate, this will make choosing the "first" protein consistent over experiments.
        protDict = OrderedDict(sorted(protDict.items()))
        scoreDict = OrderedDict()
        for key in protDict:
            '''
            if using any peptide scoring mechanism:
                scoreDict[key] = (
                    number of current connections,
                    number of original connections,
                    peptide score ~negative representation value if lower==better~)
            '''
            scoreDict[key] = (len(protDict[key]), overallScoreDict[key])
        sortedKeys = sorted(scoreDict, key=scoreDict.get, reverse=True)
        topKeys = [sortedKeys[0]]
        '''
        if using peptide scoring mechanism AND if the scores of the top two are identical:
            Find members of close-knit group of sortedKeys[0], set topKeys to it (function should return list of length 1 if no close-knit group)
        '''

        # under current conditions, the length of the 'topKeys' list variable will always be 1. However, should close-knit protein groups be included in the future, I write the code to make the inclusion easy.
        #   Basically, this for loop takes the top protein(s), converts it to a string format as described in the function description, and adds it to the protein list while simutaneously excluding their peptides
        #   from future consideration.
        for key in topKeys:
            protein = str(len(key))
            for x in key:
                protein += ('/'+str(x))
            proteins.append(protein)
            usedPeptides.update(protDict[key])
    return proteins


'''
Function: find_valid_proteins()
Purpose: This function takes a list of connections and provides a list of proteins that are considered "valid" using a series
            of IDPicker-based algorithms found in this file.
Parameters:
    'data' - see 'data' parameter description for group_nodes_with_same_edge() function.
Returns:
    dictionary 'finalDict'
        key - string representing a single protein.
        value - string representing the protein group in which you'd find the protein represented by 'key'. See 'proteins'
            return value description in reduce_cluster() function for formatting description.
    NOTE: The final output is returned as a dictionary rather than a straight list because it makes matching protein group
    to peptides in the cosine similarity score output a little easier. An alternative method is including the original
    peptide connections in the output of this function, perhaps in a peptide:set(connected protein groups) depiction. I
    just thought of this format first.
'''


def find_valid_proteins(data):

    # Part 2 of the original paper - grouping proteins/peptides with identical connections together
    data = group_nodes_with_same_edge(data)
    data = group_nodes_with_same_edge(data, first=False)

    # Part 3 of the paper - separating connections by clusters. No connections span between clusters, they are all self-directed.
    cl = separate_into_clusters(data)

    # Part 4 of the paper - reducing the number of proteins under consideration.
    finalProteins = []
    for c in cl:
        finalProteins += reduce_cluster(c)

    # Formatting the final dictionary return value.
    finalDict = {}
    for proteinGroup in finalProteins:
        proteins = proteinGroup.split('/')[1:]
        for protein in proteins:
            finalDict[protein] = proteinGroup
    return finalDict


'''
## Using test data based on the original IDPicker algorithm.
# Part 1 of the original paper - initializing your peptide/protein connections.
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
finalProteins = find_valid_proteins(testData)
print('final proteins: ' + str(finalProteins))
'''
