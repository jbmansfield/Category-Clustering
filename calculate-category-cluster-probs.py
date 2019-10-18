"""
calculate-category-cluster-probs.py

"""

import sys, re, os, math, ast

def usage():
	print """
#################################################
# USAGE INSTRUCTIONS: calculate-category-cluster-probs.py
#################################################


"""

outfilePartitions = open("paradigm-partition-probs.txt", "w")
outfilePartitions.write("LID\tLanguage\tStock\tHOST.CAT\tCategory\tPosition.binned5\tSLOTS\tSLOTS.USED\tPARTITION\tENTROPY\tPROB\tCLUST.INDEX\n")
outfileReadable = open("clustering-probabilities.readable.txt", "w")
possiblePartsIndex = open("possible-partitions.txt", "w")
possiblePartsIndex.write("MORPHS\tSLOTS\tPARTITION\tENTROPY\tPROB\n")

def main():

	# Working through a set of typological data on inflectional paradigms
	clusteringData = open("paradigm-partitions.txt", "r")
	lines = clusteringData.readlines()

	# Read the input data into this multidimensional dictionary
	occurrences = {} # occurrences[morphsN][slotsN][partition] = [(LID, HostCat), (LID, HostCat), (LID, HostCat)]
	for line in lines:
		#skip the header
		if re.search(r'^LID', line):
			continue
		line = re.sub('[\n\r]*', '', line) #chomp
		(LID, Language, Stock, Category, slotsStr, partitionStr, HostCat, PosBin, BehaveBin) = line.split('\t')
		slotsN = int(slotsStr)
		partition = ast.literal_eval(partitionStr)
		morphsN = sum(partition)
		if not occurrences.has_key(morphsN):
			occurrences[morphsN] = {}
		if not occurrences[morphsN].has_key(slotsN):
			occurrences[morphsN][slotsN] = {}
		if not occurrences[morphsN][slotsN].has_key(partitionStr):
			occurrences[morphsN][slotsN][partitionStr] = [(LID, Language, Stock, HostCat, Category, PosBin)]
		else:
			occurrences[morphsN][slotsN][partitionStr].append((LID, Language, Stock, HostCat, Category, PosBin))

	for morphsN in occurrences.keys():
		if morphsN==1:
			continue #uninformative
		for slotsN in occurrences[morphsN].keys():
			#if slotsN==1:
			#	continue #uninformative
			outfileReadable.write("##############\nmorphsN: "+str(morphsN)+"\n")
			outfileReadable.write("slotsN: "+str(slotsN)+"\n")

			# First we work out the probability of every *possible* partition, and its entropy. Thus when we get on to actually observed partitions in the next bit, we can place our observation in the distribution of outcome likelihoods, and establish whether it is in the extremes of this distribution
			possPartitions = {} # possPartitions[str(partition)] = (entropy, probability)
			sys.setrecursionlimit(5000) # higher than standard limit needed for some recursive partitioning
			grossPartitions = partition_int(morphsN)
			# remove any cardSets that have more parts than the number of slots available
			partitions = []
			for p in grossPartitions:
				if len(p) <= slotsN:
					partitions.append(p)
			for p in partitions:
				entropy = get_entropy(p)
				(partitionsPermutations, outcomePermutations) = calculate_partition_prob(p, slotsN)
				prob = float(partitionsPermutations) / outcomePermutations
				possPartitions[str(p)] = (entropy, prob)
				possiblePartsIndex.write(str(morphsN)+"\t"+str(slotsN)+"\t"+str(p)+"\t"+str(entropy)+"\t"+str(prob)+"\n")

			trialsN = 0 # the total number of observed paradigms with this many markers and this many slots
			for p in occurrences[morphsN][slotsN].keys():
				trialsN += len(occurrences[morphsN][slotsN][p])
			outfileReadable.write("trialsN: "+str(trialsN)+"\n\n")
			for partitionStr in occurrences[morphsN][slotsN].keys():
				occurrencesN = len(occurrences[morphsN][slotsN][partitionStr])
				outfileReadable.write("partition: "+partitionStr+"\n")
				outfileReadable.write("occurrencesN: "+str(occurrencesN)+"\n")

				# sum the probs with this config and all those that equal or greater entropy, i.e. more dispersion; and all those that have equal or less entropy, i.e. more clustering
				tailProb = None # decided not to use this for the moment
				dispersionProb = 0
				clusteringProb = 0 # read this as: "probability of observing this much clustering or more"
				entropy = get_entropy(ast.literal_eval(partitionStr))
				for ps, vals in possPartitions.iteritems():
					#print str(vals[1])
					if vals[0] > entropy: # i.e. all those partitions with greater entropy
						dispersionProb += vals[1] # get counted towards the prob of greater dispersion
					elif vals[0] < entropy:
						clusteringProb += vals[1] # conversely these count towards prob of more clustering
					else: #has same entropy, so counts towards both; (NB this means that sum of dispersionProb and clusteringProb is more than 1; they have overlap)
						dispersionProb += vals[1]
						clusteringProb += vals[1]
				outfileReadable.write("disp: "+str(dispersionProb)+"\n")
				outfileReadable.write("clust: "+str(clusteringProb)+"\n")
				partitionType = "unknown"
				if clusteringProb > 0.5 and dispersionProb > 0.5: # tails in each direction are over 50%, i.e. we're in medial territory
					partitionType = "Medial"
					tailProb = 1 # i.e. tails in both directions, covers all possiblities
				elif clusteringProb <= 0.5: # on the clustering side
					partitionType = "Clustered"
					tailProb = clusteringProb	# decided not to use this for the moment
				else:
					partitionType = "Dispersed" # only possibility left is dispersion
					tailProb = dispersionProb

				clustIndex = 1 - clusteringProb

				slotsUsed = len(ast.literal_eval(partitionStr))
				
				for paradigm in occurrences[morphsN][slotsN][partitionStr]:
					# outfilePartitions.write("LID\tHOST.CAT\tSLOTS\tPARTITION\tENTROPY\tPROB\tPROB.CUMUL\n")
					outfilePartitions.write(paradigm[0] + "\t" + paradigm[1] + "\t" + paradigm[2] + "\t" + paradigm[3] + "\t" + paradigm[4] + "\t" + paradigm[5] + "\t" + str(slotsN) + "\t" + str(slotsUsed) + "\t" + partitionStr + "\t" + str(entropy) + "\t" + str(possPartitions[partitionStr][1]) + "\t" + str(clustIndex) + "\n")

				(partitionsPermutations, outcomePermutations) = calculate_partition_prob(ast.literal_eval(partitionStr), slotsN)
				prob = float(partitionsPermutations) / outcomePermutations
				outfileReadable.write("Prob of occurring in a single trial: "+str(prob)+"\n")


# given a partition and a number of slots, this calculates how many ways that partition can occur, and how many possible outcomes there are, together giving the probability of the partition
def calculate_partition_prob (partition, slotsN):
	morphsN = sum(partition) # i.e. number of morphs
	groupsN = len(partition) # i.e. number of *non-empty* groups

	#number of sets with cardinality 1..morphsN
	partitionCount = {}
	for n in range(1, morphsN):
		#how many with group.card = n in partition?
		partitionCount[n] = partition.count(n)

	# Number of ways that a partition config can arise: groupConfigsN
	grossConfigs = 1
	remainingMorphs = morphsN
	for group in partition:
		# select group from remainingMoprhs
		thisGroupConfigs = nCk(remainingMorphs, group)
		grossConfigs *= thisGroupConfigs
		remainingMorphs -= group
	#print "gross " + str(grossConfigs)
	identPermutationsN = 1
	for c in partitionCount.keys():
		if partitionCount[c]>0:
			#identPermutationsN *= partitionCount[c]
			identPermutationsN *= math.factorial(partitionCount[c])

	groupConfigsN = grossConfigs / identPermutationsN
	#print "setConfigsN: " + str(setConfigsN)

	# Number of ways the groups can be arranged in the slots: slotConfigsN
	groupPermutations = math.factorial(groupsN) # how many ways can the sets permute?
	slotSelections = nCk(slotsN, groupsN) # how many ways can they select slots
	slotConfigsN = groupPermutations * slotSelections
	#print "slotConfigsN: " + str(slotConfigsN)

	# Number of permutations for this partition configuration (in one trial)
	partitionsPermutations = groupConfigsN * slotConfigsN
	#print "partitionsPermutations: " + str(partitionsPermutations)

	# Number of possible outcomes
	outcomePermutations = slotsN ** morphsN
	#print "outcomePermutations: " + str(outcomePermutations)
	
	return (partitionsPermutations, outcomePermutations)

# partitioner function by David Eppstein, http://code.activestate.com/recipes/218332/
# version here is a "variation" added by George Yoshida, which returns partitions like [4,2,1] as opposed to [1,2,4] in the original version
def partition_int(n):
    # reverse order
    if n == 0:
        yield []
        return

    for p in partition_int(n-1):
        yield p + [1]
        if p and (len(p) < 2 or p[-2] > p[-1]):
            yield p[:-1] + [p[-1] + 1]

def get_entropy(integers):
	# takes a set of integers as input, e.g. [5, 3, 1, 1]
	entropy = 0
	for i in integers:
		prob = float(i) / sum(integers)
		surprisal = -math.log(prob,2)*prob
		#print "s:"+str(surprisal)
		entropy += surprisal
	return entropy

# this neat "n choose k" function is by Nas Banov
# https://stackoverflow.com/questions/3025162/statistics-combinations-in-python
from operator import mul
from fractions import Fraction
def nCk(n,k): 
	return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

if __name__ == "__main__":
    main()