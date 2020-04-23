"""
autotyp-extract-clustering-counts.py

"""

import sys, re, os, math, ast

def usage():
	print """
#################################################
# USAGE INSTRUCTIONS: autotyp-extract-clustering-counts.py
#################################################

There are no usage notes. If still unsure, consider taking a nap.

"""

def main():
	# Input is the table of Grammatical_markers, with their slot specifications
	gmData = open("gm.langnames.csv", "r")

	# Outputs one record for each marker paradigm (can be more than one per language); for each paradigm we count its slot-sharing cardinalities, and count the total number of slots available
	outfile = open("paradigm-partitions.txt", "w")
	header = "\t".join(["LID", "Language", "Stock", "Category", "SlotsN", "Partition", "HostCat", "Position.binned5", "Behavior.binned4"])
	outfile.write(header+"\n")

	disinctfile = open("distinctive-alignment.txt", "w")
	header = "\t".join(["LID", "Language", "Stock", "SlotArray", "Atr", "P"])
	disinctfile.write(header+"\n")

	markers = {} # markers = {[mkrId][LID][featureSet] = [Number, Person, Role, G, Aditr, S, Atr]}
	mkrCount = 0
	#for line in gmData:
	for line in gmData.read().split('\n'):
		#skip the header
		if re.search(r'^LID', line):
			continue
		line = re.sub('[\n\r]*', '', line) #chomp
		if not len(line.split(','))==10:
			if not line=="":
				print "wrong number of fields in: "+line
			continue
		(LID, AUTOTYPLabel, Exponence, Roles, HostCat, Slot, PosBin, BehaveBin, Language, Stock) = line.split(',')
		#markers[mkrId] = {}
		#exclude markers with HostCat = "S" (i.e. clause attachment), these are clitics and outside the scope of our investigation. In practice only Bilinarra was turning up S marker with slots
		if HostCat=="S":
			continue
		if Slot=="NA":
			continue
		if re.search(r'^V', HostCat):
			HostCat = 'V' # need to substring these because otherwise HostCat like 'V' and 'V;N' are not identified as same
		else:
			continue # easiest to just skip everything else
		#exclude Chintang and Puma! ... because they do not have fixed morphotactic positioning for some markers
		#if LID == "2862" or LID == "2863":
		#	continue

		exponences = Exponence.split(';')
		roless = Roles.split(';')

		# hack for paradigms listed with "NA" in Roles, but where Role is effectively encoded in the AUTOTYP label
		if "NA" in roless and re.search("AGR", AUTOTYPLabel):
			if AUTOTYPLabel=="S/A-AGR":
				roless = ["S", "Atr"]
			elif AUTOTYPLabel=="O-AGR":
				roless = ["P"]
			elif AUTOTYPLabel=="O-AGR.PRO":
				roless = ["O"]
			elif AUTOTYPLabel=="P-AGR":
				roless = ["P"]
			elif AUTOTYPLabel=="S-AGR":
				roless = ["S"]
			elif AUTOTYPLabel=="A-AGR":
				roless = ["Atr"]
			elif AUTOTYPLabel=="G-AGR":
				roless = ["G"]
			elif AUTOTYPLabel=="S/O-AGR":
				roless = ["S", "P"]

		catSet = exponences + roless
		for c in catSet:
			if c in ("S","P","Atr","G","T","Aditr"):
				mkrId = str(mkrCount)+c
				markers[mkrId] = {}
				markers[mkrId]["LID"] = LID
				markers[mkrId]["AUTOTYPLabel"] = AUTOTYPLabel
				markers[mkrId]["Exponence"] = Exponence
				markers[mkrId]["Roles"] = Roles
				markers[mkrId]["HostCat"] = HostCat
				markers[mkrId]["Slot"] = Slot
				markers[mkrId]["PosBin"] = PosBin
				markers[mkrId]["BehaveBin"] = BehaveBin
				markers[mkrId]["Language"] = Language
				markers[mkrId]["Stock"] = Stock
				markers[mkrId]["Category"] = c
		mkrCount += 1


	paradigms = {}
	seenParadigms = []
	for mkrId in markers.keys():
		mkr = markers[mkrId]
		if paradigms.has_key((mkr["LID"], mkr["Category"])):
			continue
		# find other members of the paradigm, i.e. same LID, same category
		altIds = [i for i in markers.keys() if (markers[i]["LID"]==mkr["LID"] and markers[i]["Category"]==mkr["Category"])]
		alternants = {}
		for i in altIds:
			alternants[i] = markers[i]
		#plus the original one is one of the alternants
		alternants[mkrId] = markers[mkrId]
		paradigms[(mkr["LID"], mkr["Category"])] = alternants

		# e.g. paradigms[(59, Atr)] = {marker, marker, marker}

	# Work out number of available slots for each LID / HostCat
	avSlots = {}
	# e.g. avSlots = {(59, V) = (-1, 1, 2), (1439, V) = (1, 2, 5)}
	for mkrId in markers.keys():
		mkr = markers[mkrId]
		if mkr["Slot"]=="NA":
			continue

		slotSpec = mkr["Slot"] # default case where the slotSpec is simply a pos/neg integer
		# ... but there are a bunch of datapoints where the Slot value is something like "3PRS" or "-1SNPST"; these indicate slot position *relative to some other marker*. In all such cases the paradigmatic alternants are all positioned relative to the same marker, so we can essentially just extract the number element here. Also all these cases are instances of absolute category clustering. NOT TRUE
		'''if re.search(r"[^\-0-9]", mkr["Slot"]):
			slotSpec = ""
			match = re.match(r'(\-?[0-9])', mkr["Slot"])
			slotSpec = match.group(1)'''

		if avSlots.has_key((mkr["LID"], mkr["HostCat"])):
			avSlots[(mkr["LID"], mkr["HostCat"])].append(slotSpec)
		else:
			avSlots[(mkr["LID"], mkr["HostCat"])] = [slotSpec]

	## NB this "slot-counter" under-counts the slot array for some templates. In particular, we have noticed some arrays that have obvious gaps like [-1 1 3], and this shows that there could also be unlisted slots at either end of the template. The upshot of this is that, since the true number of slots is higher than what we cound, the true unlikeliness of category clustering is sometimes more extreme than the figures we calculate

	# add implied slots
	for key, slots in avSlots.iteritems():
		extras = []
		for s in slots:
			if not re.search(r"[A-Za-z]", s):
				sn = int(s)
				if sn > 1:
					for n in range(1, sn):
						extras.append(str(n))
				elif sn < 0:
					for n in range(sn, 0):
						extras.append(str(n))
		avSlots[key].extend(extras)

	slotCounts = {}
	for key, slots in avSlots.iteritems():
		uniq = set(slots)
		count = len(uniq)
		slotCounts[key] = count
		#print str(key)
		#print str(uniq)
		#print str(set(slots))

	for key, alternants in paradigms.iteritems():
		#print alternants
		# "alternants" is a set of markers from the same language that have the same Role
		(LID, Category) = key
		firstOne = alternants[alternants.keys()[0]]

		# work out the set partition
		#print alternants
		slotArray = []
		for a in alternants:
			if alternants[a]["Slot"]=="NA":
				continue
			# there are a bunch of datapoints where the Slot value is something like "3PRS" or "-1SNPST"; these indicate slot position *relative to some other marker*. In all such cases the paradigmatic alternants are all positioned relative to the same marker, so we can essentially just extract the number element here. Also all these cases are instances of absolute category clustering.
			'''elif re.search(r"[^\-0-9]", alternants[a]["Slot"]):
				#print str(key) + "\t" + alternants[a]["Slot"]
				slotSpec = ""
				match = re.match(r'(\-?[0-9])', alternants[a]["Slot"])
				slotSpec = match.group(1)
				slotArray.append(slotSpec)
			else:'''
			slotArray.append(alternants[a]["Slot"])
		if not slotArray:
			continue
		# slotArray e.g. ['1', '1', '1', '1', '1', '1', '-2']
		valueCounter = {}
		for s in slotArray:
			if valueCounter.has_key(s):
				valueCounter[s] += 1
			else:
				valueCounter[s] = 1
		partition = []
		for key, count in valueCounter.iteritems():
			partition.append(count)
		partition.sort(reverse=True)
		outString = "\t".join([])

		# look up the slotCount
		if slotCounts.has_key((firstOne["LID"], firstOne["HostCat"])):
			slotsN = slotCounts[(firstOne["LID"], firstOne["HostCat"])]
		else:
			print "Couldn't find slots number for: "+firstOne["LID"]+", "+firstOne["HostCat"]
			continue # this means no slots were found for this paradigm

		#"LID", "AUTOTYPLabel", "SlotsN", "SetCards", "HostCat", "Position.binned5", "Behavior.binned4"
		# output only those that have SlotsN > 1 and sum(SetCards) > 1; otherwise there clustering probability is uninteresting because it can only be fully clustered
		# now including single-slot verbs, though they won't be included in stats test later
		#if (slotsN > 1 and sum(partition) > 1):
		if (sum(partition) > 1):
			outString = "\t".join([firstOne["LID"], firstOne["Language"],firstOne["Stock"],firstOne["Category"], str(slotsN), str(partition), firstOne["HostCat"], firstOne["PosBin"], firstOne["BehaveBin"],])
			outfile.write(outString+"\n")

		# if this is an Atr paradigm, and there is also a P paradigm for the same language, then print their slot matrices
		if (Category=='Atr' and paradigms.has_key((LID, "P"))):
			#print firstOne["Language"]
			#print paradigms[(LID, "P")]
			#now we need to know what goes in what slot
			slots = sorted(list(set(avSlots[(LID, "V")])))
			AtrHits = []
			PHits = []
			for s in slots:
				# how many Atr markers in this slot?
				hitIds = [i for i in alternants if alternants[i]["Slot"]==s]
				hitcount = len(hitIds)
				AtrHits.append(hitcount)

				pAlternants = paradigms[(LID, "P")]
				hitIds = [i for i in pAlternants if pAlternants[i]["Slot"]==s]
				hitcount = len(hitIds)
				PHits.append(hitcount)
			AtrHitsStr = [str(i) for i in AtrHits]
			PHitsStr = [str(i) for i in PHits]
			if slotsN > 1: # we're excluding verbs with only one slot
				distring = "\t".join([LID, firstOne["Language"], firstOne["Stock"], ",".join(slots), ",".join(AtrHitsStr), ",".join(PHitsStr)])
				disinctfile.write(distring+"\n")

if __name__ == "__main__":
	main()