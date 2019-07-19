import pyodbc as db
import random
import pandas as pd
import numpy as np
import math
import sys
import statistics as stat
import c03_globalVariables as glob

filePopulationHistory = open("C:\\Users\\u374441\\Danfoss\\PS WF - MSA - Documents\\Projects\\2019 03 - Changeover\\populationHistory.txt", "w", encoding="utf-8")


# sort array by fitness
def udf_sortByFitness(lFitness):
	''' sorts a povided fitness array by fitness for each member. 
	It further assigns a selection probability based on the distance to the worst member.
	Further it calculates a running sum and sorts accordingly.
	E.g. the further a member is away from the population's worst the higher the chances of later selection will be.
	
	INPUT: Fitness Array, with minimum colums [member], [fitness]
	SIDE EFFECTS: none	
	RETURNS: sorted fitness array with columns [selectionProb_cumSum],[member]  '''

	dfFitness = pd.DataFrame(lFitness, columns=['member','fitness']) # load member and calc. fitness into dataFrame
	dfFitness['distanceToWorst'] = dfFitness['fitness'] - (np.nanmax(dfFitness['fitness'].values)+1) # calculate distance to worst fitness
	dfFitness['selectionProb'] = dfFitness['distanceToWorst'] / dfFitness['distanceToWorst'].sum() # calculate selection probability as % of distanceToWorst
	dfFitness['selectionProb'] = dfFitness['selectionProb'].abs() # ensure it is a positive number
	dfFitness = dfFitness.sort_values('selectionProb', ascending=False) # sort descing by selection probability
	dfFitness['selectionProb_cumSum'] = dfFitness['selectionProb'].cumsum(axis=0) # create cummulative sum of selection probability
	dfFitness=dfFitness.round({'selectionProb':5, 'selectionProb_cumSum':5}) # round probability
	lFitness_sorted = dfFitness[['selectionProb_cumSum','member']].values.tolist() # convert dataFrame to regular python list

	return lFitness_sorted


# ACTIVE calculate fitness of a population
def udf_calcFitness3(dPopulation, dWcList, dMaterialFamily, dTimeMatrix, dMaterialCO,lMinFitness, dMachineConfig):
	
	'''
	INPUT:
	:param dPopulation:	 		>dict; population
	:param dWcList: 			>dict; with orders and material numbers
	:param dMaterialFamily:		>dict; mapping of material to family index
	:param dTimeMatrix: 		>dict; family change over times
	:param dMaterialCO: 		>dict; material change over times
	:param lMinFitness: 		>list; minimum fitness recorded in each run
	:param dMachineConfig: 		>dict; illegal machine configuration
	
	SIDE EFFECTS: 
	none

	RETURN: 
	:return lFitness: 			>list; with all member fitnesses
	:return dMembers: 			>dict; for members, genome, name and fitness
	:return lMinFitness:		>list; minimum fitness per run
	:return fMinFitness_run:	>float; minimum fitness for this run
	:return fIllegalPerc: 		>float; illegal percentage for this run
	
	SUMMARY:
	> Calculates the fitness per member of a population/array based on changeover times between materials (default to family).
	Applys penalty terms for illegal configurations and uneven population distribution
	Loops (outer to inner):
	0: member in population
	1: machines in member, as per break points
	2: genes in Machine in Member



	'''
	
	# initialize all values
	sMaterial1 = ''
	sMaterial2 = ''
	sFamily1 = ''
	sFamily2 = ''
	sQuantity1 = 0
	sQuantity2 = 0
	sCycleTime1 = 0
	sCycleTime2 = 0
	sChangeovertime = 0
	fFitness = 0
	lFitness = []
	sMemberName = ''
	dMembers = {}
	fMinFitness_run = 100000000
	iTotalRuns = 0.0
	iIllegalRuns = 0.0
	fIllegalPerc = 0.0

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< MEMBER LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for index0, (key,member) in enumerate(dPopulation.items()) :	# for every member in a population


		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< MACHINE LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++

		iPreviousBreak = 0
		fFitness = 0
		iCountIndex1 = 0
		fFitnessBalance = []
		iMemberRuns = 0
		bisIllegal = False

		iTotalRuns +=1
		iIllegalConfigMultiplier = 1

		for k in range(0, glob.iNumberMachines):
			
			if k == glob.iNumberMachines-1:
				iNextBreak = len(member["genome"])

			else:
				iNextBreak = member["breaker"][k]

			fFitnessM = 0 # reset fitness
			lGenomeW0 = [x for x in member["genome"][iPreviousBreak:iNextBreak] if x != 0]

			#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< GENE LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++
			for index1, gene in enumerate(lGenomeW0):		# for every gene in a member

				iMemberRuns +=1

				if index1 < len(lGenomeW0)-1:				# if not the last material >> last material does not have changeover Times: BUT has processing times--> adj. needed
					sMaterial1 = dWcList[gene]['material']	# get family information
					sMaterial2 = dWcList[lGenomeW0[index1+1]]['material']
					sQuantity1 = dWcList[gene]['quantity']
					sQuantity2 = dWcList[lGenomeW0[index1+1]]['quantity']
					sFamily1 = dMaterialFamily[sMaterial1]['family']
					sFamily2 = dMaterialFamily[sMaterial2]['family']
					sCycleTime1 = dMaterialFamily[sMaterial1]['cycleTime']
					sCycleTime2 = dMaterialFamily[sMaterial2]['cycleTime']

					# take material changeover times if available, otherwise family, otherwise error
					if (str(sMaterial1)+"-"+str(sMaterial2)) in dMaterialCO:
						sChangeovertime = dMaterialCO[str(sMaterial1)+"-"+str(sMaterial2)]
					elif (str(sFamily1)+"-"+str(sFamily2)) in dTimeMatrix:
						sChangeovertime = dTimeMatrix[str(sFamily1)+"-"+str(sFamily2)]
					else: print("ERROR")

					fFitnessM += (sQuantity1*sCycleTime1)+(sChangeovertime) #calculate overall fitness for every pair
					
					# if machineNbr is in the MachineConfig (contains illegal configs), add a penalty term and set illegal flag
					if (k+1) in dMachineConfig[sFamily1] or (k+1) in dMachineConfig[sFamily2]:
						iIllegalConfigMultiplier += 50
						bisIllegal = True

				# for last member of the machine
				else: 
					sMaterial1 = dWcList[gene]['material']	
					sQuantity1 = dWcList[gene]['quantity']
					sCycleTime1 = dMaterialFamily[sMaterial1]['cycleTime']
					sFamily1 = dMaterialFamily[sMaterial1]['family']
					fFitnessM += (sQuantity1*sCycleTime1)  

					if (k+1) in dMachineConfig[sFamily1]:
						iIllegalConfigMultiplier += 50
						bisIllegal = True

				iCountIndex1 += 1


			fFitness += fFitnessM
			fFitnessBalance.append(fFitnessM)
			iPreviousBreak = iNextBreak

		# add stDev of distribution as penalty term > the more uneven a population is, the higher the penalty
		fFitness = (fFitness+stat.stdev(fFitnessBalance))*iIllegalConfigMultiplier

		if bisIllegal == True:
			iIllegalRuns +=1
			
		else: 
			if glob.bDebug1 == True: print("--------------------------------------------new member: ", key)
			if glob.bDebug1 == True: print("Fitness: ", fFitness, iIllegalConfigMultiplier)
			if glob.bDebug1 == True: print("Others: ", fFitnessBalance)

		sMemberName = key # set the memberName based on previous input and iterations
		lFitness.append([sMemberName, fFitness]) # create fitness array

		#set minimum fitness if it is lower than the previous fitness
		if fFitness < lMinFitness[0]:
			lMinFitness[0] = fFitness
			lMinFitness[1] = sMemberName
			lMinFitness[2] = member["genome"]
			lMinFitness[3] = member["breaker"]
			lMinFitness[4] = "Generation: "+str(glob.iGenerationCount-1)

		if fFitness < fMinFitness_run:
			fMinFitness_run = fFitness

		# create member array with name, fitness and genome
		dMembers[sMemberName] = {}
		dMembers[sMemberName]['fitness'] = fFitness
		dMembers[sMemberName]['genome'] = member["genome"]
		dMembers[sMemberName]["illegal"] = iIllegalConfigMultiplier

	return lFitness, dMembers, lMinFitness, fMinFitness_run, fIllegalPerc

# select best parents from pool
def udf_selectParentsFromPool(dMembers, lFitness_sorted, dPopulation):
	
	'''
	INPUT:
	:param dMembers:			>dict; Member array with fitness, genome, name
	:param lFitness_sorted:		>list; sorted fitness list
	:param dPopulation:			>dict; active population array

	SIDE EFFECTS:
	none

	RETURNS:
	lPopulation_new, lPopulation_new_names, dPopulation_new
	:return lPopulation_new:		>list; genomes of newly seleted population
	:return lPopulation_new_name:	>list; parallel array with member names
	:return dPopulation_new:		>dict; array of new population !may contain duplicates

	SUMMARY:
	Selects parents from population randomly based on sorted fitness list (distance to worst)
	!KingPrevails: guarantess that the fittest individual is selected as a parent 
	'''
	
	lPopulation_new = []
	lPopulation_new_names = []
	dPopulation_new={}
	iStartParents = 0

	# guarantees that the fittest member survives
	if glob.bKingPrevails == True:
		dPopulation_new[0]={}
		dPopulation_new[0]["member"]=lFitness_sorted[0][1]
		dPopulation_new[0]["genome"]=dMembers[lFitness_sorted[0][1]]['genome']
		dPopulation_new[0]["fitness"]=dMembers[lFitness_sorted[0][1]]['fitness']
		dPopulation_new[0]["breaker"]=dPopulation[lFitness_sorted[0][1]]["breaker"]
		iStartParents = 1


	for i in range(iStartParents, glob.limPopulationSize):	#perform POLPULATIONSIZE number of iterations
		fRand = random.uniform(0.0, 1.0)	# create random number

		for member in lFitness_sorted:		# if RAND is smaller than the cummulative sum, select the member into the pool of Parents
			if fRand <= member[0]:
				lPopulation_new.append(dMembers[member[1]]['genome'])
				lPopulation_new_names.append(member[1])
				dPopulation_new[i]={}
				dPopulation_new[i]["member"]=member[1]
				dPopulation_new[i]["genome"]=dMembers[member[1]]['genome']
				dPopulation_new[i]["fitness"]=dMembers[member[1]]['fitness']
				dPopulation_new[i]["breaker"]=dPopulation[member[1]]["breaker"]
				break

	return lPopulation_new, lPopulation_new_names, dPopulation_new

# perform [p]artially [m]apped [x]CrossOver
def udf_matingPMX(lPopulation_new, iChildCounter, lPopulation_new_names, dPopulation_new, dMembers, fMutationRate):

	'''
	! DOES NOT WORK WITH CURRENT GENOME SETUP ! 
	INPUT:
	:param lPopulation_new:			>list; new population (selected parents)
	:param iChildCounter:			>int; member naming counter
	:param lPopulation_new_names:	>list; parallel array with member names
	:param dPopulation_new:			>dict; new population (selected parents)
	:param dMembers:				>dict; full members array
	:param fMutationRate:			>float; mutation rate

	SIDE EFFECTS:
	none

	RETURNS: return lPopulation_offspring, iChildCounter, lPopulation_offspring_names, dPopulation_offspring
	:return iChildCounter: 			>int; member naming counter
	:return dPopulation_offspring:	>dict; newly created offspring

	SUMMARY:
	'''

	lPopulation_offspring = []
	lPopulation_offspring_names = []
	sParents= ""
	dPopulation_offspring={}

	for index_m, mother in enumerate(lPopulation_new[::2]):	# only choose every second member of the array (first is mother, second is father)

		fRandXO1 = random.randint(0,math.floor(len(mother)/2))			# create two crossover points randomly
		fRandXO2 = random.randint(math.floor(len(mother)/2),len(mother))
		
		# since only every 2nd item, index for skipped item needs to be calculated
		iFather = (index_m)+(index_m)+1
		iMother = (index_m)+(index_m)
		lChild1 = []
		lChild2 = []
		lMapMother = []
		lMapFather = []
		sParents= "("+lPopulation_new_names[iMother] + ","+lPopulation_new_names[iFather]+")"
		fFitnessMother = dMembers[dPopulation_new[iMother]["member"]]["fitness"]
		fFitnessFather = dMembers[dPopulation_new[iFather]["member"]]["fitness"]


		lMotherBreak = dPopulation_new[iMother]["breaker"]
		lFatherBreak = dPopulation_new[iFather]["breaker"]
		lChild1Break = []
		lChild2Break = []

		# catch if rands are the same and adjust to be at least 1 apart
		if fRandXO1 == fRandXO2 and fRandXO1 > 1:
			fRandXO1 = fRandXO1 -1
		elif fRandXO1 == fRandXO2 and fRandXO1 <=1:
			fRandXO2 = fRandXO2 + 1

		# create mapping arrays
		lMapMother.extend(mother[fRandXO1:fRandXO2])
		lMapFather.extend(lPopulation_new[iFather][fRandXO1:fRandXO2])


		for index_gM, gene in enumerate(mother):	# iterate over mother genes for the first child


			if index_gM < fRandXO1 or index_gM >= fRandXO2:	# if the index is outside the crossover zone, perform mapping

				if gene not in lMapFather:	# check if the gene is in the map, if not: append
					lChild1.append(gene)

				else: 
					bGeneFound = False
					iMapIndex = 0
					sOppParentGene = gene

					# gene search:
					# perform as long as not found:
					# check what gene maps to in the fathers map
					# check if new gene is still in fathers map
					# if not, append. If yes, repeat. 

					while bGeneFound == False:
						if sOppParentGene == lMapFather[iMapIndex]:	
							sOppParentGene = lMapMother[iMapIndex]

							if sOppParentGene not in lMapFather:
								bGeneFound = True
								lChild1.append(sOppParentGene)
								break
							else:
								iMapIndex = 0
						else:
							iMapIndex += 1

			else:
				lChild1.append(lPopulation_new[iFather][index_gM])


		# same as for mother, but inverse
		geneFather = lPopulation_new[iFather]

		for index_gF, gene in enumerate(geneFather):

			if index_gF < fRandXO1 or index_gF >= fRandXO2:

				if gene not in lMapMother:
					lChild2.append(gene)
				else: 
					bGeneFound = False
					iMapIndex = 0
					sOppParentGene = gene

					while bGeneFound == False:

						if sOppParentGene == lMapMother[iMapIndex]:
							sOppParentGene = lMapFather[iMapIndex]

							if sOppParentGene not in lMapMother:
								bGeneFound = True

								lChild2.append(sOppParentGene)
								break
							else:

								iMapIndex = 0
						else:
							iMapIndex += 1
						
			else:

				lChild2.append(mother[index_gF])

		# creating new offspring arrays
		iChildCounter += 1
		lPopulation_offspring.append(lChild1)
		lPopulation_offspring_names.append("child"+str(iChildCounter))
		dPopulation_offspring["child"+str(iChildCounter)]={}
		dPopulation_offspring["child"+str(iChildCounter)]["genome"]=lChild1
		dPopulation_offspring["child"+str(iChildCounter)]["breaker"]=lChild1Break

		iChildCounter += 1
		lPopulation_offspring.append(lChild2)
		lPopulation_offspring_names.append("child"+str(iChildCounter))
		dPopulation_offspring["child"+str(iChildCounter)]={}
		dPopulation_offspring["child"+str(iChildCounter)]["genome"]=lChild2
		dPopulation_offspring["child"+str(iChildCounter)]["breaker"]=lChild2Break



	return lPopulation_offspring, iChildCounter, lPopulation_offspring_names, dPopulation_offspring

#mutation by swaping
def udf_mutateSwap(fMutationRate,lPopulation_offspring, dPopulation_offspring):
	lchild_mutate = []
	lPopulation_offspringMutated=[]

	for j,child in dPopulation_offspring.items():		#iterate over all children
		if random.uniform(0.0, 1.0) < fMutationRate:		#calculate random number; if lower than mutation rate, MUTATE
			fRandMutate1 = random.randint(0,math.floor(len(child["genome"])/2))	#calculate two random mutation points
			fRandMutate2 = random.randint(math.floor(len(child["genome"])/2),len(child["genome"])-1)

			#check if mutation points are not the same
			if fRandMutate1 == fRandMutate2 and fRandMutate1 > 1:
				fRandMutate1 = fRandMutate1 -1
			elif fRandMutate1 == fRandMutate2 and fRandMutate1 <=1:
				fRandMutate2 = fRandMutate2 + 1

			#exchange mutated genes

			lchild_mutate = dPopulation_offspring[j]["genome"]
			
			lchild_mutate[fRandMutate1], lchild_mutate[fRandMutate2] = lchild_mutate[fRandMutate2], lchild_mutate[fRandMutate1]
			dPopulation_offspring[j]["genome"] = lchild_mutate



		lPopulation_offspringMutated.append(lchild_mutate)

	return lPopulation_offspringMutated, dPopulation_offspring


def udf_mutateRegen(lList,iGeneLen):
	lchild_mutate = []
	for i,item in enumerate(lList):		#iterate over all children
		if random.uniform(0.0, 1.0) < glob.fMutationRate/4:		#calculate random number; if lower than mutation rate, MUTATE
		#	print("MUTATING Breakers")

			iNewBreak = random.randint(1,iGeneLen)

			while iNewBreak in lList:
				 iNewBreak = random.randint(1,iGeneLen)
				 #print("newBreaker already in list, searched new")

			#print(lList)
			lList[i] = iNewBreak
			#print(lList)

	return	lList


def udf_cataclysm(dPopulation,lGenome):

	iRunIndex = 0
	for key, member in dPopulation.items():

		if random.uniform(0.0,1.0) < glob.iDeletionProb:
			glob.iChildCounter +=1
			sNewKey = "member"+str(glob.iChildCounter)
			dPopulation[sNewKey] = dPopulation.pop(key)
			dPopulation[sNewKey]["genome"], dPopulation[sNewKey]["breaker"] = udf_makeNewMember(lGenome)
			

	return dPopulation


def udf_makeNewMember(lGenome_0):

	
	#print(i)
	lNewMember = []
	lBreakGenome = []
	#fill member randomly based on available genes - reduce available genes
	
	lEmptyAppend = [i*0 for i in range(0, (glob.iNumberMachines-1)*len(lGenome_0))]

	lGenome = lGenome_0+lEmptyAppend
	lGenesAvailable = lGenome.copy()

	#fill member randomly based on available genes - reduce available genes
	for gene in range(0, len(lGenome)):
		sOrder= random.choice(lGenesAvailable)
		lNewMember.append(sOrder)
		lGenesAvailable.remove(sOrder)


	for j in range(1, glob.iNumberMachines):
		
		iRandBreaker = (j*len(lGenome_0))
		lBreakGenome.append(iRandBreaker)

	lBreakGenome.sort()

	return lNewMember, lBreakGenome


def udf_cloneMutate(lPopulation_new, lPopulation_new_names, dPopulation_new, dMembers, dMaterialFamily, dMachineConfig, dWcList, lGenome):
	#print(lPopulation_new)
	lPopulation_offspring = []
	lPopulation_offspring_names = []
	sParents= ""
	dPopulation_offspring={}

	lLoopList = [x["genome"] for i,x in dPopulation_new.items()]
	#print(lLoopList)

	for index_m, mother in enumerate(lLoopList[::2]):	#only choose every second member of the array (first is mother, second is father)
		#print("New Pairing")
		fRandXO1 = random.randint(0,math.floor(len(mother)/2))			#create two crossover points randomly
		fRandXO2 = random.randint(math.floor(len(mother)/2),len(mother))
	
		#fRandXO1 = 3
		#fRandXO2 = 6	
		
		# since only every 2nd item, index for skipped item needs to be calculated
		iFather = (index_m)+(index_m)+1
		iMother = (index_m)+(index_m)
		lChild1 = []
		lChild2 = []
		geneFather = dPopulation_new[iFather]["genome"][:]
		#sParents= "("+lPopulation_new_names[iMother] + ","+lPopulation_new_names[iFather]+")"
		fFitnessMother = dMembers[dPopulation_new[iMother]["member"]]["fitness"]
		fFitnessFather = dMembers[dPopulation_new[iFather]["member"]]["fitness"]


		#print("fitness (Mum, Dad): ", fFitnessMother, fFitnessFather)
		#print("breakers (Mum, Dad): ", lMotherBreak, lFatherBreak)

		#catch if rands are the same and adjust to be at least 1 apart
		if fRandXO1 == fRandXO2 and fRandXO1 > 1:
			fRandXO1 = fRandXO1 -1
		elif fRandXO1 == fRandXO2 and fRandXO1 <=1:
			fRandXO2 = fRandXO2 + 1

		lGenome_dom = []
		lGenome_sub = []

		lMotherBreak = dPopulation_new[iMother]["breaker"][:]
		lFatherBreak = dPopulation_new[iFather]["breaker"][:]
		lChild1Break = []
		lChild2Break = []
		lBreak_dom = []
		lBreak_sub = []
		bDominantMother = "" 

		#print("MOTHER: ", mother, " | ", lMotherBreak)
		#print("FATHER: ", geneFather, " | ", lFatherBreak)


		if fFitnessMother > fFitnessFather:
			bDominantMother = "mother"
			lGenome_dom = mother[:]
			lGenome_sub = geneFather[:]
			lBreak_dom = lMotherBreak[:]
			lBreak_sub = lFatherBreak[:]
			fCloneDom_cut = (fFitnessMother)/(fFitnessMother+fFitnessFather)

			#print("Mother is Dominant, XO: ", fXObreaker_cut)

		elif fFitnessMother < fFitnessFather:
			bDominantMother = "father"
			lGenome_dom = geneFather[:]
			lGenome_sub = mother[:]
			lBreak_dom = lFatherBreak[:]
			lBreak_sub = lMotherBreak[:]
			fCloneDom_cut = (fFitnessFather)/(fFitnessFather+fFitnessMother)
			#print("Father is Dominant, XO: ", fXObreaker_cut)
		else: 
			bDominantMother = "none"
			fCloneDom_cut = 0.5
			lGenome_dom = geneFather[:]	 # set, but irrelevant
			lGenome_sub = mother[:]
			lBreak_dom = lFatherBreak[:]
			lBreak_sub = lMotherBreak[:]

		#child1
		if random.uniform(0.0, 1.0) <= fCloneDom_cut:
			lChild1 = lGenome_dom[:]
			lChild1Break = lBreak_dom[:]
		else: 
			lChild1 = lGenome_sub[:]
			lChild1Break = lBreak_sub[:]


		#child2
		if random.uniform(0.0, 1.0) <= fCloneDom_cut:
			lChild2 = lGenome_dom[:]
			lChild2Break = lBreak_dom[:]
		else: 
			lChild2 = lGenome_sub[:]
			lChild2Break = lBreak_sub[:]
		#print("MOTHER: ", mother)
		#print("Child1: ", lChild1)

		lChild1Break.sort()
		lChild2Break.sort()

		#print("Child1_UM: ", lChild1, " | ", lChild1Break)
		#print("Child2_UM: ", lChild2, " | ", lChild2Break)
		#print("Child1_U_Ill: ", udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList))
		#print("------------------CHILD", index_m,": ",lChild1)
		#print("Genome BEFORE: ", lChild1)
		#print("Illegals BEFORE: ", udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList))



		#>>>>>>>>>>>>>>>>>>>>>>> ENSURE LEGAL SOLUTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<
		if glob.bCorrectChild == True:
			

			for iC1, gene in enumerate(lChild1):	# iterate over mother genes for the first child

				iForceAllocation = 0
				lAllowedMachines = []

				if lChild1[iC1] == 0:
					continue

				sMaterial1 = dWcList[lChild1[iC1]]['material']	# get family information
				sFamily1 = dMaterialFamily[sMaterial1]['family']
				#print(iC1, lChild1[iC1], lChild1Break ,dMachineConfig[sFamily1])

				iMyMachine = min([i  if iC1 < k else len(lChild1Break) for i,k in enumerate(lChild1Break)])+1 #I love list comprehension!

				for k in range(1, glob.iNumberMachines+1):
					if k not in dMachineConfig[sFamily1]:
						lAllowedMachines.append(k)

				#print("Child1 Current Machine, Allowed: ", iMyMachine, lAllowedMachines)

				lAllowedPositions = []
				iUpperBound = 0
				iLowerBound = 0

				for i,m in enumerate(lAllowedMachines):
					#set for First Machine
					if m == 1:
						iLowerBound = 0
						iUpperBound = lChild1Break[m-1]
					elif m==glob.iNumberMachines:
						iLowerBound = lChild1Break[m-2]
						iUpperBound = len(lChild1)

					else:
						iLowerBound = lChild1Break[m-2]
						iUpperBound = lChild1Break[m-1]


					lAllowedPositions.extend(list(range(iLowerBound, iUpperBound)))
					#set for last Machine



				#print("ALLOWED POSITIONS: ", lAllowedPositions)
				if iMyMachine not in lAllowedMachines:
					iForceAllocation = glob.iForceAllocation_G
				#print(glob.fMutationRate)
				if random.uniform(0.0, 1.0) < max(glob.fMutationRate, iForceAllocation):
					fRandXO2 = random.choice(lAllowedPositions)
					lChild1[iC1], lChild1[fRandXO2] = lChild1[fRandXO2], lChild1[iC1]


			for iC2, gene2 in enumerate(lChild2):	# iterate over mother genes for the first child

				iForceAllocation = 0
				lAllowedMachines = []


				if lChild1[iC2] == 0:
					continue

				sMaterial1 = dWcList[lChild1[iC2]]['material']	# get family information
				sFamily1 = dMaterialFamily[sMaterial1]['family']
				#print(iC1, lChild1[iC1], lChild1Break ,dMachineConfig[sFamily1])

				iMyMachine = min([i  if iC2 < k else len(lChild2Break) for i,k in enumerate(lChild2Break)])+1 #I love list comprehension!

				for k in range(1, glob.iNumberMachines+1):
					if k not in dMachineConfig[sFamily1]:
						lAllowedMachines.append(k)

				#print("Child1 Current Machine, Allowed: ", iMyMachine, lAllowedMachines)

				lAllowedPositions = []
				iUpperBound = 0
				iLowerBound = 0

				for i,m in enumerate(lAllowedMachines):
					#set for First Machine
					if m == 1:
						iLowerBound = 0
						iUpperBound = lChild2Break[m-1]
					elif m==glob.iNumberMachines:
						iLowerBound = lChild2Break[m-2]
						iUpperBound = len(lChild2)

					else:
						iLowerBound = lChild2Break[m-2]
						iUpperBound = lChild2Break[m-1]


					lAllowedPositions.extend(list(range(iLowerBound, iUpperBound)))
					#set for last Machine

				#print("ALLOWED POSITIONS: ", lAllowedPositions)
				if iMyMachine not in lAllowedMachines:
					iForceAllocation = 1

				#print(glob.fMutationRate)

				if random.uniform(0.0, 1.0) < max(glob.fMutationRate, iForceAllocation):
					fRandXO2 = random.choice(lAllowedPositions)
					lChild2[iC1], lChild2[fRandXO2] = lChild2[fRandXO2], lChild2[iC1]



		#print("Parent breakers (Mum, Dad): ", lMotherBreak, lFatherBreak)
		#print("Children breakers (1,2): ", lChild1Break, lChild2Break)
		
		lChild1Check = []
		lChild1Check = udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList)
		#print("Genome AFTER: ", lChild1)
		#print("Illegals AFTER: ", udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList))

		lChild2Check = []
		#lChild2Check = udf_identifyIllegals(lChild2, lChild2Break, dMaterialFamily, dMachineConfig, dWcList)
		if glob.bDebug1 == True: print("Child 1 Check: ","child"+str(glob.iChildCounter)," :: " ,lChild1Check)

		#print(lChild1Check)

		#print("Child1_M: ", lChild1, " | ", lChild1Break)
		#print("Child2_M: ", lChild2, " | ", lChild2Break)
		##print("Child1_M_Ill: ", udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList))


		lChild1 = udf_listSortByBreak(lChild1, lChild1Break, 0)
		lChild2 = udf_listSortByBreak(lChild2, lChild2Break, 0)

		#print("Child1_S: ", lChild1, " | ", lChild1Break)
		#print("Child2_S: ", lChild2, " | ", lChild2Break)
		#print("Child1_S_Ill: ", udf_identifyIllegals(lChild1, lChild1Break, dMaterialFamily, dMachineConfig, dWcList))


		#print("MOTHER: ", mother)	

		glob.iChildCounter += 1
		lPopulation_offspring.append(lChild1)
		lPopulation_offspring_names.append("child"+str(glob.iChildCounter))
		dPopulation_offspring["child"+str(glob.iChildCounter)]={}
		dPopulation_offspring["child"+str(glob.iChildCounter)]["genome"]=lChild1
		dPopulation_offspring["child"+str(glob.iChildCounter)]["breaker"]=lChild1Break

		glob.iChildCounter += 1
		lPopulation_offspring.append(lChild2)
		lPopulation_offspring_names.append("child"+str(glob.iChildCounter))
		dPopulation_offspring["child"+str(glob.iChildCounter)]={}
		dPopulation_offspring["child"+str(glob.iChildCounter)]["genome"]=lChild2
		dPopulation_offspring["child"+str(glob.iChildCounter)]["breaker"]=lChild2Break



	return lPopulation_offspring, lPopulation_offspring_names, dPopulation_offspring


def udf_identifyIllegals(lMemberGenome, lBreaker, dMaterialFamily, dMachineConfig, dWcList):

	lIllegals = []
	iPreviousBreak = 0


	for k in range(0, glob.iNumberMachines):

		if k == glob.iNumberMachines-1:
			iNextBreak = len(lMemberGenome)

		else:
			iNextBreak = lBreaker[k]
				#iNextBreak = len(member["genome"])
		
		#print("Breakers: ", lBreaker, iPreviousBreak, iNextBreak)
		for i in lMemberGenome[iPreviousBreak:iNextBreak]:

			if i == 0:
				lIllegals.append(0)
				continue

			sMaterial1 = dWcList[i]['material']	# get family information
			sFamily1 = dMaterialFamily[sMaterial1]['family']

			if (k+1) in dMachineConfig[sFamily1]:
				lIllegals.append(1)
				
			else: lIllegals.append(0)

		iPreviousBreak = iNextBreak


	return lIllegals


def udf_listSortByBreak(lList, lBreaker, sSortBy):

	iPreviousBreak = 0
	iNextBreak = 0
	iCountSorts = 0

	lTransferList = lList
	lList = []

	for k in range(0, glob.iNumberMachines):
			
		if k == glob.iNumberMachines-1:
			iNextBreak = len(lTransferList)

		else:
			iNextBreak = lBreaker[k]

		lList_SortBys = []
		lList_Sorted = []
		
		lSortList = lTransferList[iPreviousBreak:iNextBreak]
		lList_Sorted = lTransferList[iPreviousBreak:iNextBreak]
		
		for i,j in enumerate(lSortList):
			if j == sSortBy:
				lList_SortBys.append(sSortBy)
				lList_Sorted.remove(sSortBy)

		lList.extend(lList_Sorted)
		lList.extend(lList_SortBys)
		
		iPreviousBreak = iNextBreak

	return lList


def udf_printMachinesCMD(lList, lBreaker, memberName):
	iPreviousBreak = 0
	iNextBreak = 0

	print(" #################################### PRINTING MACHINES:",memberName,"####################################")


	for k in range(0, glob.iNumberMachines):
		print("Machine #",str(k+1),":   ",end="")

		if k == glob.iNumberMachines-1:
			iNextBreak = len(lList)

		else: iNextBreak = lBreaker[k]

		for i,j  in enumerate(lList[iPreviousBreak:iNextBreak]):
			
			if j != 0:
				print(j, "->", end="")

		print(" |")

		iPreviousBreak = iNextBreak

	print("--------------------------------------")


def udf_printMachinesFamCMD(lList, lBreaker, memberName, dMaterialFamily, dWcList):
	iPreviousBreak = 0
	iNextBreak = 0

	print(" #################################### PRINTING MACHINES:",memberName,"####################################")


	for k in range(0, glob.iNumberMachines):
		print("Machine #",str(k+1),":   ",end="")

		if k == glob.iNumberMachines-1:
			iNextBreak = len(lList)

		else: iNextBreak = lBreaker[k]

		for i,j  in enumerate(lList[iPreviousBreak:iNextBreak]):
			
			if j != 0:
				sMaterial1 = dWcList[j]['material']	
				sFamily1 = dMaterialFamily[sMaterial1]['family']
				print(sFamily1, "->", end="")

		print(" |")

		iPreviousBreak = iNextBreak

	print("--------------------------------------")



######################################### RETIRED FUNCTIONS #########################################

# LEGACY_DONOT USE! calculate fitness of a population
def udf_calcFitness_LEGACY(lPopulation, dWcList, dMaterialFamily, dTimeMatrix, lMinFitness, lPopulation_names):
	'''
	LEGACY, DO NOT USE
	Calculates the Fitness of a given population based on Machine changeover and material processing times.
	'''
	# initialize all values
	sMaterial1 = ''
	sMaterial2 = ''
	sFamily1 = ''
	sFamily2 = ''
	sQuantity1 = 0
	sQuantity2 = 0
	sCycleTime1 = 0
	sCycleTime2 = 0
	sChangeovertime = 0
	fFitness = 0
	lFitness = []
	sMemberName = ''
	dMembers = {}

	for index0,member in enumerate(lPopulation):	# for every member in a population

		fFitness = 0 #reset fitness
		for index1, gene in enumerate(member):		# for every gene in a member
			
			if index1 < len(member)-1:				# if not the last material >> last material does not have changeover Times: BUT has processing times--> adj. needed
				sMaterial1 = dWcList[gene]['material']	# get family information
				sMaterial2 = dWcList[member[index1+1]]['material']
				sQuantity1 = dWcList[gene]['quantity']
				sQuantity2 = dWcList[member[index1+1]]['quantity']
				sFamily1 = dMaterialFamily[sMaterial1]['family']
				sFamily2 = dMaterialFamily[sMaterial2]['family']
				sCycleTime1 = dMaterialFamily[sMaterial1]['cycleTime']
				sCycleTime2 = dMaterialFamily[sMaterial2]['cycleTime']

				sChangeovertime = dTimeMatrix[str(sFamily1)+"-"+str(sFamily2)] # set changeover time by family

				fFitness += (sQuantity1*sCycleTime1)+(sQuantity2*sCycleTime2)+(sChangeovertime) # calculate overall fitness for every pair

		sMemberName = lPopulation_names[index0] # set the memberName based on previous input and iterations


		lFitness.append([sMemberName, fFitness]) # create fitness array

		# like the appendix, not really needed! 
		if fFitness <= lMinFitness[0]:
			lMinFitness[0] = fFitness
			lMinFitness[1] = sMemberName
			lMinFitness[2] = member

		# create member array with name, fitness and genome
		dMembers[sMemberName] = {}
		dMembers[sMemberName]['Fitness'] = fFitness
		dMembers[sMemberName]['genome'] = member

	return lFitness, dMembers, lMinFitness

# LEGACY_DONOT USE! calculate fitness of a population
def udf_calcFitness2_LEGACY(lPopulation, dPopulation, dWcList, dMaterialFamily, dTimeMatrix, lMinFitness, lPopulation_names, iNumberMachines, dMachineConfig):
	'''
	LEGACY, DO NOT USE
	Calculates the Fitness of a given population based on Machine changeover and material processing times.
	'''

	#initialize all values
	sMaterial1 = ''
	sMaterial2 = ''
	sFamily1 = ''
	sFamily2 = ''
	sQuantity1 = 0
	sQuantity2 = 0
	sCycleTime1 = 0
	sCycleTime2 = 0
	sChangeovertime = 0
	fFitness = 0
	lFitness = []
	sMemberName = ''
	dMembers = {}
	fMinFitness_run = 100000000
	iTotalRuns = 0.0
	iIllegalRuns = 0.0
	fIllegalPerc = 0.0

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< MEMBER LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for index0, (key,member) in enumerate(dPopulation.items()) :	# for every member in a population

		#print(index0, member["genome"], member["breaker"])

		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< MACHINE LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++

		iPreviousBreak = 0
		fFitness = 0
		iCountIndex1 = 0
		fFitnessBalance = []
		iMemberRuns = 0
		bisIllegal = False

		#if glob.bDebug1 == True: print("--------------------------------------------new member: ", key)

		#print(">>>> ", index0, member)
		iTotalRuns +=1
		iIllegalConfigMultiplier = 1

		for k in range(0, iNumberMachines):
			
			if k == iNumberMachines-1:
				iNextBreak = len(member["genome"])

			else:
				iNextBreak = member["breaker"][k]
				#iNextBreak = len(member["genome"])


			fFitnessM = 0 #reset fitness
			

			#print("-----------------------new sub-genome")
			
			#print("Sub-Member: ", key, member["genome"][iPreviousBreak:iNextBreak], member["breaker"])

			#+++++++++++++++++++++++++++++++++++++++++++++++++++++++< GENE LOOP >+++++++++++++++++++++++++++++++++++++++++++++++++++++++
			for index1, gene in enumerate(member["genome"][iPreviousBreak:iNextBreak]):		# for every gene in a member
				
					

				iMemberRuns +=1
				#print(member["genome"][iPreviousBreak:iNextBreak])
				if index1 < len(member["genome"][iPreviousBreak:iNextBreak])-1:				# if not the last material >> last material does not have changeover Times: BUT has processing times--> adj. needed
					sMaterial1 = dWcList[gene]['material']	# get family information
					sMaterial2 = dWcList[member["genome"][iCountIndex1+1]]['material']
					sQuantity1 = dWcList[gene]['quantity']
					sQuantity2 = dWcList[member["genome"][iCountIndex1+1]]['quantity']
					sFamily1 = dMaterialFamily[sMaterial1]['family']
					sFamily2 = dMaterialFamily[sMaterial2]['family']
					sCycleTime1 = dMaterialFamily[sMaterial1]['cycleTime']
					sCycleTime2 = dMaterialFamily[sMaterial2]['cycleTime']

					#print("ID:",iCountIndex1,":", sCycleTime1,sQuantity1,sCycleTime2,sQuantity2, sChangeovertime)

					sChangeovertime = dTimeMatrix[str(sFamily1)+"-"+str(sFamily2)] # set changeover time by family

					fFitnessM += (sQuantity1*sCycleTime1)+(sChangeovertime) #calculate overall fitness for every pair
					
					if (k+1) in dMachineConfig[sFamily1] or (k+1) in dMachineConfig[sFamily2]:
						iIllegalConfigMultiplier = 50
						bisIllegal = True
						#print("E1 Illegal config found", gene,member["genome"][iCountIndex1+1] , k+1, dMachineConfig[sFamily1], dMachineConfig[sFamily2])

				else: 
					#print("Adding last material processing w/o CO: ", gene)
					sMaterial1 = dWcList[gene]['material']	
					sQuantity1 = dWcList[gene]['quantity']
					sCycleTime1 = dMaterialFamily[sMaterial1]['cycleTime']
					sFamily1 = dMaterialFamily[sMaterial1]['family']
					fFitnessM += (sQuantity1*sCycleTime1)  

					if (k+1) in dMachineConfig[sFamily1]:
						iIllegalConfigMultiplier = 50
						bisIllegal = True
						#print("E2 Illegal config found", gene, k, dMachineConfig[sFamily1])
				



				iCountIndex1 += 1

				#print("FitnessM: ", fFitnessM)

			fFitness += fFitnessM
			fFitnessBalance.append(fFitnessM)

			iPreviousBreak = iNextBreak
			
			#print("Fitness All: ", fFitness)

		fFitness = (fFitness+stat.stdev(fFitnessBalance))*iIllegalConfigMultiplier

		if bisIllegal == True:
			iIllegalRuns +=1
			
			#if glob.bDebug1 == True: print("Illegal Stuff found")
			#if glob.bDebug1 == True: print("Fitness: ", fFitness, iIllegalConfigMultiplier)
			#if glob.bDebug1 == True: print("Others: ", member["genome"], member["breaker"])

		else: 
			if glob.bDebug1 == True: print("--------------------------------------------new member: ", key)
			if glob.bDebug1 == True: print("Fitness: ", fFitness, iIllegalConfigMultiplier)
			if glob.bDebug1 == True: print("Others: ", fFitnessBalance)

		#print(fFitnessBalance, stat.stdev(fFitnessBalance))
		# add stDev of distribution as penalty term
		
		#print("Fitness: ", fFitness)
		#if iIllegalConfigMultiplier >1:
		#	print(key,iIllegalConfigMultiplier, fFitness ,member["genome"], member["breaker"])

		sMemberName = key # set the memberName based on previous input and iterations


		lFitness.append([sMemberName, fFitness]) # create fitness array



		#set minimum fitness if it is lower than the previous fitness
		if fFitness < lMinFitness[0]:
			lMinFitness[0] = fFitness
			lMinFitness[1] = sMemberName
			lMinFitness[2] = member["genome"]
			lMinFitness[3] = member["breaker"]
			lMinFitness[4] = "Generation: "+str(glob.iGenerationCount-1)

		if fFitness < fMinFitness_run:
			fMinFitness_run = fFitness

		#create member array with name, fitness and genome
		dMembers[sMemberName] = {}
		dMembers[sMemberName]['fitness'] = fFitness
		dMembers[sMemberName]['genome'] = member["genome"]
		dMembers[sMemberName]["illegal"] = iIllegalConfigMultiplier
		#print("-----member runs: ",iMemberRuns)

	if glob.bDebug1 == True: print("runs: ", iTotalRuns, iIllegalRuns, iIllegalRuns/iTotalRuns)
	fIllegalPerc = iIllegalRuns/iTotalRuns
	return lFitness, dMembers, lMinFitness, fMinFitness_run, fIllegalPerc

# LEGACY simple mating with simple crossover, exchanging whole parts of strings. NOT USED, as this leads to duplication
def udf_simpleMating(lPopulation_new):
	lPopulation_offspring=[]

	for index, mother in enumerate(lPopulation_new[::2]):
		fRandXO1 = random.randint(0,len(mother)/2)
		fRandXO2 = random.randint(len(mother)/2,len(mother))

		iFather = (index)+(index)+1
		lChild1 = []
		lChild2 = []


		lChild1.extend(mother[0:fRandXO1])
		lChild1.extend(lPopulation_new[iFather][fRandXO1:fRandXO2])
		lChild1.extend(mother[fRandXO2:])

		lChild2.extend(lPopulation_new[iFather][0:fRandXO1])
		lChild2.extend(mother[fRandXO1:fRandXO2])
		lChild2.extend(lPopulation_new[iFather][fRandXO2:])
		print('----------------> ' +str(index))
		print('Mother: ', mother)
		print('Father: ', lPopulation_new[iFather])
		print("Child1: ", lChild1)
		print("Child2: ", lChild2)

		lPopulation_offspring.append(lChild1)
		lPopulation_offspring.append(lChild2)

