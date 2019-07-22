# --------CHANGE OVERS

######################################### 0 Intro ######################################### 
# Naming Convention: first letter of variable indicates the type
# a = array
# b = binary / boolean
# c = code, for .py files only
# d = dictionary
# f = float
# g = graph
# i = integer
# l = list
# lim = limit
# s = string
# file = file generated or modified through code

import pyodbc as db
import random
import pandas as pd
import numpy as np
import c02_geneticAlgorithmFunctions as gak
import c03_globalVariables as glob
import sys
import matplotlib.pyplot as plt
import math
import datetime


fMinFitness = 100000000
iGenerationCount = 0
lMinFitness = [10000000, 'START', [],[],""]
lMinFitness_history = [100000000]
lFitness_history=[]
lIllegal_history=[]
######################################### 1 DATA IMPORT ######################################### 

# LOAD orders from SQL
# -> 2 look up structure:
# 		1) materials family
# 		2) changeOver matrix
#		3) Product Machine restrictions
#		4) WC List

# TO DO:


### 1.1 Get Material Family Data
# for now based on excel; check c13_ImportFromSQL.py for SQL import cide


dFamilyCO = {}
dMaterialFamily = {}
dWcList = {}
dMachineConfig = {}
dMaterialCO ={}

glob.lFamilyAtlas_0 = []
glob.lMaterialAtlas_0 = []

#import from Excel
dfWCImport = pd.read_excel("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\03_co_setup_alternative.xlsx", sheet_name="order")
dfFamilies = pd.read_excel("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\03_co_setup_alternative.xlsx", sheet_name="families")
dfFamilyCO = pd.read_excel("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\03_co_setup_alternative.xlsx", sheet_name="familyCO")
dfMachineConfig = pd.read_excel("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\03_co_setup_alternative.xlsx", sheet_name="notOnMachine")
dfMaterialCO = pd.read_excel("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\03_co_setup_alternative.xlsx", sheet_name="materialCO")

#fill WC List
for index, row in dfWCImport.iterrows():
	if index >= glob.iBreakImport:
		break

	dWcList[row.orderNumber]={}
	dWcList[row.orderNumber]['material'] = row.materialCode
	dWcList[row.orderNumber]['quantity'] = row.quantity
	dWcList[row.orderNumber]['priority'] = row.priority



#Create TimeMatrix dictionary from Query Results
for index, row in dfFamilyCO.iterrows():
	dFamilyCO[row.relID]= row['time']
	glob.lFamilyAtlas_0.append(row["familyAtlas"])



#Create materialFamily dictionary from Query Results
for index, row in dfFamilies.iterrows():
	dMaterialFamily[row.Material] = {}
	dMaterialFamily[row.Material]['family'] = row.materialFamily
	dMaterialFamily[row.Material]['cycleTime'] = row.cycleTime
#Create MachineConfig >> ILLEGAL MACHINE CONFIG, machines the family is not allowed on
for index, row in dfMachineConfig.iterrows():
	dMachineConfig[row.family] = [int(x) for x in str(row.notOnMachine).split(",")]

#Create Material changeover time mapping
for index, row in dfMaterialCO.iterrows():
	dMaterialCO[row.materialRel] = row["timeCO"]
	glob.lMaterialAtlas_0.append(row["materialAtlas"])


#open file to track usage history
filePopulationHistory = open("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\90_populationHistory.txt", "w", encoding="utf-8")
fileFitnessHistory_runs = open("C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\91_fitnessRuns.txt", "a", encoding="utf-8")

######################################### 2 GA SETUP ######################################### 
# TO DO
# > use a more intelligent filling for initial population 
### 2.1 Iterate over WC list and populate arrays
# initialize population randomly from list
lGenome = []
lPopulation = []
dPopulation = {}
lPopulation_names =[]
glob.lGenome_0 = []

# create original genome with all orders contained
for order in dWcList.keys():
	glob.lGenome_0.append(order)

# create list of 0s to fill fill the genome to a length of machines x genes to represent all possible machines in one genome
lEmptyAppend = [i*0 for i in range(0, (glob.iNumberMachines-1)*len(glob.lGenome_0))]
lGenome = glob.lGenome_0+lEmptyAppend


# from the filled Genome, create n = limPopulationSize initial parents
for i in range(0,glob.limPopulationSize):
	lPopulation.append([])
	lGenesAvailable = lGenome.copy()

	lPopulation_names.append("member"+str(i))
	lBreakGenome = []
	lCOlookup = []

	# fill member randomly based on available genes - reduce available genes
	for gene in range(0, len(lGenome)):
		sOrder= random.choice(lGenesAvailable)
		lPopulation[i].append(sOrder)
		lGenesAvailable.remove(sOrder)

	# craete breaking points / machine seperators based on the lenght of a single machine (length of the origin genome)
	for j in range(1, glob.iNumberMachines):
		iRandBreaker = (j*len(glob.lGenome_0))
		lBreakGenome.append(iRandBreaker)

	# legacy sort from when the breaker was generated randomly
	lBreakGenome.sort()
	gak.udf_listSortByBreak(lPopulation[i], lBreakGenome, 0)


	# populate the Population dictionary
	dPopulation["member"+str(i)] = {}
	dPopulation["member"+str(i)]["genome"] = lPopulation[i]
	dPopulation["member"+str(i)]["breaker"] = lBreakGenome	

# write the first population to the history file
filePopulationHistory.write("#"+str(iGenerationCount)+".1------------------------ Original Population ------------------------"+"\n")
for i,w in enumerate(lPopulation):
	filePopulationHistory.write(lPopulation_names[i]+": "+str(w)+"\n")



######################################### 3 GA Algorithm ######################################### 
# ! Arrays ending on "_names" are parallel arrays to track member names
# iterate until break point reached (see below)
iBreakLoop = glob.iBreakGeneration


while iGenerationCount < iBreakLoop:

	fIllegalPerc = 0.0

	iGenerationCount += 1
	print("--------------------------------- GENERATION: "+str(iGenerationCount)+"---------------------------------")


	# execute function to calculate fitness of population
	# determine randomly if a cataclysm should occur; cataclsym = "kills off" the population and fills it with newly created one
	if random.uniform(0.0, 1.0) < glob.iCataclysmicProb and glob.bCataclysm == True:

		print("<<<<<<<<<<<<<<<<<<< CATACLYSM TIME <<<<<<<<<<<<<<<<<<<")

		dPopulation = gak.udf_cataclysm(dPopulation, glob.lGenome_0)
		# Add runs to the overall counter after cataclysm
		glob.iCataclysmicProb = glob.iCataclysmicProb/2
		iBreakLoop += glob.iBreakGeneration

	# calculte fitness for each member in the population
	
	lFitness, dMembers, lMinFitness, fMinFitness_run, fIllegalPerc = gak.udf_calcFitness3(dPopulation, dWcList, dMaterialFamily, dFamilyCO, dMaterialCO, lMinFitness, dMachineConfig, iGenerationCount)
	lFitness_history.append(fMinFitness_run)
	lIllegal_history.append(fIllegalPerc)

	# if the fitness is lower then the previous fintness lever, update the minimum fitness
	if lMinFitness[0] <= fMinFitness:
		fMinFitness = lMinFitness[0]

	# append calculated fitness for new lowest level
	lMinFitness_history.append(fMinFitness)

	# create table and calculate selection probabilities
	lFitness_sorted = gak.udf_sortByFitness(lFitness)


	# initialize population arrays
	lPopulation_new = []
	lPopulation_new_names = []
	dPopulation_new ={}

	# select parents randomly to form new population
	lPopulation_new, lPopulation_new_names, dPopulation_new = gak.udf_selectParentsFromPool(dMembers, lFitness_sorted, dPopulation)


	# Mating time - execute mating functions and initialize offspring arrays
	lPopulation_offspring = []
	lPopulation_offspring_names = []
	dPopulation_offspring ={}

	# lPopulation_offspring, glob.iChildCounter, lPopulation_offspring_names, dPopulation_offspring = gak.udf_matingPMX(lPopulation_new, glob.iChildCounter, lPopulation_new_names, dPopulation_new, dMembers, glob.fMutationRate)
	lPopulation_offspring, lPopulation_offspring_names, dPopulation_offspring = gak.udf_cloneMutate(lPopulation_new, lPopulation_new_names, dPopulation_new, dMembers, dMaterialFamily, dMachineConfig, dWcList, lGenome)


	# Mutating Time - execute swap-mutate function
	lPopulation_offspringMutated = []
	lPopulation_offspringMutated, dPopulation_offspring = gak.udf_mutateSwap(glob.fMutationRate,lPopulation_offspring, dPopulation_offspring)
	

	# overwrite previous population with new selected and offspring
	lPopulation = lPopulation_new + lPopulation_offspringMutated
	lPopulation_names = lPopulation_new_names + lPopulation_offspring_names

	# recreate the population array with the selected parents from previous iteration
	dPopulation={}
	for i,member in dPopulation_new.items():
		# avoid double entries, which are technically possible due to selection method
		if member["member"] not in dPopulation:
			dPopulation[member["member"]]={}
			dPopulation[member["member"]]["genome"]=member["genome"]
			dPopulation[member["member"]]["breaker"]=member["breaker"]

	# deconstruct newly created parent array and the mutated offspring array
	dPopulation = {**dPopulation, **dPopulation_offspring}

	# calculate starting point for trailing average
	iAvgStart = len(lMinFitness_history)-glob.iPastAverage
	if iAvgStart < 5:
		iAvgStart = 0

	# break the while loop if no lower fitness could be found for "iAvgStart" number of generations
	if sum(lMinFitness_history[(iAvgStart):(len(lMinFitness_history))]) < fMinFitness:
		break



# close file
filePopulationHistory.close()

# terminal end prints
print("===============================================================================================")
print("RESULT: ", lMinFitness[0])
print("MEMBER: ", lMinFitness[1])
print( lMinFitness[4])

print(np.corrcoef(lFitness_history, lIllegal_history)[1])

# print machines in termial
gak.udf_printMachinesCMD(lMinFitness[2], lMinFitness[3], lMinFitness[1])
print("__________________________________________")

# print machines with familes not materials
gak.udf_printMachinesFamCMD(lMinFitness[2], lMinFitness[3], lMinFitness[1], dMaterialFamily, dWcList)


######################################### 4 Graphing it #########################################
# set min and max for the y axes
y1Min = math.floor(min(lFitness_history)/1000)*1000
y1Max = math.ceil(min(lFitness_history)/1000)*2000
y2Min = math.floor(min(lIllegal_history))-0.1
y2Max = math.ceil(min(lIllegal_history))+0.1

# set parameters for saving the plot
dateTime = datetime.datetime.now()
iMilliseconds = int(round(dateTime.timestamp() * 1000))

sFileNamePlot = str(iMilliseconds)+"__RESULT_"+str(math.floor(lMinFitness[0]))+"__orders_"+str(len(glob.lGenome_0))+"--machines_"+str(glob.iNumberMachines)+"--Runs_"+str(glob.iBreakGeneration)+"--popSize_"+str(glob.limPopulationSize)+"--mut_"+str(glob.fMutationRate)+"--King_"+str(glob.bKingPrevails)+"--fAlloc_"+str(glob.iForceAllocation_G)+"--CAT_"+str(glob.bCataclysm)+"_"+str(glob.iCataclysmicProb)+"_"+str(glob.iDeletionProb)
sPlotPath = "C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\99_Output\\"+sFileNamePlot+".png"

# create subplot
gFitness, ax1 = plt.subplots()

# set options
color = "tab:blue"
ax1.set_ylabel("Fitness")
ax1.set_xlabel("Runs")
ax1.set_ylim(y1Min, y1Max)
ax1.plot(lFitness_history, color=color)

# create twin axis and set 
ax2 = plt.twinx()
color ="tab:green"
ax2.set_ylabel("Illegal Percentage")
ax2.set_ylim(y2Min, y2Max)
ax2.plot(lIllegal_history, color=color, linestyle="--")

gFitness.tight_layout()

#save and plot
plt.savefig(sPlotPath)
plt.show()

fileFitnessHistory_runs.write(str(lMinFitness[0])+"\n")

############################################## THE END ##############################################