# --------CHANGE OVERS

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


		

######################################### 0 GLOBAL  VARIABLES ######################################### 


print(glob.iCataclysmicProb)

# fMutationRate = 0.2
# glob.limPopulationSize = 10	# has to be multiple of 2
# glob.iChildCounter = glob.limPopulationSize
# glob.iGenerationCount = 0
# glob.iBreakGeneration = 5
# glob.iPastAverage = 5
# glob.iNumberMachines = 4
# glob.iDeletionProb = 0.25
# glob.iCataclysmicProb = 0.50
# glob.iBreakImport = 10


fMinFitness = 100000000
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


#Create TimeMatrix dictionary from Query Results
for index, row in dfFamilyCO.iterrows():
	dFamilyCO[row.relID]= row['time']
	glob.lFamilyAtlas_0.append(row["familyAtlas"])



#Create materialFamily dictionary from Query Results
for index, row in dfFamilies.iterrows():
	dMaterialFamily[row.Material] = {}
	dMaterialFamily[row.Material]['family'] = row.materialFamily
	dMaterialFamily[row.Material]['cycleTime'] = row.cycleTime

for index, row in dfMachineConfig.iterrows():
	dMachineConfig[row.family] = [int(x) for x in str(row.notOnMachine).split(",")]


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


for order in dWcList.keys():
	glob.lGenome_0.append(order)

lEmptyAppend = [i*0 for i in range(0, (glob.iNumberMachines-1)*len(glob.lGenome_0))]

lGenome = glob.lGenome_0+lEmptyAppend



for i in range(0,glob.limPopulationSize):
	lPopulation.append([])
	lGenesAvailable = lGenome.copy()
	#print(i)
	lPopulation_names.append("member"+str(i))
	lBreakGenome = []
	lCOlookup = []

	#fill member randomly based on available genes - reduce available genes
	for gene in range(0, len(lGenome)):
		sOrder= random.choice(lGenesAvailable)
		lPopulation[i].append(sOrder)
		lGenesAvailable.remove(sOrder)


	for j in range(1, glob.iNumberMachines):
		
		iRandBreaker = (j*len(glob.lGenome_0))
			#print("iBreaker in List, restarting")

		lBreakGenome.append(iRandBreaker)

	lBreakGenome.sort()

	#print(lCOlookup)
	#print(gak.udf_listSortByBreak(lCOlookup, lBreakGenome, 0))



	dPopulation["member"+str(i)] = {}
	dPopulation["member"+str(i)]["genome"] = gak.udf_listSortByBreak(lPopulation[i], lBreakGenome, 0)
	dPopulation["member"+str(i)]["breaker"] = lBreakGenome	




#for i,j in dPopulation.items():
#	gak.udf_printMachinesCMD(j["genome"], j["breaker"], i)


filePopulationHistory.write("#"+str(glob.iGenerationCount)+".1------------------------ Original Population ------------------------"+"\n")
for i,w in enumerate(lPopulation):
	filePopulationHistory.write(lPopulation_names[i]+": "+str(w)+"\n")


#print(lPopulation)
#print("---------------------")
#print(dPopulation)

### 2.2 Estimate fitness of each member



######################################### 3 GA Algorithm ######################################### 
#! Arrays ending on "_names" are parallel arrays to track member names
#iterate until break point reached (see below)
while True:

	glob.iGenerationCount += 1
	print("--------------------------------- GENERATION: "+str(glob.iGenerationCount)+"---------------------------------")


	#execute function to calculate fitness of population

	########## THE CATACLYSM ##########
	if random.uniform(0.0, 1.0) < glob.iCataclysmicProb and glob.bCataclysm == True:

		print("<<<<<<<<<<<<<<<<<<< CATACLYSM TIME <<<<<<<<<<<<<<<<<<<")

		#for i,j in dPopulation.items():
		#	print(i,j)
	
		dPopulation = gak.udf_cataclysm(dPopulation, glob.lGenome_0)
		glob.iBreakGeneration +=100

	#print("-------------- New post apocalyptic society --------------")
	fIllegalPerc=0.0

	lFitness, dMembers, lMinFitness, fMinFitness_run, fIllegalPerc = gak.udf_calcFitness3(lPopulation, dPopulation, dWcList, dMaterialFamily, dFamilyCO, dMaterialCO, lMinFitness, lPopulation_names, glob.iNumberMachines, dMachineConfig)

	lFitness_history.append(fMinFitness_run)
	lIllegal_history.append(fIllegalPerc)


	#print(fMinFitness)

	#print(lFitness)
	#print(dMembers)
	#print(lMinFitness)


	#if the fitness is lower then the previous fintness lever, update the minimum fitness
	if lMinFitness[0] <= fMinFitness:
		fMinFitness = lMinFitness[0]

	#append calculated fitness for new lowest level
	lMinFitness_history.append(fMinFitness)

	#create table and calculate selection probabilities
	lFitness_sorted = gak.udf_sortByFitness(lFitness)


	# for i,j in dMembers.items():
	# 	print(i, j["fitness"], j["illegal"])



	#initialize population arrays
	lPopulation_new = []
	lPopulation_new_names = []
	dPopulation_new ={}

	#select parents randomly to form new population
	lPopulation_new, lPopulation_new_names, dPopulation_new = gak.udf_selectParentsFromPool(dMembers, glob.limPopulationSize, lFitness_sorted, dPopulation)


	#Mating time - execute mating functions and initialize offspring arrays
	lPopulation_offspring = []
	lPopulation_offspring_names = []
	dPopulation_offspring ={}

	#lPopulation_offspring, glob.iChildCounter, lPopulation_offspring_names, dPopulation_offspring = gak.udf_matingPMX(lPopulation_new, glob.iChildCounter, lPopulation_new_names, dPopulation_new, dMembers, glob.fMutationRate)
	lPopulation_offspring, lPopulation_offspring_names, dPopulation_offspring = gak.udf_cloneMutate(lPopulation_new, lPopulation_new_names, dPopulation_new, dMembers, dMaterialFamily, dMachineConfig, dWcList, lGenome)


	#Mutating Time - execute swap-mutate function
	lPopulation_offspringMutated = []
	lPopulation_offspringMutated, dPopulation_offspring = gak.udf_mutateSwap(glob.fMutationRate,lPopulation_offspring, dPopulation_offspring)
	

	#overwrite previous population with new selected and offspring
	lPopulation = lPopulation_new + lPopulation_offspringMutated
	lPopulation_names = lPopulation_new_names + lPopulation_offspring_names
	#dPopulation = {**dPopulation_new, **dPopulation_offspring}

	dPopulation={}
	for i,member in dPopulation_new.items():


		if member["member"] not in dPopulation:
			dPopulation[member["member"]]={}
			dPopulation[member["member"]]["genome"]=member["genome"]
			dPopulation[member["member"]]["breaker"]=member["breaker"]






	dPopulation = {**dPopulation, **dPopulation_offspring}
		#for i, j in dPopulation.items():
	#	print(i, j)
	# if lMinFitness[1] in dPopulation:

	# 	print("#E7 PREF CHILD: ", lMinFitness[1], dPopulation[lMinFitness[1]]["genome"])
	


	#break the while loop if generation limit is reached
	if glob.iGenerationCount == glob.iBreakGeneration:
		break


	#calculate starting point for trailing average
	iAvgStart = len(lMinFitness_history)-glob.iPastAverage
	if iAvgStart < 5:
		iAvgStart = 0


	#break the while loop if no lower fitness could be found for "iAvgStart" number of generations
	if sum(lMinFitness_history[(iAvgStart):(len(lMinFitness_history))]) < fMinFitness:
		break

#	for i,j in dPopulation.items():
#		print(i)
			
#close file
filePopulationHistory.close()
print("===============================================================================================")
print("RESULT: ", lMinFitness[0])
print("MEMBER: ", lMinFitness[1])
print( lMinFitness[4])

#print(lFitness_history)
#print(lFitness_history)
#print(lIllegal_history)

print(np.corrcoef(lFitness_history, lIllegal_history)[1])

gak.udf_printMachinesCMD(lMinFitness[2], lMinFitness[3], lMinFitness[1])
print("__________________________________________")
gak.udf_printMachinesFamCMD(lMinFitness[2], lMinFitness[3], lMinFitness[1], dMaterialFamily, dWcList)
############################### GRAPH IT

y1Min = math.floor(min(lFitness_history)/1000)*1000
y1Max = math.ceil(min(lFitness_history)/1000)*1100

y2Min = math.floor(min(lIllegal_history))-0.1
y2Max = math.ceil(min(lIllegal_history))+0.1


dt = datetime.datetime.now()
iMilliseconds = int(round(dt.timestamp() * 1000))

sFileNamePl0t = str(iMilliseconds)+"__RESULT_"+str(math.floor(lMinFitness[0]))+"__orders_"+str(len(glob.lGenome_0))+"--machines_"+str(glob.iNumberMachines)+"--Runs_"+str(glob.iBreakGeneration)+"--popSize_"+str(glob.limPopulationSize)+"--mut_"+str(glob.fMutationRate)+"--King_"+str(glob.bKingPrevails)+"--fAlloc_"+str(glob.iForceAllocation_G)+"--CAT_"+str(glob.bCataclysm)+"_"+str(glob.iCataclysmicProb)+"_"+str(glob.iDeletionProb)
sPlotPath = "C:\\Users\\u374441\\Desktop\\desktopWorkfiles\\201905 Changeovers\\99_Output\\"+sFileNamePl0t+".png"

fig, ax1 = plt.subplots()

color = "tab:blue"
ax1.set_ylabel("Fitness")
ax1.set_xlabel("Runs")
ax1.set_ylim(y1Min, y1Max)
ax1.plot(lFitness_history, color=color)

ax2=plt.twinx()
color ="tab:green"
ax2.set_ylabel("Illegal Percentage")
ax2.set_ylim(y2Min, y2Max)
ax2.plot(lIllegal_history, color=color, linestyle="--")

fig.tight_layout()


plt.savefig(sPlotPath)
plt.show()


fileFitnessHistory_runs.write(str(lMinFitness[0])+"\n")

############################################## THE END ##############################################