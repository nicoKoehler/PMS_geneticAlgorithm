#GLOBAL VARIABLES

global fMutationRate 
global limPopulationSize 
global iChildCounter 
global iBreakGeneration 
global iPastAverage 
global iNumberMachines
global iDeletionProb 
global iCataclysmicProb 
global iBreakImport
#global bMutateChild
global bDebug1 
global bDebug2
global bDebug3
global iForceAllocation_G
#global lGenome_0
#global lMaterialAtlas_0
#global lFamilyAtlas_0
global bKingPrevails

fMutationRate = 0.12	# best @ 1% 
limPopulationSize = 250	# has to be multiple of 2
iBreakGeneration = 10
iBreakImport = 10000
iForceAllocation_G = 0.5

iPastAverage = 5
iNumberMachines = 4
iDeletionProb = 0.5
iCataclysmicProb = 0.02
iChildCounter = limPopulationSize

bKingPrevails = True
bCataclysm = True
bCorrectChild = True
bDebug1 = False
bDebug2 = False
bDebug3 = False

