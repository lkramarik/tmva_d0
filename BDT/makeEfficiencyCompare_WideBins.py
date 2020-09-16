import os, sys
import shutil
import datetime
import time
import getpass
import array
from multiprocessing import Pool

nCores = 6
		
#Funkce t
def func(input, ptmin, ptmax, ntrees, depth, bdtresponse,weight):
	print(input)
	command="root -l -q -b \'efficiencyScripts/efficiencyCompare.cpp+(\"%s\",%f,%f,%f,%f,%f,\"%s\")\' " % (str(input),ptmin,ptmax,ntrees,depth,bdtresponse,weight)
	# command="root -l -q -b \'efficiencyCompare.cpp+(\"%s\",%f,%f,%f,%f,%f,\"%s\")\' " % (str(input),ptmin,ptmax,ntrees,depth,bdtresponse,weight)
	print(command)
	os.system(command)
	return True

#hlavni
if __name__ == '__main__':
    #Smaze a znovu vytvori slozku kde jsou ulozeny vysledky
	workDir = os.getcwd()
	print workDir
	folder='wideBins'

	directory = "/home/lukas/work/tmva_d0/BDT/finalAnalysis/wideBins/effStudy"
	# if os.path.exists(directory):
	# 	shutil.rmtree(directory)
	# os.makedirs(directory)

    #Vytvori pool n valekn pro beh kodu a nastaveni param
	pool = Pool(processes=int(nCores))

	ptmin=[1,2,3]
	ptmax=[2,3,5]
	nTrees=[100,150,400]
	depth=[3,3,3]
	bdtResponse=[0.7552,0.64516,0.53154]

	print(ptmin)
	print(ptmax)
	print(nTrees)
	print(depth)
	print(bdtResponse)

	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx','ntpTMVA_full_D0.toyMc.0303']
	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx','ntpTMVA_full_D0.toyMc.0303.2303']
	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx.weight','ntpTMVA_full_D0.toyMc.0303.2303.weight']
	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx','ntp_full_D0.toyMc.hijing.1407.weight']
	# inputFileNames=['ntp_fullEvent_full_production.vtx.3M.geant.primaryTag.1307','ntp_full_D0.toyMc.1507.2.weight']
	# inputFileNames=['ntp_full_D0.toyMc.1507.2.weight', 'ntp_fullEvent_full_production.vtx.3M.geant.primaryTag.1307']
	# inputFileNames=['ntp_full_D0.toyMc.1507.2.weight', 'ntp_fullEvent_full_production.3M.vtx.2107.weight']

	# inputFileNames=['ntp_full_D0.toyMc.2407.2.tpceffDca1.newHijingInputs', 'ntp_fullEvent_full_production.3M.vtx.2107']
	inputFileNames=['ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ', 'ntp_fullEvent_full_production.vtx.3M.1308.AllD0']
	# inputFileNames=['ntp_FS_hijing_nonPrimary_DCA1cm_goodEvent', 'ntp_fullEvent_full_production.vtx.3M.1308.goodEvts']
	legend=['FastSim','HIJING']
	weight=['','']
	# weight=['weight','']
	print(inputFileNames)

	for j in range(2):
		for i in range(3):
			pool.apply_async(func, args=(inputFileNames[j], ptmin[i], ptmax[i], nTrees[i], depth[i], bdtResponse[i],weight[j]))
			time.sleep(1)

	pool.close()
	pool.join()

	for j in range(2):
		command1="root -l -q -b \'/home/lukas/work/tmva_d0/BDT/efficiencyScripts/plotTogetherEfficiencyCompare_wideBins.cpp+(\"%s\")\' " % (str(inputFileNames[j]))
		print(command1)
		os.system(command1)
		directory="finalAnalysis/%s/effStudy/%s" % (str(folder),str(inputFileNames[j]))
		if os.path.exists(directory):
			shutil.rmtree(directory)
		os.makedirs(directory)
		command2="cp finalAnalysis/%s/final_result_SIM.root finalAnalysis/%s/effStudy/%s/." % (str(inputFileNames[j]),str(folder),str(inputFileNames[j]))
		os.system(command2)

	command3="root -l -q -b \'/home/lukas/work/tmva_d0/BDT/efficiencyScripts/plotTogetherEfficiencyCompareImages.cpp++(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\")\' " % (str(inputFileNames[0]), str(inputFileNames[1]),str(folder),str(legend[0]),str(legend[1]))
	os.system(command3)



