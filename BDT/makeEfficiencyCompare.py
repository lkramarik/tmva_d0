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
	directory = "/home/lukas/work/tmva_d0/BDT/finalAnalysis/narrowBins/effStudy"
	# if os.path.exists(directory):
	# 	shutil.rmtree(directory)
	# os.makedirs(directory)

    #Vytvori pool n valekn pro beh kodu a nastaveni param
	pool = Pool(processes=int(nCores))

	ptmin=[]
	ptmax=[]
	nTrees=[]
	depth=[]
	bdtResponse=[]

	folder='narrowBins'
	nPtBins = 16
	for i in range(16):
		ptmin.append(1.+i*0.25)
		ptmax.append(1.+(i+1.)*0.25)

	for i in range(16):
		if ptmin[i]>=1 and ptmax[i]<=2:
			nTrees.append(100)
			depth.append(3)
			bdtResponse.append(0.7552)

		if ptmin[i]>=2 and ptmax[i]<=3:
			nTrees.append(150)
			depth.append(3)
			bdtResponse.append(0.64516)

		if ptmin[i]>=3 and ptmax[i]<=5:
			nTrees.append(400)
			depth.append(3)
			bdtResponse.append(0.53154)

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

	# inputFileNames=['ntp_full_D0.toyMc.2407.2.tpceffDca1.newHijingInputs.weight', 'ntp_fullEvent_full_production.3M.vtx.2107.weight']
	# inputFileNames=['ntp_full_D0.toyMc.dca1cm.HIJING.0608', 'ntp_fullEvent_full_production.3M.vtx.2107.weight']

	# inputFileNames=['ntp_full_D0.toyMc.hijing.nonPrimaryInputs.dca1cm.1208.weight',
	# 				'ntp_fullEvent_full_production.3M.vtx.2107.weight']

	inputFileNames=['ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ',
					'ntp_fullEvent_full_production.vtx.3M.1308.AllD0']
	weight=['','']
	legend=['FastSim','HIJING']



	# inputFileNames=['ntp_full_D0.toyMc.hijing.nonPrimaryInputs.dca1cm.1208.2', 'ntp_FS_hijing_nonPrimary_DCA1cm_goodEvent',
	# 				'ntp_FS_hijing_nonPrimary_DCA1cm', 'ntp_full_D0.toyMc.dca1cm.HIJING.0608.weight', 'ntp_fullEvent_full_production.vtx.3M.1507']
	# weight=['','','','']


	# inputFileNames=['ntp_full_D0.toyMc.2307.tpceffdca1', 'ntp_fullEvent_full_production.3M.vtx.2107']
	# inputFileNames=['ntp_full_D0.toyMc.2407.tpceffDca1.newHijingInputs.weight', 'ntp_fullEvent_full_production.3M.vtx.2107']

	# inputFileNames=['ntp_full_D0.toyMc.2307.tpceffdca1.weight', 'ntp_fullEvent_full_production.3M.vtx.2107']
	# inputFileNames=['ntp_full_D0.toyMc.1507.2.weight', 'ntp_fullEvent_full_production.vtx.3M.1507']

	# inputFileNames=['ntp_fullEvent_full_production.vtx.3M.geant.primaryTag.1307','ntp_full_D0.toyMc.hijing.1407.weight.TPCdata.weight']
	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx','ntp_full_D0.toyMc.hijing.1407.weight.TPCdata.weight']
	# inputFileNames=['ntp_fullEvent_full_production.1M.vtx','ntp_full_D0.toyMc.0303.0707']
	# weight=['','']
	# weight=['weight','']
	print(inputFileNames)

	# for j in range(4):
	for j in range(2):
		for i in range(16):
			pool.apply_async(func, args=(inputFileNames[j], ptmin[i], ptmax[i], nTrees[i], depth[i], bdtResponse[i],weight[j]))
			time.sleep(1)

	pool.close()
	pool.join()

	# for j in range(4):
	for j in range(2):
		command1="root -l -q -b \'/home/lukas/work/tmva_d0/BDT/efficiencyScripts/plotTogetherEfficiencyCompare.cpp+(\"%s\",16,%f,%f)\' " % (str(inputFileNames[j]),ptmin[0],ptmax[len(ptmax)-1])
		print(command1)
		os.system(command1)
		directory="finalAnalysis/narrowBins/effStudy/%s" % (str(inputFileNames[j]))
		if os.path.exists(directory):
			shutil.rmtree(directory)
		os.makedirs(directory)
		command2="cp finalAnalysis/%s/final_result_SIM.root finalAnalysis/narrowBins/effStudy/%s/." % (str(inputFileNames[j]),str(inputFileNames[j]))
		os.system(command2)

	command3="root -l -q -b \'/home/lukas/work/tmva_d0/BDT/efficiencyScripts/plotTogetherEfficiencyCompareImages.cpp++(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\")\' " % (str(inputFileNames[0]), str(inputFileNames[1]),str(folder),str(legend[0]),str(legend[1]))
	os.system(command3)



