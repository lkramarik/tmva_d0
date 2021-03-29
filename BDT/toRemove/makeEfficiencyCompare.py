import os, sys
import shutil
import datetime
import time
import getpass
import array
from multiprocessing import Pool

nCores = 1
		
def func(input, ptmin, ptmax, weight):
	print(input)
	command=f"root -l -q -b \'efficiencyScripts/makeEfficiencyOneBin.cpp+(\"{str(input)}\",{ptmin},{ptmax},\"{str(weight)}\")\' "
	print(command)
	os.system(command)
	return True

#hlavni
if __name__ == '__main__':
    #Smaze a znovu vytvori slozku kde jsou ulozeny vysledky
	workDir = os.getcwd()
	directory = "finalAnalysis/narrowBins/"

    #Vytvori pool n valekn pro beh kodu a nastaveni param
	pool = Pool(processes=int(nCores))

	ptmin=[]
	ptmax=[]

	nPtBins = 2
	for i in range(nPtBins):
		ptmin.append(1.+i*0.25)
		ptmax.append(1.+(i+1.)*0.25)

	print(ptmin)
	print(ptmax)




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

	# inputFileNames=['ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ',
	# 				'ntp_fullEvent_full_production.vtx.3M.1308.AllD0']
	#
	inputFileNames=['ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ']
	weight=['','']
	legend=['FastSim','HIJING']

	print(inputFileNames)

	# for j in range(len(inputFileNames)):
	# 	for i in range(nPtBins):
	# 		pool.apply_async(func, args=(inputFileNames[j], ptmin[i], ptmax[i], weight[j]))
	# 		time.sleep(1)

	pool.close()
	pool.join()

	for j in range(len(inputFileNames)):
		command1=f"root -l -q -b \'efficiencyScripts/plotTogetherEfficiencyCompare.cpp+(\"{str(inputFileNames[j])}\",{nPtBins},{ptmin[0]},{ptmax[len(ptmax)-1]})\'"
		print(command1)
		os.system(command1)

		directoryWithFilename=f"{directory}/{str(inputFileNames[j])}"
		if os.path.exists(directoryWithFilename):
			shutil.rmtree(directoryWithFilename)
		os.makedirs(directoryWithFilename)
		command2=f"cp finalAnalysis/{str(inputFileNames[j])}/final_result_SIM.root {directoryWithFilename}/."
		os.system(command2)

	# command3="root -l -q -b \'/home/lukas/work/tmva_d0/BDT/efficiencyScripts/plotTogetherEfficiencyCompareImages.cpp++(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\")\' " % (str(inputFileNames[0]), str(inputFileNames[1]),str(folder),str(legend[0]),str(legend[1]))
	# os.system(command3)



