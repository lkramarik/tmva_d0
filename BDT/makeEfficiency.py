import os, sys
import shutil
import datetime
import time
import getpass
import array
from multiprocessing import Pool

nCores = 6

#############################################
def func(input, ptmin, ptmax, weight):
	print(input)
	command=f"root -l -q -b \'efficiencyScripts/makeEfficiencyOneBin.cpp+(\"{str(input)}\",{ptmin},{ptmax},\"{str(weight)}\")\' "
	print(command)
	os.system(command)
	return True

#############################################
def makeSetup(anaSetup):
	pool = Pool(processes=int(nCores))

	print(inputFileNames)
	for j in range(len(inputFileNames)):
		for i in range(nPtBins):
			pool.apply_async(func, args=(inputFileNames[j], ptmin[i], ptmax[i], weight[j]))
			time.sleep(1)

	pool.close()
	pool.join()

	for j in range(len(inputFileNames)):
		command1=f"root -l -q -b \'efficiencyScripts/plotTogetherEfficiencyCompare.cpp+(\"{str(inputFileNames[j])}\",{nPtBins},{ptmin[0]},{ptmax[len(ptmax)-1]},{int(anaSetup)})\'"
		print(command1)
		os.system(command1)

		directoryWithFilename=f"finalAnalysis/{directory}/{str(inputFileNames[j])}"
		if os.path.exists(directoryWithFilename):
			shutil.rmtree(directoryWithFilename)
		os.makedirs(directoryWithFilename)
		command2=f"cp finalAnalysis/{str(inputFileNames[j])}/final_result_SIM.root {directoryWithFilename}/."
		os.system(command2)

#############################################
if __name__ == '__main__':
	workDir = os.getcwd()

	# inputFileNames=['ntp_FS_hijing_goodEvent_hotspot_1401']
	# inputFileNames=["ntp.FS.1cm.all.newE.1509.goodPid",
	# 				"ntp.FS.data.hft2tof2.global.3009.goodPid",
	# 				"ntp.FS.1cm.primaries.newE.1509.goodPid",
	# 				"ntp.FS.hft2.global.1multEdge.0810.goodPid",
	# 				"ntp.FS.hft2.global.4multEdge.0810.goodPid",
	# 				"ntp_FS_hijing_goodEvent_hotspot_1401.goodPid",
	# 				"ntp_FS_hijing_goodevent_1401.goodPid"]
	#
	# inputFileNames=["ntp_fullEvent_full_production.vtx.3M.1308.goodEvent",
	# 				"ntp_fullEvent_full_production.vtx.3M.1308.goodEvent.hotspot"]

	# inputFileNames=["ntp_FS_hijing_goodevent_1401.goodPid",
	# 				"ntp_fullEvent_full_production.vtx.3M.1308.goodEvent.goodPid"]
	#
	inputFileNames=["ntp_FS_hijing_goodEvent_hotspot_1401.goodPid","ntp_fullEvent_full_production.vtx.3M.1308.goodEvent.hotspot"]
	#
	weight=['',
			'']

	# weight=['weight',
	# 		'weight',
	# 		'weight',
	# 		'weight',
	# 		'weight',
	# 		'',
	# 		'']

	legend=['FastSim','HIJING']

	#narrow bins
	directory = "narrowBins"
	ptmin=[]
	ptmax=[]

	nPtBins = 16
	for i in range(nPtBins):
		ptmin.append(1.+i*0.25)
		ptmax.append(1.+(i+1.)*0.25)

	print(ptmin)
	print(ptmax)

	makeSetup(0)

	#wide bins
	directory = "wideBins"
	ptmin=[1,2,3]
	ptmax=[2,3,5]

	nPtBins = 3

	print(ptmin)
	print(ptmax)

	makeSetup(1)



# Ratio
	command3=f"root -l -q -b \'efficiencyScripts/plotTogetherEfficiencyCompareImages.cpp++(\"{str(inputFileNames[0])}\",\"{str(inputFileNames[1])}\",\"{str(directory)}\",\"{str(legend[0])}\",\"{str(legend[1])}\")\'"
	os.system(command3)



