import os, sys
import shutil
import datetime
import time
import getpass
import array
from multiprocessing import Pool

nCores = 3

def func(command,folder):
    print(command)
    os.chdir(folder)
    os.system(str(command))
    return True

#hlavni
if __name__ == '__main__':
    # simFolder='/home/lukas/work/D0-fullEvent/analyse/ntp/'
    simFolder='/home/lukas/work/tmva_d0/sim/'


    # simFile='ntp_full_D0.toyMc.data.0408.dca1cm.root'
    # simFile='ntp_FS_hijing_goodEvent_hotspot_3101.root'

    # simFile='ntp_FS_hijing_goodevent_hotspot_0302.root'
    # simFile='ntp_FS_data_global_HS_dca1_0302.root'
    # simFile='ntp_FS_hijing_goodevent_0302.root'
    # simFile='ntp_FS_data_primaries_HS_dca1_0302.root'
    # simFile='ntp_FS_hijing_goodevent_0302.root'
    # simFile='ntp_FS_hijing_goodevent_hotspot_0302.root'

    # simFile='ntp_FS_data_primaries_HS_dca1_0302.root'

    simFile=['ntp_FS_hijing_goodevent_0302.root', 'ntp_FS_hijing_goodevent_hotspot_0302.root']
    # simFile=['ntp_fullEvent_full_production.3M.vtx.0109.vtxMaxMult.goodevent.hotspot.root', 'ntp_fullEvent_full_production.3M.vtx.0109.vtxMaxMult.goodevent.root']

    # simFile='ntp_fullEvent_full_production.3M.vtx.0109.vtxMaxMult.goodevent.root'
    # simFile='ntp_fullEvent_full_production.3M.vtx.0109.vtxMaxMult.goodevent.hotspot.root'


    # TString inputPlot[] = {"ntp_FS_data_primaries_HS_dca1_0302.root",
    #                        "ntp_FS_data_global_HS_dca1_0302.root",
    #                        "ntp_FS_hijing_goodevent_hotspot_0302.root",
    #                        "ntp_FS_hijing_goodevent_0302.root"};

    # simFile='ntp_fullEvent_full_production.vtx.3M.1308.goodEvent.hotspot.root'
    # simFile='ntp_FS_hijing_nonPrimary_DCA1cm_goodEvent.root'
    # simFile='ntp_FS_hijing_nonPrimary_DCA1cm.root'
    # simFile='ntp_FS_hijing_nonPrimary_DCA1cm.root'
    # simFile='ntp_full_D0.toyMc.2407.tpceffDca1.newHijingInputs.root'
    # simFile='ntp_full_D0.toyMc.dca1cm.HIJING.0608.root'
    # simFile='ntp_full_D0.toyMc.data.0408.dca1cm.root'
    # simFile='/home/lukas/work/tmva_d0/sim/ntp_full_D0.toyMc.2407.tpceffDca1.newHijingInputs.root'

    folders=['/pt_1_2/n150_d3/','/pt_2_3/n150_d3/','/pt_3_5/n400_d3/']
    ptmin=[1,2,3]
    ptmax=[2,3,5]
    workDir = os.getcwd()

    pool = Pool(processes=int(nCores))


    for jFile in range(2):
        for j in range(3):
            print(str(workDir)+str(folders[j]))
            os.chdir(str(workDir)+str(folders[j]))
            print(os.getcwd())
            os.system('cp  ../../TMVAClassificationApplicationSIM.C ./')
            os.system('cp ../../tmvaCuts.h ./')
            outFile="out_local_SIM_%s.root" % (str(os.path.splitext(simFile[jFile])[0]))
            print(outFile)
            command="root -l -q -b \'%s/TMVAClassificationApplicationSIM.C++(\"%s\",\"%s\",%f,%f)\' " % (str(workDir)+str(folders[j]),str(simFolder)+str(simFile[jFile]),str(workDir)+str(folders[j])+str(outFile),ptmin[j],ptmax[j])
            pool.apply_async(func, args=(command,str(workDir)+str(folders[j])))
            time.sleep(1)

    pool.close()
    pool.join()
