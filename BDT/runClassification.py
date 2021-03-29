import os, sys
import shutil
import datetime
import time
import getpass
import array
from multiprocessing import Pool

nCores = 6

def train(ptmin,ptmax,nTrees,maxDepth):
    os.chdir(workDir)
    folder=f"/pt_{ptmin}_{ptmax}/n{nTrees}_d{maxDepth}/"
    print(str(folder))
    os.system(f"mkdir -p {str(workDir)+str(folder)}")
    # os.system(f"rm -r {str(workDir)+str(folder)}*")

    os.system(f"cp -r ./analyse {str(workDir)+str(folder)}")
    os.system(f"cp TMVAClassification.C {str(workDir)+str(folder)}")
    os.system(f"cp TMVAClassificationApplication.C {str(workDir)+str(folder)}")
    os.system(f"cp tmvaCuts.h {str(workDir)+str(folder)}")

    # Classification
    # command=f"root -l -q -b \'TMVAClassification.C+({ptmin},{ptmax},{nTrees},{maxDepth})\' "
    # print(command)
    # func(command,str(workDir)+str(folder))
    #
    # #ROC curve
    # filename=f"TMVA_bdt_d0_pt_{ptmin:.1f}_{ptmax:.1f}.root"
    # print(filename)
    # command=f"root -l -q -b \'makeTxtFiles.C++(\"{filename}\")\'"
    # func(command,str(workDir)+str(folder))

    # Application
    command=f"root -l -q -b \'TMVAClassificationApplication.C+(\"./../../files_to_run.list\",\"out_local.root\",{ptmin},{ptmax})\' "
    func(command,str(workDir)+str(folder))


    os.system(f"rm -r {str(workDir)+str(folder)}/analyse/signals*")
    command=f"root -l -q -b \'project_bdt.C++({ptmin},{ptmax},{nTrees},{maxDepth})\' "
    analyse="analyse"
    func(command,str(workDir)+str(folder)+str(analyse))
    os.system(f"cp significance_pt* ../../")


def func(command,folder):
    print(command)
    os.chdir(folder)
    os.system(str(command))
    return True

#hlavni
if __name__ == '__main__':
    simFolder='/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/'
    simFile='ntp.D0.2608.root'
    fileAdress='%s%s' % (str(simFolder),str(simFile))
    outFileName="out_"+simFile
    print(outFileName)
    text_file = open("files_to_run.list", "w")
    n = text_file.write(str(fileAdress))
    text_file.close()

    #current ana setup:
    ptmin=[1, 2, 3]
    ptmax=[2, 3, 5]
    nTrees=[150, 150, 400]
    maxDepth=[3, 3, 3]

    # ptmin=[1]
    # ptmax=[2]
    # nTrees=[150]
    # maxDepth=[3]

    workDir = os.getcwd()

    pool = Pool(processes=int(nCores))

    # for j in range(len(ptmin)):
    #     for nT in nTrees:
    #         for depth in maxDepth:
    #             pool.apply_async(train, args=(ptmin[j], ptmax[j], nT, depth))

    for j in range(len(ptmin)):
        pool.apply_async(train, args=(ptmin[j], ptmax[j], nTrees[j], maxDepth[j]))

    pool.close()
    pool.join()
