#!/usr/bin/env bash
ptmin=3
ptmax=5

folderPt=pt_${ptmin}_${ptmax}
mkdir -p ${folderPt}
cd ${folderPt}

for nTrees in  300 100 400
#for nTrees in 800
do
  for treeDepth in 3
	do
	    folderTMVA=n${nTrees}_d${treeDepth}
         mkdir -p ${folderTMVA}
        cd ${folderTMVA}
        cp ../../TMVAClassification.C ./
        cp ../../TMVAClassificationApplication.C ./
        cp ../../TMVAClassificationApplicationSIM.C ./
        cp ../../TMVAClassificationApplicationMixed.C ./
        cp ../../tmvaCuts.h ./
        root -l -q -b TMVAClassification.C++\($ptmin,$ptmax,$nTrees,$treeDepth\)
        root -l -q -b TMVAClassificationApplication.C++\(\"./../../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)
        root -l -q -b TMVAClassificationApplicationSIM.C++\(\"out_local_SIM.root\",${ptmin},${ptmax}\)

        cp ./../../files_to_run.list ./
        #root -l -q TMVAClassificationApplicationMixed.C++\(\"out_local_mix.root\",${ptmin},${ptmax}\)
        cp -r ../../analyse ./
        cd analyse
        rm -r signals*
        root -l -q -b project_bdt.C++\(${ptmin},${ptmax},${nTrees},${treeDepth}\)
        #root -l project_bdt_mixed.cpp++\(${ptmin},${ptmax}\)
        cp significance_pt* ../../
        cd ../../
    done
done






