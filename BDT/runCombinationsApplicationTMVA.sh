#!/usr/bin/env bash
ptmin=3
ptmax=4

folderPt=pt_${ptmin}_${ptmax}
cd ${folderPt}

#for nTrees in 100 200 300 400
for nTrees in 500
#for nTrees in 300 400 500 600 700
do
    for treeDepth in 3
	do
	    folderTMVA=n${nTrees}_d${treeDepth}
        cd ${folderTMVA}
        cp ../../TMVAClassificationApplication.C ./
        cp ../../TMVAClassificationApplicationSIM.C ./
        cp ../../tmvaCuts.h ./
        root -l -q -b TMVAClassificationApplication.C++\(\"./../../files_to_run_app.list\",\"out_local.root\",${ptmin},${ptmax}\)
        cp ./../../files_to_run_app.list ./
        root -l -q TMVAClassificationApplicationSIM.C++\(\"out_local_SIM.root\",${ptmin},${ptmax}\)

        #root -l -q TMVAClassificationApplicationMixed.C++\(\"out_local_mix.root\",${ptmin},${ptmax}\)
        cp -r ../../analyse ./
        cd analyse
        rm -r signals*
        root -l -q -b project_bdt.C++\(${ptmin},${ptmax},${nTrees},${treeDepth}\)
        #root -l project_bdt_mixed.cpp++\(${ptmin},${ptmax}\)
        cp significance_pt* ../../
        cd ../..
    done
done






