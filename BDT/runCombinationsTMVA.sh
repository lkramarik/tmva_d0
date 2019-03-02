#!/usr/bin/env bash
ptmin=1
ptmax=2

folderPt=pt_${ptmin}_${ptmax}
mkdir -p ${folderPt}
cd ${folderPt}

for nTrees in 200 300 400
do
    for treeDepth in 2 3 4
	do
	    folderTMVA=n${nTrees}_d${treeDepth}
        mkdir -p ${folderTMVA}
        cd ${folderTMVA}
        cp ../../TMVAClassification.C ./
        cp ../../TMVAClassificationApplication.C ./
        cp ../../TMVAClassificationApplicationMixed.C ./
        cp ../../tmvaCuts.h ./
        #root -l TMVAClassification.C++\($ptmin,$ptmax\)
        root -l -q TMVAClassification.C++\($ptmin,$ptmax,$nTrees,$treeDepth\)
        root -l -q TMVAClassificationApplication.C++\(\"./../../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)
        #root -l -q TMVAClassificationApplicationMixed.C++\(\"out_local_mix.root\",${ptmin},${ptmax}\)
        cp -r ../../analyse ./
        cd analyse
        rm -r signals*
        root -l -q project_bdt.C++\(${ptmin},${ptmax},${nTrees},${treeDepth}\)
        #root -l project_bdt_mixed.cpp++\(${ptmin},${ptmax}\)
        cp significance_pt* ../../
        cd ../..
    done
done






