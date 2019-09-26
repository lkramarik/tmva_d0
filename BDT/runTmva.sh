#!/bin/bash
ptmin=2
ptmax=3
treeDepth=2
nTrees=300

folderPt=pt_${ptmin}_${ptmax}
folderTMVA=n${nTrees}_d${treeDepth}
mkdir -p $folderPt
mkdir -p $folderPt/$folderTMVA
cd $folderPt/$folderTMVA
#cd $folderPt/ ##mixed

cp -i ../../TMVAClassification.C ./
cp -i ../../TMVAClassificationApplication.C ./
cp -i ../../TMVAClassificationApplicationSIM.C ./
#cp -i ../TMVAClassificationApplicationMixed.C ./
cp -i ../../tmvaCuts.h ./

root -l -q TMVAClassification.C++\($ptmin,$ptmax,$nTrees,$treeDepth\)
root -l -q TMVAClassificationApplication.C++\(\"./../../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)
root -l -q TMVAClassificationApplicationSIM.C++\(\"out_local_SIM.root\",${ptmin},${ptmax}\)

#root -l -q TMVAClassificationApplicationMixed.C++\(\"out_local_mix.root\",${ptmin},${ptmax}\)
#
#cp -r ../analyse ./
cp -r ../../analyse ./
cd analyse
rm -r signals*
root -l -b -q project_bdt.C++\(${ptmin},${ptmax},${nTrees},${treeDepth}\)
#root -l project_bdt_mixed.cpp++\(${ptmin},${ptmax}\)
