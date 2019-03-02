#!/bin/bash
ptmin=3
ptmax=5
treeDepth=3
nTrees=350

folderPt=pt_${ptmin}_${ptmax}
folderTMVA=n${nTrees}_d${treeDepth}
mkdir -p $folderPt
#mkdir -p $folderPt/$folderTMVA
#cd $folderPt/$folderTMVA
cd $folderPt/ ##mixed

#cp -i ../../TMVAClassification.C ./
#cp -i ../../TMVAClassificationApplication.C ./
cp -i ../TMVAClassificationApplicationMixed.C ./
#cp -i ../../tmvaCuts.h ./

#cp -i ../TMVAGui.C ./
#cp -i ../tmvaglob.C ./

#root -l TMVAClassification.C++\($ptmin,$ptmax\)
#root -l -q TMVAClassification.C++\($ptmin,$ptmax,$nTrees,$treeDepth\)
#root -l -q TMVAClassificationApplication.C++\(\"./../../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)
root -l -q TMVAClassificationApplicationMixed.C++\(\"out_local_mix.root\",${ptmin},${ptmax}\)
#
cp -r ../analyse ./
#cp -r ../../analyse ./
cd analyse
rm -r signals*
#root -l project_bdt.C++\(${ptmin},${ptmax},${nTrees},${treeDepth}\)
root -l project_bdt_mixed.cpp++\(${ptmin},${ptmax}\)
