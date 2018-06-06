#!/bin/bash
ptmin=3
ptmax=5

folderPt=pt_${ptmin}_$ptmax
mkdir -p $folderPt
cd $folderPt

cp -i ../TMVAClassification.C ./
cp -i ../TMVAClassificationApplication.C ./
cp -i ../tmvaCuts.h ./

cp -n ../TMVAGui.C ./
cp -n ../tmvaglob.C ./

#root -l TMVAClassification.C++\($ptmin,$ptmax\)
#root -l -q TMVAClassification.C++\($ptmin,$ptmax\)
root -l -q TMVAClassificationApplication.C++\(\"./../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)

cp -r ../analyse ./
cd analyse
root -l project_bdt.C++\(${ptmin},${ptmax}\)
