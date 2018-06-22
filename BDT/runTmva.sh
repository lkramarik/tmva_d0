#!/bin/bash
ptmin=2
ptmax=3
#sim=""

folderPt=pt_${ptmin}_${ptmax}
mkdir -p $folderPt
cd $folderPt

cp -i ../TMVAClassification.C ./
cp -i ../TMVAClassificationApplication.C ./
cp -i ../tmvaCuts.h ./

cp -i ../TMVAGui.C ./
cp -i ../tmvaglob.C ./

#root -l TMVAClassification.C++\($ptmin,$ptmax\)
root -l -q TMVAClassification.C++\($ptmin,$ptmax\)
root -l -q TMVAClassificationApplication.C++\(\"./../files_to_run.list\",\"out_local.root\",${ptmin},${ptmax}\)
#
cp -r ../analyse ./
cd analyse
rm -r signals*
root -l project_bdt.C++\(${ptmin},${ptmax}\)
