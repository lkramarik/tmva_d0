#!/bin/bash
ptmin=2
ptmax=3

folderPt=pt_${ptmin}_$ptmax
mkdir -p $folderPt
cd $folderPt

cp -i ../TMVAClassification.C ./
cp -i ../TMVAClassificationApplication.C ./
cp -i ../tmvaCuts.h ./

cp -n ../TMVAGui.C ./
cp -n ../tmvaglob.C ./

root -l -q TMVAClassification.C++\($ptmin,$ptmax\)
root -l TMVAClassificationApplication.C++\(\"nic\", \"out_local.root\",$ptmin,$ptmax\)