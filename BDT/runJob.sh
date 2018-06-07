#!/bin/bash
#!/bin/csh
ptmin=3
ptmax=5
analyzer="lkramarik"
productionId=`date +%F_%H-%M`
echo ${productionId}
mkdir -p pt_${ptmin}_${ptmax}/workDir
mkdir -p workDir/${productionId}
cp /gpfs01/star/pwg/lkramarik/tmva_pm/files_to_run.list   pt_${ptmin}_${ptmax}/workDir/${productionId}
list="files_to_run.list"

cp -r /gpfs01/star/pwg/lkramarik/tmva_d0/BDT/TMVAClassificationApplication.C  pt_${ptmin}_${ptmax}/workDir/${productionId}
cp -r /gpfs01/star/pwg/lkramarik/tmva_d0/BDT/pt_${ptmin}_${ptmax}/weights  pt_${ptmin}_${ptmax}/workDir/${productionId}
cp -r /gpfs01/star/pwg/lkramarik/tmva_pm/submit   pt_${ptmin}_${ptmax}/workDir/${productionId}


path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

baseFolder="/gpfs01/star/pwg/lkramarik/tmva_pm/BDT/pt_"${ptmin}"_"${ptmax}"/workDir/"
echo ${baseFolder}

cd pt_${ptmin}_${ptmax}/workDir/${productionId}
echo pt_${ptmin}_${ptmax} > ptrange
mkdir -p production
mkdir -p report
mkdir -p csh
mkdir -p list
mkdir -p log
mkdir -p err
mkdir -p jobs
mkdir -p jobs/log
mkdir -p jobs/err

echo ${baseFolder}${productionId}

#star-submit-template -template submit/submit.xml -entities basePath=${baseFolder}${productionId},prodId=${productionId}