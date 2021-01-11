#include "/home/lukas/work/tmva_d0/BDT/BDTCutEstimate.cpp"
#include "/home/lukas/work/tmva_d0/BDT/analysisSetup.h"

void plotTogetherEfficiencyCompare(TString input,
        const Int_t nBins,
        Double_t ptmin,
        Double_t ptmax){

    gROOT->ProcessLine(".L analyse/FitD0Peak.cpp+");
    setCurrentFolder(input);

    Double_t ptMin[nBins];
    Double_t ptMax[nBins];
    Double_t nTrees[nBins];
    Double_t depth[nBins];
    Double_t bdtResponse[nBins];

    Double_t ptstep = (ptmax-ptmin)/(float)nBins;
    for (int j = 0; j < nBins; ++j) {
        ptMin[j] = 1+j*ptstep;
        ptMax[j] = 1+(j+1)*ptstep;
    }

    cout<<analysisSetup::nPtBins<<endl;
    Float_t ptminFilename=2, ptmaxFilename=5;
    Int_t analysisBin=-1;
    for (int i = 0; i < nBins; ++i) {
        for (int iPtAnalysis = 0; iPtAnalysis < analysisSetup::nPtBins; ++iPtAnalysis) {
            if (ptMin[i]>=analysisSetup::ptBinMin[iPtAnalysis] && ptMax[i]<=analysisSetup::ptBinMax[iPtAnalysis]) {
                bdtResponse[i] = analysisSetup::bdtCut[iPtAnalysis];
                depth[i] = analysisSetup::maxDepth[iPtAnalysis];
                nTrees[i] = analysisSetup::nTrees[iPtAnalysis];
            }
        }
    }

    for (int i = 0; i < nBins; ++i) {
        cout<<ptMin[i]<<" "<<ptMax[i]<<" "<<bdtResponse[i]<<" "<<nTrees[i]<<" "<<depth[i]<<endl;
    }

    plotTogetherSIM(nBins, ptMin, ptMax, nTrees, depth, bdtResponse);
}