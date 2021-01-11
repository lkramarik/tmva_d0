#include "TH1.h"
#include "TRatioPlot.h"
#include "/home/lukas/work/tmva_d0/BDT/BDTCutEstimate.cpp"
#include "/home/lukas/work/tmva_d0/BDT/analysisSetup.h"

//_________________________________________________________________________________________________
void makeEfficiencyOneBin(TString input="ntp_fullEvent_full_production.1M.vtx",
                       Float_t ptMin=1.,
                       Float_t ptMax=2,
                       TString weight="") {
    gROOT->ProcessLine(".L analyse/FitD0Peak.cpp+");

    setCurrentFolder(input); //defined in BDTCutEstimate

    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    TString prefix = "out_local_SIM_";

    precutsTMVA = Form("k_pt>%1.2f && pi1_pt>%1.2f && "
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && "
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

    TString folder = Form("finalAnalysis/%s", input.Data());
    gSystem->Exec(Form("mkdir -p %s", folder.Data()));
//    gSystem->Exec(Form("mkdir -p finalAnalysis/%s/effStudy/%s", folderDate.Data(), input.Data()));
//                gSystem->Exec(Form("cp /home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/%s%s finalAnalysis/%s/effStudy/.", ptMin[i], ptMax[i], (int) nTrees[i], (int) depth[i], prefix.Data(), input[j].Data(), folderDate.Data()));

    Float_t ptminFilename=2, ptmaxFilename=5;
    Int_t analysisBin=-1;
    float_t  nTrees, depth, bdtResponse;
    for (int iPt = 0; iPt < analysisSetup::nPtBins; ++iPt) {
        if (ptMin>=analysisSetup::ptBinMin[iPt] && ptMax<=analysisSetup::ptBinMax[iPt]){
            ptminFilename=analysisSetup::ptBinMin[iPt];
            ptmaxFilename=analysisSetup::ptBinMax[iPt];
            bdtResponse=analysisSetup::bdtCut[iPt];
            depth=analysisSetup::maxDepth[iPt];
            nTrees=analysisSetup::nTrees[iPt];
            analysisBin=iPt;
        }
    }
    if (analysisBin==-1){
        cout<<"selected pT doesnot correspond to any analysis pT bin, ending."<<endl;
        return;
    }

    cout<<"projecting from file corresponding to pT bin: "<<ptminFilename<<" "<<ptmaxFilename;

    TString inputFileAdressName = Form("/home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/%s%s.root", ptminFilename, ptmaxFilename, (int) nTrees, (int) depth, prefix.Data(), input.Data());

    mBDTCut = Form("BDTresponse>=%.3f", analysisSetup::bdtCut[analysisBin]);
    mCuts.push_back("mcEtas>0");

//    mCuts.push_back("mcEtas>0 && hft>0 && tpc>0 && pid>0");

    project_bdt_oneCut_SIM(ptMin, ptMax, nTrees, depth, bdtResponse, inputFileAdressName, "ntp_signal", weight);
//                project_bdt_oneCut_SIM(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i], folder + prefix + input[j]);

    mCuts.clear();
    mCuts.shrink_to_fit();
}