#include "/home/lukas/work/tmva_d0/BDT/BDTCutEstimate.cpp"
#include "/home/lukas/work/tmva_d0/BDT/analysisSetup.h"

void fit(TGraphErrors* grToFit);
void correctYield();

void makeFinalCalculation() {
    gROOT->ProcessLine(".L /home/lukas/work/tmva_d0/BDT/analyse/FitD0Peak.cpp++");

//    setCurrentFolder("finalBinsTest");
    setCurrentFolder("finalBins");
//    simulationFileName="out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root";
//    simulationFileName="out_local_SIM_ntp.FS.1cm.all.newE.1509.root";
//    simulationFileName="out_local_SIM_ntp_full_D0.toyMc.data.0408.dca1cm.goodPid.root";
//    simulationFileName="out_local_SIM_ntp.FS.1cm.primaries.newE.1509.root";
    simulationFileName="out_local_SIM_ntp_FS_data_global_HS_dca1_0302.goodPid_dca1.root";
    dataFileName="out_local.root";

    gSystem->Exec(Form("mkdir -p finalAnalysis/%s", folderDate.Data()));
    gSystem->Exec(Form("mkdir finalAnalysis/%s/img", folderDate.Data()));
    gSystem->Exec(Form("mkdir finalAnalysis/%s/img/root", folderDate.Data()));

    precutsTMVA = Form(
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<1 &&"
                       "pi1_dca>%f && pi1_dca<1 &&"
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

//    for (int i = 0; i < 1; ++i) {
    for (int i = 0; i < analysisSetup::nPtBins; ++i) {
        mCuts.push_back(Form("D_pt>=%.3f && D_pt<%.3f", analysisSetup::ptBinMin[i], analysisSetup::ptBinMax[i]));
        mBDTCut=Form("BDTresponse>=%.3f", analysisSetup::bdtCut[i]);

        doRawYield(analysisSetup::ptBinMin[i], analysisSetup::ptBinMax[i], analysisSetup::nTrees[i], analysisSetup::maxDepth[i], analysisSetup::bdtCut[i]);
        mCuts.clear();
        mCuts.shrink_to_fit();
    }

    plotTogether(analysisSetup::nPtBins, analysisSetup::ptBinMin, analysisSetup::ptBinMax, analysisSetup::nTrees, analysisSetup::maxDepth, analysisSetup::bdtCut);

    drawSignificanceData(analysisSetup::nPtBins, analysisSetup::ptBinMin, analysisSetup::ptBinMax, analysisSetup::nTrees, analysisSetup::maxDepth, analysisSetup::bdtCut);
    cout<<"Work is done. Everything is in folder "<<folderDate<<endl;

}


void makeFinalCalculationNarrow() {
    gROOT->ProcessLine(".L /home/lukas/work/tmva_d0/BDT/analyse/FitD0Peak.cpp++");

    setCurrentFolder("finalNarrowBins");
    simulationFileName="out_local_SIM_ntp_full_D0.toyMc.data.0408.dca1cm.goodPid.root";
    dataFileName="out_local_ntp.D0.2608.root";

    gSystem->Exec(Form("mkdir -p finalAnalysis/%s", folderDate.Data()));
    gSystem->Exec(Form("mkdir finalAnalysis/%s/img", folderDate.Data()));
    gSystem->Exec(Form("mkdir finalAnalysis/%s/img/root", folderDate.Data()));

    Double_t ptMin[16];
    Double_t ptMax[16];
    Double_t nTrees[16];
    Double_t depth[16];
    Double_t bdtResponse[16];

    for (int j = 0; j < 16; ++j) {
        ptMin[j] = 1+j*0.25;
        ptMax[j] = 1+(j+1)*0.25;
    }
    const int nBins=sizeof(ptMin)/ sizeof(Double_t);

    TString inputSim[nBins];

    for (int k = 0; k < nBins; ++k) {
        if (ptMin[k]>=1 && ptMax[k]<=2){
            nTrees[k]=100;
            depth[k]=3;
            bdtResponse[k]=0.7552;
        }

        if (ptMin[k]>=2 && ptMax[k]<=3){
            nTrees[k]=150;
            depth[k]=3;
            bdtResponse[k]=0.64516;
        }

        if (ptMin[k]>=3 && ptMax[k]<=5){
            nTrees[k]=400;
            depth[k]=3;
            bdtResponse[k]=0.53154;
        }
    }

    precutsTMVA = Form(
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && "
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

    for (int i = 0; i < nBins; ++i) {
        Float_t ptmin=2, ptmax=5;
        if (ptMin[i]>=1 && ptMax[i]<=2){
            ptmin=1;
            ptmax=2;
        }

        if (ptMin[i]>=2 && ptMax[i]<=3){
            ptmin=2;
            ptmax=3;
        }

        if (ptMin[i]>=3 && ptMax[i]<=5){
            ptmin=3;
            ptmax=5;
        }
        TString inputFileAdressName = Form("/home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/%s", ptmin, ptmax, (int) nTrees[i], (int) depth[i], simulationFileName.Data());
//        mCuts.push_back("refMult>10");

        mBDTCut=Form("BDTresponse>=%.3f", bdtResponse[i]);

        project_bdt_oneCut_SIM(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i], inputFileAdressName, "ntp_signal", "");

        mCuts.clear();
        mCuts.shrink_to_fit();
    }

    plotTogetherSIM(nBins, ptMin, ptMax, nTrees, depth, bdtResponse);
    cout<<"Work is done. Everything is in folder "<<folderDate<<endl;
}




