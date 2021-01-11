#include "efficiencyCalculation.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TSystem.h"
#include "TROOT.h"
#include "StAnaCuts.h"
#include "TLorentzVector.h"


#include "ntpFastSimBDT.h"
void makeEff(TString, TString, TString, TString);

void testMake() {
    gROOT->ProcessLine(".L efficiencyCalculation.cpp++");
    gROOT->LoadMacro("ntpFastSimBDT.C++");

//    TFile *data = new TFile("/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root","read");
//    TNtuple *ntp = (TNtuple*) data->Get("ntp_signal");
//    ntp->MakeClass("ntpFastSimBDT");

//    TString coeffHistoName[] = {"_coeff0.95", "", "_coeff1.05"};
//    TString fileNameSuffix[] = {"nHitsFit15", "nHitsFit17", "nHitsFit20"};

    TString coeffHistoName[] = {""};
    TString fileNameSuffix[] = {"goodPid"};

    TString folderTpcName[] = {"/home/lukas/work/embedding_dAu/analyse/tpc_eff_dca1/"};

//    TString folderTpcName[] = {"/home/lukas/work/embedding_dAu/analyse/tpc_eff_nHitsFit15/", "/home/lukas/work/embedding_dAu/analyse/tpc_eff_nHitsFit17/", "/home/lukas/work/embedding_dAu/analyse/tpc_eff_nHitsFit20/"};

//    TString ntpFileName = "/home/lukas/work/tmva_d0/sim/ntpTMVA_full_D0.toyMc.0303.root";

    TString ntpFileName[] = {
//                             "/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.1003.105perc.root",
                             "/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntp_full_D0.toyMc.data.0408.dca1cm.root",
                             "/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntp_full_D0.toyMc.data.0408.dca1cm.root",
                             "/home/lukas/work/tmva_d0/BDT/pt_3_5/n400_d3/out_local_SIM_ntp_full_D0.toyMc.data.0408.dca1cm.root"};

//                                "/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_3_5/n400_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root"};

//    TString ntpFileName[] = {"/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.105perc.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.105perc.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_3_5/n400_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.105perc.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.95perc.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.95perc.root",
//                             "/home/lukas/work/tmva_d0/BDT/pt_3_5/n400_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.95perc.root",
//                             };

    for (int k = 0; k < 3; ++k) {
        cout<<"----------"<<ntpFileName[k]<<"----------"<<endl;
        for (int j = 0; j < 1; ++j) { //nHitsFit
            cout<<"----------"<<fileNameSuffix[j]<<"----------"<<endl;
            cout<<"----------"<<folderTpcName[j]<<"----------"<<endl;

            for (int i = 0; i < 1; ++i) { //coefficients
                cout<<"----------"<<coeffHistoName[i]<<"----------"<<endl;

                makeEff(coeffHistoName[i], folderTpcName[j], ntpFileName[k], fileNameSuffix[j]);
            }
        }
    }

}

//__________________________________________________________________________________________________________
void makeEff(TString coeffHistoName, TString folderTpcName, TString ntpFileName, TString fileNameSuffix){
    efficiencyCalculation eff;
    eff.setHistoSuffix(coeffHistoName);
    eff.setFolderTpc(folderTpcName);
    eff.setTpcGraphs();
    eff.setTofMatch(1);

    //PIDs
    int bbcMin=0;
    int bbcMax=950;
    int nTof=0;
    float nsigma=3;
    float tofInvBeta=0.03;
    float ptTrackCut=0.;

    TString name;
//    TString cutComb=Form("bbc%i_%i_nHft%i_nsigma%.1f_tof%.2f_pt%.1f/",bbcMin, bbcMax, nTof, nsigma, tofInvBeta, ptTrackCut);
//    name="/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/tofPidEff_"+cutComb+"rootFiles/";

    eff.setFolderTofPid("pid_files/bbc0_950_nHft-1_nsigma300.0_tof999.00_pt0.0/");
    eff.setTofPid();

//    name="/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/pidtest/tpc_"+cutComb+"rootFiles/";
    eff.setFolderTpcPid("pid_files/bbc0_950_nHft0_nsigma3.0_tof0.03_pt0.0/");
    eff.setTpcPid();

    if (eff.isTpcReconstructed("pion", 1, 0, 0.3)) cout<<"yes"<<endl;
    if (eff.isTofmatched("pion", 1, 0.3)) cout<<"yes"<<endl;
    else cout<<"no"<<endl;



    // set input and output
    TFile *f  = new TFile(ntpFileName, "READ");
    auto* ntp = (TNtuple*)f->Get("ntp_signal");

    ntpFastSimBDT* tEvent = new ntpFastSimBDT(ntp);
    long int nEvents = tEvent->fChain->GetEntries();

    TString shortName = ntpFileName.ReplaceAll("root", 4, "", 0);
    cout<<shortName<<endl;

    TString nameSuffixFile = coeffHistoName+fileNameSuffix;
    TFile *dataOu = new TFile(shortName + nameSuffixFile+ ".root","recreate");

//    TNtuple *ntpOut= new TNtuple("ntp_signal","D Meson Tree","D_mass:D_decayL:D_cosThetaStar:cosTheta:D_pt:D_ptSIM:pi1_pt:k_pt:pi1_dca:k_dca:dcaDaughters:dcaD0ToPv:hft:pid:etas:tpc");
    auto* ntpOut = new TNtuple("ntpOut", "ntpOut", "D_mass:D_pt:D_ptSIM:BDTresponse:dcaD0ToPv:dcaDaughters:k_dca:pi1_dca:D_decayL:cosTheta:D_cosThetaStar:pid:hft:etas:tpc:k_pt:pi1_pt");
//    ntpOut->SetDirectory(0);

    const int nNtVars = ntpOut->GetNvar();
    float ntVar[nNtVars];

    for(long int iEvent=0; iEvent<nEvents; iEvent+=1) {
        tEvent->GetEntry(iEvent);
        if (!((iEvent + 1) % ((nEvents) / 5)))
            cout << "________ entries progress = " << iEvent*100 / static_cast<float>(nEvents) << "%" << endl;

        float tpc=0., pid=0.;

        if(eff.isTpcReconstructed("pion", 1, 0, tEvent->pi1_pt) && eff.isTpcReconstructed("kaon", 1, 0, tEvent->k_pt)){
            tpc=1;
        }

        //hybrid TOF
        bool goodTofPion= false;
        bool goodTofKaon= false;

        if (eff.isTofmatched("pion", 0, tEvent->pi1_pt)){
            if (eff.isTofPid("pion", tEvent->pi1_pt)) goodTofPion= true;
            else goodTofPion= false;
        }
        else goodTofPion=true;

        if (eff.isTofmatched("kaon", 0, tEvent->k_pt)){
            if (eff.isTofPid("kaon", tEvent->k_pt)) goodTofKaon= true;
            else goodTofKaon= false;
        }
        else goodTofKaon=true;

        //TPC PID
        bool goodTpcPidPion = false;
        if (eff.isTpcPid("pion", tEvent->pi1_pt)) goodTpcPidPion = true;

        bool goodTpcPidKaon = false;
        if (eff.isTpcPid("kaon", tEvent->k_pt)) goodTpcPidKaon = true;

        if (goodTofPion && goodTofKaon && goodTpcPidPion && goodTpcPidKaon) pid=1.;


        int ii = 0;
        //this stuf are not changing
        ntVar[ii++]=tEvent->D_mass;
        ntVar[ii++]=tEvent->D_pt;
        ntVar[ii++]=tEvent->D_ptSIM;
        ntVar[ii++]=tEvent->BDTresponse;
        ntVar[ii++]=tEvent->dcaD0ToPv;
        ntVar[ii++]=tEvent->dcaDaughters;
        ntVar[ii++]=tEvent->k_dca;
        ntVar[ii++]=tEvent->pi1_dca;
        ntVar[ii++]=tEvent->D_decayL;
        ntVar[ii++]=tEvent->cosTheta;
        ntVar[ii++]=tEvent->D_cosThetaStar;
        ntVar[ii++]=pid;
        ntVar[ii++]=tEvent->hft;
        ntVar[ii++]=tEvent->etas;
        ntVar[ii++]=tpc;
        ntVar[ii++]=tEvent->k_pt;
        ntVar[ii++]=tEvent->pi1_pt;

        //        ntVar[ii++]=tEvent->D_mass;
//        ntVar[ii++]=tEvent->D_decayL;
//        ntVar[ii++]=tEvent->D_cosThetaStar;
//        ntVar[ii++]=tEvent->cosTheta;
//        ntVar[ii++]=tEvent->D_pt;
//        ntVar[ii++]=tEvent->D_ptSIM;
//        ntVar[ii++]=tEvent->pi1_pt;
//        ntVar[ii++]=tEvent->k_pt;
//        ntVar[ii++]=tEvent->pi1_dca;
//        ntVar[ii++]=tEvent->k_dca;
//        ntVar[ii++]=tEvent->dcaDaughters;
//        ntVar[ii++]=tEvent->dcaD0ToPv;
//        ntVar[ii++]=tEvent->hft;
//        ntVar[ii++]=tEvent->pid;
//        //...and these are
//        ntVar[ii++]=tEvent->etas;
//        ntVar[ii++]=tpc;

        ntpOut->Fill(ntVar);
    }


    dataOu->cd();
    ntpOut->Write("ntp_signal", TObject::kOverwrite);

    dataOu->Close();
    delete tEvent;



}
