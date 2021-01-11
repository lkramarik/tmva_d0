#include "topoComparison.h"
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
#include "tmvaCuts.h"


void projectTopo() {
    gROOT->ProcessLine(".L topoComparison.cpp++");

    Double_t ptMin[]={1,2,3};
    Double_t ptMax[]={2,3,5};

    Double_t nTrees[]={100,150,400};
    Double_t depth[]={3,3,3};

//    Double_t bdtResponseCut[]={0.6, 0.5, 0.4};
//    Double_t bdtResponseCut[]={0.7552, 0.645159, 0.531541};
    Double_t bdtResponseCut[]={0., 0., 0.};

    const int nBins=sizeof(ptMin)/ sizeof(Double_t);

    for (int i = 0; i < nBins; ++i) {
        topoComparison* topo = new topoComparison("hijing");

//        TString fileNameSim7 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_full_D0.toyMc.hijing.nonPrimaryInputs.dca1cm.1208.2.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF7 = new TFile(fileNameSim7 ,"r");
//        auto* ntpSim7 = (TNtuple*)simF7 -> Get("ntp_signal");
////        ntpSim->SetDirectory(0);
//        topo->addNtpToCompare(ntpSim7, "FastSim HIJING primary", "k_pt>0.15 && pi1_pt>0.15");


//        TString fileNameSim = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_FS_hijing_nonPrimary_DCA1cm_goodEvent.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF = new TFile(fileNameSim ,"r");
//        auto* ntpSim = (TNtuple*)simF -> Get("ntp_signal");
////        ntpSim->SetDirectory(0);
//        topo->addNtpToCompare(ntpSim, "FastSim HIJING nonP goodE", "k_pt>0.15 && pi1_pt>0.15");

//
        TString fileNameSim4 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_FS_hijing_nonPrimary_DCA1cm.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
        auto* simF4 = new TFile(fileNameSim4 ,"r");
        auto* ntpSim4 = (TNtuple*)simF4 -> Get("ntp_signal");
        topo->addNtpToCompare(ntpSim4, "FastSim+HIJING", "k_pt>0.15 && pi1_pt>0.15");

        TString fileNameSimaf = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
        auto* simaf = new TFile(fileNameSimaf ,"r");
        auto* ntpSimaf = (TNtuple*)simaf -> Get("ntp_signal");
        topo->addNtpToCompare(ntpSimaf, "FastSim+HIJING D^{0} weight", "k_pt>0.15 && pi1_pt>0.15");
//
//        TString fileNameSimS = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_full_D0.toyMc.dca1cm.HIJING.0608.weight.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simS = new TFile(fileNameSimS ,"r");
//        auto* ntpSimS = (TNtuple*)simS -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSimS, "FastSim HIJING old", "k_pt>0.15 && pi1_pt>0.15", "weight");

//        TString fileNameSim7 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_fullEvent_full_production.vtx.3M.1308.recoVertexD0.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF7 = new TFile(fileNameSim7 ,"r");
//        auto* ntpSim7 = (TNtuple*)simF7 -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSim7, "HIJING reco vertex", "k_pt>0.15 && pi1_pt>0.15");


//        TString fileNameSim3 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_fullEvent_full_production.vtx.3M.1308.AllD0.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF3 = new TFile(fileNameSim3 ,"r");
//        auto* ntpSim3 = (TNtuple*)simF3 -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSim3, "HIJING All", "k_pt>0.15 && pi1_pt>0.15");

//        TString fileNameSim34 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_fullEvent_full_production.vtx.3M.1308.recoVertexD0.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF34 = new TFile(fileNameSim34 ,"r");
//        auto* ntpSim34 = (TNtuple*)simF34 -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSim34, "HIJING refMult<15", "k_pt>0.15 && pi1_pt>0.15 && refMult<15");
//
//
//
//        TString fileNameSim345 = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_fullEvent_full_production.vtx.3M.1308.recoVertexD0.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF345 = new TFile(fileNameSim345 ,"r");
//        auto* ntpSim345 = (TNtuple*)simF345 -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSim345, "HIJING refMult>20", "k_pt>0.15 && pi1_pt>0.15 && refMult>20");
//
//        TString fileNameSim345a = Form("/home/lukas/work/tmva_d0/BDT/pt_%i_%i/n%i_d%i/out_local_SIM_ntp_fullEvent_full_production.vtx.3M.1308.recoVertexD0.root", int(ptMin[i]), int(ptMax[i]), int(nTrees[i]), int(depth[i]) );
//        auto* simF345a = new TFile(fileNameSim345a ,"r");
//        auto* ntpSim345a = (TNtuple*)simF345a -> Get("ntp_signal");
//        topo->addNtpToCompare(ntpSim345a, "HIJING 15<refMult<20", "k_pt>0.15 && pi1_pt>0.15 && refMult>15 && refMult<20");



        topo->setText(Form("%i < D^{0} p_{T} < %i GeV/c", (int)ptMin[i], (int)ptMax[i]));
        topo->setOutFileName(Form("topologyProjections_pt_%i_%i.root",(int)ptMin[i], (int)ptMax[i]));

        topo->addCut(Form("D_pt >= %f && D_pt < %f", ptMin[i], ptMax[i]));
//        topo->addCut("D_mass<2");
//        topo->addCut("D_mass>1.7");

        TCut simToBeSame = "mcEtas>0";
        topo->addCut(simToBeSame);

//        TCut refCut = "refMult>30";
//        topo->addCut(refCut);

        TCut detectorCuts = "hft>0 && tpc>0";
        topo->addCut(detectorCuts);

        TCut precutsTMVA = Form(
                            "k_pt>%1.2f && pi1_pt>%1.2f && "
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
//        topo->addCut(precutsTMVA);
//        topo->addCut(Form("BDTresponse>%f", bdtResponseCut[i]));
//
        topo->setFolderName(Form("pT_%i_%i/", (int)ptMin[i], (int)ptMax[i]));
//
//        topo->addVarToCompare("D_pt", 0, 5, "D_pT");
//        topo.addVarToCompare("BDTresponse", 0, 1, "BDT response");
//        topo->addVarToCompare("k_dca", 0, 0.15, "DCA_{K} [cm]");
//        topo->addVarToCompare("pi1_dca", 0, 0.15, "DCA_{#pi} [cm]");

        topo->addVarToCompare("k_dca", 0, 0.15, "DCA_{K} [cm]", 30, 1);
        topo->addVarToCompare("pi1_dca", 0, 0.15, "DCA_{#pi} [cm]", 30, 1);

//        topo->addVarToCompare("k_dca", 0, 0.05, "DCA_{K} [cm]", 30, 1);
//        topo->addVarToCompare("pi1_dca", 0, 0.05, "DCA_{#pi} [cm]", 30, 1);



//        topo->addVarToCompare("k_dca", 0, 0.2, "DCA_{K} [cm]");
//        topo->addVarToCompare("pi1_dca", 0, 0.2, "DCA_{#pi} [cm]");
//        topo.addVarToCompare("pi1_dca", 0, 0.2, "DCA_{#pi} [cm]");
//    topo.addVarToCompare("k_pt", 0, 3, "p_{TK} [GeV]");
//    topo.addVarToCompare("pi1_pt", 0, 3, "p_{T#pi} [GeV]");
        topo->addVarToCompare("dcaDaughters", 0, 0.015, "DCA daughters [cm]", 30);
        topo->addVarToCompare("D_decayL", 0, 0.1, "decay length [cm]", 30);
//        topo.addVarToCompare("D_decayL", 0, 0.05, "D_decayL [cm]");
        topo->addVarToCompare("dcaD0ToPv", 0, 0.015, "DCA_{D^{0}} [cm]", 30);
        topo->addVarToCompare("cosTheta", 0.5, 1, "cos(#theta)", 30);
        topo->addVarToCompare("D_cosThetaStar", 0.5, 1, "cos(#theta*)", 30);

        topo->project();

//        dataF->Close();
//        simF1->Close();

        delete topo;
    }





}
