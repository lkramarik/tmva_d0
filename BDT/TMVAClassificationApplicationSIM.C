#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include <math.h>
#if not defined(__CINT__) || defined(__MAKECINT__)
#endif
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "tmvaCuts.h"

using namespace TMVA;

void TMVAClassificationApplicationSIM(TString output = "out_local_SIM.root", float ptmin = 2, float ptmax = 3) {
    cout<<ptmin<<" "<<ptmax<<endl;

    TFile* sim = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.1605.root" ,"r");
    TNtuple* ntp = (TNtuple*)sim -> Get("ntp_signal");

    TFile* Dplus_file = new TFile (output, "RECREATE");

    TMVA::Tools::Instance();
    std::map<std::string,int> Use;

    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;

    // --- Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    Float_t k_pt, pi1_pt, k_dca, pi1_dca, pi2_dca, dcaDaughters, cosTheta, D_decayL, dcaD0ToPv, primary, D_cosThetaStar, diffRemovedPrimary, pid;
    reader->AddVariable("k_dca", &k_dca );
    reader->AddVariable("pi1_dca", &pi1_dca );
    reader->AddVariable("dcaDaughters", &dcaDaughters );
    reader->AddVariable("cosTheta", &cosTheta  );
    reader->AddVariable("D_decayL", &D_decayL );
    reader->AddVariable("dcaD0ToPv", &dcaD0ToPv );
    reader->AddVariable("D_cosThetaStar", &D_cosThetaStar);

    TString dir    = "dataset/weights/";
    TString prefix = "TMVAClassification";

    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
        if (it->second) {
            TString methodName = TString(it->first) + TString(" method");
            TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            reader->BookMVA( methodName, weightfile );
        }
    }

    UInt_t nbin = 100;
    TH1F *histBdt(0);
    if (Use["BDT"]) histBdt = new TH1F("MVA_BDT","MVA_BDT", nbin, -1, 1 );
//    if (Use["BDTD"]) histBdtD = new TH1F("MVA_BDTD","MVA_BDTD", nbin, -0.8, 0.8 );
//    if (Use["BDTG"]) histBdtG = new TH1F("MVA_BDTG","MVA_BDTG", nbin, -1.0, 1.0 );

    Float_t D_mass;
    Float_t D_theta;
    Float_t flag;
    Float_t D_pt;

    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

    TStopwatch sw;
    sw.Start();
    TNtuple* ntp_range = new TNtuple("ntp_signal", "ntp_signal", "D_mass:D_pt:BDTresponse:dcaD0ToPv:dcaDaughters:pid");

    ntp->SetBranchAddress("k_pt", &k_pt);
    ntp->SetBranchAddress("pi1_pt", &pi1_pt);
    ntp->SetBranchAddress("k_dca", &k_dca);
    ntp->SetBranchAddress("pi1_dca", &pi1_dca);
    ntp->SetBranchAddress("D_theta", &D_theta);
    ntp->SetBranchAddress("D_decayL", &D_decayL);
    ntp->SetBranchAddress("D_mass", &D_mass);
    ntp->SetBranchAddress("D_pt", &D_pt);
    ntp->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
    ntp->SetBranchAddress("cosTheta", &cosTheta);
    ntp->SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp->SetBranchAddress("pid", &pid);

    for (Long64_t ievt = 0; ievt < ntp->GetEntries(); ievt++) {
        if (ievt % 1000000 == 0) std::cout << "--- ... Processing signal, event: " << ievt << std::endl;
        ntp->GetEntry(ievt);
        if (D_pt < ptmax && D_pt > ptmin) {
            if  (k_pt>tmvaCuts::minPt && pi1_pt>tmvaCuts::minPt &&
                 D_decayL>tmvaCuts::decayLength && D_decayL<0.2 &&
                 dcaDaughters<tmvaCuts::dcaDaughters &&
                 k_dca>tmvaCuts::kDca && k_dca<0.2 &&
                 pi1_dca>tmvaCuts::pDca && pi1_dca<0.2 &&
                 dcaD0ToPv < tmvaCuts::dcaV0ToPv &&
                 cosTheta > tmvaCuts::cosTheta) {

                if (Use["BDT"]) {
                    float valueMVA = reader->EvaluateMVA("BDT method");
                    histBdt->Fill(valueMVA);
                    ntp_range->Fill(D_mass, D_pt, valueMVA, dcaD0ToPv, dcaDaughters,pid);
                }
            }
        }
    }

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    Dplus_file->cd();
    histBdt->Write();
    ntp_range->Write();
    Dplus_file->Close();

    delete reader;
    std::cout << "==> TMVAClassificationApplicationSIM is done!" << endl << std::endl;
} 
