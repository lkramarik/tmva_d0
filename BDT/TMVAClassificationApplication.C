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

void TMVAClassificationApplication( const char* inputF = "./../../files_to_run.list", TString output = "out_local.root", float ptmin = 2, float ptmax = 3) {
    cout<<ptmin<<" "<<ptmax<<endl;
    cout<<inputF<<endl;
//    /* input file from input list
    TChain *ntp[2] = {new TChain("ntp_signal","ntp_signal"), new TChain("ntp_background","ntp_background")};
    std::string line;
    std::ifstream infile(inputF);
    TString lineS;
    while (std::getline(infile, line)) {
        cout<<line<<endl;
        lineS = line;
//        theTree -> Add(lineS);
        ntp[0] -> Add(lineS);
        ntp[1] -> Add(lineS);
    }
//    */

//    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.dcaDaughers007.costheta.US.1006.root";
//    TFile* data = new TFile(input ,"r");
//    TFile* sim = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.1605.root" ,"r");
//    TNtuple* ntp[2] = {(TNtuple*)sim -> Get("ntp_signal"), (TNtuple*)data -> Get("ntp_background")};

    TFile* Dplus_file = new TFile (output, "RECREATE");

    TMVA::Tools::Instance();
    std::map<std::string,int> Use;

    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;

    // --- Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    Float_t D_pt, D_mass, refMult, pi1_pt, k_pt, k_dca, pi1_dca, dcaDaughters, cosTheta, D_decayL, dcaD0ToPv, D_cosThetaStar, primary, diffRemovedPrimary;

    reader->AddVariable("k_dca", &k_dca );
    reader->AddVariable("pi1_dca", &pi1_dca );
    reader->AddVariable("dcaDaughters", &dcaDaughters );
    reader->AddVariable("cosTheta", &cosTheta );
    reader->AddVariable("D_decayL", &D_decayL );
    reader->AddVariable("dcaD0ToPv", &dcaD0ToPv );
    reader->AddVariable("D_cosThetaStar", &D_cosThetaStar);

    TString dir    = "dataset/weights/";
//    TString dir    = "weights/";
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

    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

    TStopwatch sw;
    sw.Start();
    TString vars="D_mass:D_pt:BDTresponse:dcaD0ToPv:dcaDaughters:primary:refMult:k_dca:pi1_dca:cosTheta:D_cosThetaStar:D_decayL:precuts";
    TNtuple* ntp_range[2] = {new TNtuple("ntp_signal", "ntp_signal", vars), new TNtuple("ntp_background", "ntp_background", vars)};

    float hodnoty[4] = {0};
    for (int i = 0; i < 2; ++i) {
        ntp[i]->SetBranchAddress("k_pt", &k_pt);
        ntp[i]->SetBranchAddress("pi1_pt", &pi1_pt);
        ntp[i]->SetBranchAddress("k_dca", &k_dca);
        ntp[i]->SetBranchAddress("pi1_dca", &pi1_dca);
        ntp[i]->SetBranchAddress("D_decayL", &D_decayL);
        ntp[i]->SetBranchAddress("D_mass", &D_mass);
        ntp[i]->SetBranchAddress("D_pt", &D_pt);
        ntp[i]->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
        ntp[i]->SetBranchAddress("cosTheta", &cosTheta);
        ntp[i]->SetBranchAddress("primary", &primary);
        ntp[i]->SetBranchAddress("refMult", &refMult);
//        ntp[i]->SetBranchAddress("diffRemovedPrimary", &diffRemovedPrimary);
        ntp[i]->SetBranchAddress("dcaDaughters", &dcaDaughters);
        ntp[i]->SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);

//        for (Long64_t ievt = 0; ievt < 10000; ievt++) {
        for (Long64_t ievt = 0; ievt < ntp[i]->GetEntries(); ievt++) {
            float valueMVA=0., precuts=0.;
            if (ievt % 1000000 == 0 && i == 0) std::cout << "--- ... Processing signal, event: " << ievt << std::endl;
            if (ievt % 1000000 == 0 && i == 1) std::cout << "--- ... Processing background, event: " << ievt << std::endl;
            ntp[i]->GetEntry(ievt);

            if (D_pt < ptmax && D_pt >= ptmin) {
                if  (k_pt>tmvaCuts::minPt && pi1_pt>tmvaCuts::minPt &&
                     D_decayL>tmvaCuts::decayLength && D_decayL<0.2 &&
                     dcaDaughters<tmvaCuts::dcaDaughters &&
                     k_dca>tmvaCuts::kDca && k_dca<1. &&
                     pi1_dca>tmvaCuts::pDca && pi1_dca<1. &&
                     dcaD0ToPv < tmvaCuts::dcaV0ToPv &&
                     cosTheta > tmvaCuts::cosTheta) {
                    precuts=1.;
                }

                    if (Use["BDT"]) {
                        valueMVA = reader->EvaluateMVA("BDT method");
                        histBdt->Fill(valueMVA);
                    }

                ntp_range[i]->Fill(D_mass, D_pt, valueMVA, dcaD0ToPv, dcaDaughters, primary, refMult, k_dca, pi1_dca, cosTheta, D_cosThetaStar, D_decayL, precuts);
//                        ntp_range[i]->Fill(D_mass, D_pt, valueMVA, dcaD0ToPv, dcaDaughters, primary, diffRemovedPrimary);

            }
        }
    }

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    Dplus_file->cd();
    histBdt->Write();
    ntp_range[0] -> Write(ntp_range[0]->GetName(), TObject::kOverwrite);
    ntp_range[1] -> Write(ntp_range[1]->GetName(), TObject::kOverwrite);
    Dplus_file -> Close();

    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
