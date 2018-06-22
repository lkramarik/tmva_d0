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

void TMVAClassificationApplication( const char* inputF = "./../files_to_run.list", TString output = "out_local.root", float ptmin = 2, float ptmax = 3) {
    cout<<ptmin<<" "<<ptmax<<endl;

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

//    TString input = "/home/lukas/work/dmesons/Dmaker_dAu/res_analyse/ntp/ntp_lukas_1704.root";
//    TFile* data = new TFile(input ,"r");
//    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_signal"), (TNtuple*)data -> Get("ntp_background")};

    TFile* Dplus_file = new TFile (output, "RECREATE");

    #ifdef __CINT__
        gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
    #endif

    TMVA::Tools::Instance();
    std::map<std::string,int> Use;

    // --- Cut optimisation
    Use["Cuts"]            = 0;
    Use["CutsD"]           = 0;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;

    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    //
    // --- Support Vector Machine
    Use["SVM"]             = 0;
    //
    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    //
    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 0;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;

    // --- Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    Float_t k_pt, pi1_pt, k_dca, pi1_dca, pi2_dca, dcaDaughters, cosTheta, D_decayL, dcaD0ToPv;
    reader->AddVariable("k_dca", &k_dca );
    reader->AddVariable("pi1_dca", &pi1_dca );
    reader->AddVariable("dcaDaughters", &dcaDaughters );
    reader->AddVariable("cosTheta", &cosTheta  );
    reader->AddVariable("D_decayL", &D_decayL );
    reader->AddVariable("dcaD0ToPv", &dcaD0ToPv );

    TString dir    = "weights/";
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
    TNtuple* ntp_range[2] = {new TNtuple("ntp_signal", "ntp_signal", "D_mass:D_pt:BDTresponse:dcaD0ToPv:dcaDaughters"), new TNtuple("ntp_background", "ntp_background", "D_mass:D_pt:BDTresponse:dcaD0ToPv:dcaDaughters")};

    float hodnoty[4] = {0};
    for (int i = 0; i < 2; ++i) {
        ntp[i]->SetBranchAddress("k_pt", &k_pt);
        ntp[i]->SetBranchAddress("pi1_pt", &pi1_pt);
        ntp[i]->SetBranchAddress("k_dca", &k_dca);
        ntp[i]->SetBranchAddress("pi1_dca", &pi1_dca);
        ntp[i]->SetBranchAddress("D_theta", &D_theta); //D_theta only in pico, no cos, need to copy in the event loop
        ntp[i]->SetBranchAddress("D_decayL", &D_decayL);
        ntp[i]->SetBranchAddress("D_mass", &D_mass);
        ntp[i]->SetBranchAddress("D_pt", &D_pt);
        ntp[i]->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
        ntp[i]->SetBranchAddress("cosTheta", &cosTheta);
        ntp[i]->SetBranchAddress("dcaDaughters", &dcaDaughters);

//        for (Long64_t ievt = 0; ievt < 10000; ievt++) {
        for (Long64_t ievt = 0; ievt < ntp[i]->GetEntries(); ievt++) {
            if (ievt % 100000 == 0 && i == 0) std::cout << "--- ... Processing signal, event: " << ievt << std::endl;
            if (ievt % 100000 == 0 && i == 1) std::cout << "--- ... Processing background, event: " << ievt << std::endl;
            ntp[i]->GetEntry(ievt);
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
                        ntp_range[i]->Fill(D_mass, D_pt, valueMVA, dcaD0ToPv, dcaDaughters);
                    }
                }
            }
        }
    }

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    Dplus_file->cd();
    histBdt->Write();
    ntp_range[0] -> Write();
    ntp_range[1] -> Write();
    Dplus_file -> Close();

    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
