#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH2.h"
#include "TH1.h"
#include "tmvaCuts.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using tmvaCuts::PtBins;
using tmvaCuts::totalNumberOfEvents;
void makeSignificance(TString, int, int);

void tmvaD0(TString myMethodList = "CutsGA", int ptBin = 1, int pass = 1) {
    cout << "ptBin " << ptBin << endl;
    cout << "pass " << pass << endl;
    float ptmin = PtBins[ptBin];
    float ptmax = PtBins[ptBin + 1];
    cout<<ptmin<<" "<<ptmax<<endl;

    TMVA::Tools::Instance();

    std::map<std::string, int> Use;
    Use["Cuts"] = 1; //1
    Use["CutsD"] = 0; //1
    Use["CutsPCA"] = 0;
    Use["CutsGA"] = 0;
    Use["CutsSA"] = 0;

    std::cout << std::endl;
    std::cout << "==> Start TMVA D0" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector <TString> mlist = TMVA::gTools().SplitString(myMethodList, ',');
        for (UInt_t i = 0; i < mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName("tmvaD0.root");
    TFile *outputFile = TFile::Open(outfileName, "RECREATE");

    TMVA::Factory *factory = new TMVA::Factory("tmvaD0", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

//    factory->AddVariable("D_decayL", 'F');
    dataloader->AddVariable("dcaDaughters", 'F');
    dataloader->AddVariable("dcaD0ToPv", 'F');
    //   factory->AddVariable("cosThetaStar", 'F');
//     factory->AddVariable("cosTheta", 'F');
    dataloader->AddVariable("k_dca", 'F');
    dataloader->AddVariable("pi1_dca", 'F');
//    factory->AddVariable("deltaPt := abs(k_pt-pi1_pt)", 'F');
    dataloader->AddSpectator("D_pt", 'F');
    dataloader->AddSpectator("D_mass", 'F');

    // Read training and test data
    TFile *inputSignal = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.Large.root");
//    TFile *inputSignal = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.root");
    TTree *signal = (TTree *) inputSignal->Get("ntp_signal");
    std::cout << "--- TMVA D0 : Using input signal file: " << inputSignal->GetName() << std::endl;
    dataloader->AddSignalTree(signal, 1);

    TFile *inputBackground = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.picoD0AnaMaker.0802.0415.root");
    TTree *backgroundSameSign = (TTree *) inputBackground->Get("ntp_background");
    dataloader->AddBackgroundTree(backgroundSameSign, 1);

    int nEntriesSignalTree = signal->GetEntries();
    TH1F *hMcPt = (TH1F *) inputSignal->Get("hMcPt");
    int const nOriginalSignalEntriesMCPt = hMcPt->Integral(hMcPt->FindBin(ptmin), hMcPt->FindBin(ptmax)); // Number of simulated D0/D0bar in this pT bin before efficiency (and detector?)
    int const nOriginalSignalEntries = nOriginalSignalEntriesMCPt; // Number of simulated D0/D0bar in this pT bin before efficiency

    TString signalWeightExpression = TString::Format("8*((%f/%f)*2.*0.8*3.14*D_pt*2*(%f)*2.*exp(-1.45-1.73*D_pt)*0.00389*weight)", (float) totalNumberOfEvents, (float) nOriginalSignalEntries, ptmax - ptmin);
    int nEventsBackgroundTree = backgroundSameSign->GetEntries();
    TString backgroundWeightExpression = "1";

    dataloader->SetSignalWeightExpression(signalWeightExpression);
    dataloader->SetBackgroundWeightExpression(backgroundWeightExpression);

    TCut mycuts = Form("D_mass > 1.8 && D_mass < 1.95 && D_pt>%1.2f && D_pt<%1.2f && k_pt>%1.2f && pi1_pt>%1.2f && "
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && "
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f",
                       ptmin, ptmax, tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength[pass], tmvaCuts::dcaDaughters[pass],
                       tmvaCuts::kDca[pass], tmvaCuts::pDca[pass],
                       tmvaCuts::dcaV0ToPv[pass]);

    // check input pt distribution and yield
    TH1F *hPtSignal = new TH1F("hPtSignal", "hPtSignal", 100, 0, 10);
    TH1F *hPtBackgroundSameSign = new TH1F("hPtBackgroundSameSign", "hPtBackgroundSameSign", 100, 0, 10);
    backgroundSameSign->Draw("D_pt>>hPtBackgroundSameSign", backgroundWeightExpression * mycuts, "e");
    signal->Draw("D_pt>>hPtSignal", signalWeightExpression * mycuts, "e");

    outputFile->cd();
    hPtSignal->Write();
    hPtBackgroundSameSign->Write();
    int nSignal = hPtSignal->GetEntries();
    int nBackground = hPtBackgroundSameSign->GetEntries();

    cout << "------------------------------------------------------------------\n";
    cout << "To be used for significance estimation:"<< "\n";
    cout << "Signal counts passed cuts: " << nSignal << "\n";
    cout << "Background counts passed cuts: " << nBackground << endl;
    cout << "------------------------------------------------------------------\n";

    int maxNSignal = 800000;
    int maxNBackground = 10000000;
    if (nSignal > maxNSignal) nSignal = maxNSignal;
    if (nBackground > maxNBackground) nBackground = maxNBackground;

    // Tell the factory how to use the training and testing events
    dataloader->PrepareTrainingAndTestTree(mycuts, nSignal/2, nBackground/2, nSignal/2, nBackground/2);

    if (Use["Cuts"])
        factory->BookMethod(dataloader, TMVA::Types::kCuts, "Cuts", "!H:!V:FitMethod=MC:EffSel:SampleSize=2000000:VarProp=FSmart:CutRangeMin[0]=-10:CutRangeMax[0]=10"); //20 000 000

    if (Use["CutsD"])
        factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsD", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate");

    if (Use["CutsPCA"])
        factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsPCA", "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA");

    if (Use["CutsGA"])
        factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsGA", "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95");

    if (Use["CutsSA"])
        factory->BookMethod(dataloader, TMVA::Types::kCuts, "CutsSA", "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVA D0 is done!" << std::endl;

    delete factory;
    delete dataloader;

    cout << "Nentries in hMcPt in the pT range: " << nOriginalSignalEntriesMCPt << endl;

    // Launch the GUI for the root macros
//    if (!gROOT->IsBatch()) TMVA::TMVAGui(outfileName);

    cout << "------------------------------------------------------------------\n";
    cout << "To be used for significance estimation:"<< "\n";
    cout << "Signal counts passed cuts: " << nSignal << "\n";
    cout << "Background counts passed cuts: " << nBackground << endl;
    cout << "------------------------------------------------------------------\n";

    if (Use["Cuts"]) {
        makeSignificance(outfileName, nSignal, nBackground);
    }
}

void makeSignificance(TString outfileName, int nSignal, int nBackground) {
    TFile *outputFile = TFile::Open(outfileName, "UPDATE");
    TH1F *effS = static_cast<TH1F *>(outputFile->Get("dataset/Method_Cuts/Cuts/MVA_Cuts_effS"));
    TH1F *effB = static_cast<TH1F *>(outputFile->Get("dataset/Method_Cuts/Cuts/MVA_Cuts_effB"));
    TH1F* hSignificance = new TH1F("hSignificance", "Estimated significance", effS->GetNbinsX(), effS->GetMinimum(), effS->GetMaximum());
    Double_t NeffS, NeffB;
    for (int i = 1; i < effS->GetNbinsX()+1; ++i) {
        NeffS = effS->GetBinContent(i);
        NeffB = effB->GetBinContent(i);

        //Standard TMVA significance calculation:
        NeffS *= (double)nSignal;
        NeffB *= (double)nBackground;
        NeffB += NeffS;
        NeffB = sqrt(NeffB);
        if(NeffB==0)         hSignificance->SetBinContent(i, 0);
        else hSignificance->SetBinContent(i, NeffS/NeffB);
    }
    hSignificance->SetStats(0);
    hSignificance->Draw();
    hSignificance->Write("hSignificance");
    outputFile->Close();
    cout<<"Significance plot created for estimated Nsignal and Nbackground"<<endl;
}