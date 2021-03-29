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

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
//#include "TInterpreter.h"
//#include "TMVAGui.C"
#include "TMVA/TMVAGui.h"
#include "TMVA/DataLoader.h"

#include "tmvaCuts.h"
#include<fstream>
using namespace std;

void TMVAClassification(float ptmin = 2, float ptmax = 3, float nTrees = 350, float treeDepth = 4) {
   TString inputSignalStr = "/home/lukas/work/tmva_d0/sim/ntp_FS_data_global_HS_dca1_0302.root";

    cout<<ptmin<<" "<<ptmax<<endl;
    const char* inputF = "./../../files_to_run.list";
   TMVA::Tools::Instance();
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( Form("TMVA_bdt_d0_pt_%.1f_%.1f.root", ptmin, ptmax));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   dataloader->AddVariable("k_dca", "Kaon DCA" , "cm", 'F' );
   dataloader->AddVariable("pi1_dca", "Pion DCA" , "cm", 'F' );
   dataloader->AddVariable("dcaDaughters", "Daughters DCA" ,"cm", 'F' );
   dataloader->AddVariable("cosTheta","cos(#theta)" ,"", 'F' );
   dataloader->AddVariable("D_decayL", "Decay length" ,"cm", 'F' );
   dataloader->AddVariable("dcaD0ToPv", "D^{0} DCA" ,"cm", 'F' );
   dataloader->AddVariable("D_cosThetaStar", "cos(#theta^{*})" ,"", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

    Double_t backgroundWeight = 1;

//    TFile style, one file
//   TFile *inputBackground = new TFile("/home/lukas/work/dmesons/Dmaker_dAu/res_analyse/ntp/ntp_lukas_1704.root");
//   TTree *backgroundSameSign = (TTree *) inputBackground->Get("ntp_background");
//   factory->AddBackgroundTree(backgroundSameSign, backgroundWeight);

//   TChain style, more files
    TChain *backgroundSameSign = new TChain("ntp_background","ntp_background");
    std::string line;
    std::ifstream infile(inputF);
    TString lineS;
    while (std::getline(infile, line)) {
        cout<<line<<endl;
        lineS = line;
        backgroundSameSign -> Add(lineS);
    }
    dataloader->AddBackgroundTree(backgroundSameSign, backgroundWeight);

//   TFile *inputBackgroundSide = new TFile("/home/lukas/work/tmva_d0/ntp_2401_sideband.root");
//   TTree *backgroundSideBand = (TTree *) inputBackgroundSide->Get("ntp_sideband");
//   factory->AddBackgroundTree(backgroundSideBand, backgroundWeight);

   TFile *inputSignal = new TFile(inputSignalStr);
   TTree *signal = (TTree *) inputSignal->Get("ntp_signal");
   std::cout << "--- TMVA D0       : Using input signal file: " << inputSignal->GetName() << std::endl;
   Double_t signalWeight = 1.;
   dataloader->AddSignalTree(signal, signalWeight);

   int nEntriesSignalTree = signal->GetEntries();
   int const totalNumberOfEvents = 100e6;
   TH1F*  hMcPt = (TH1F*)inputSignal->Get("hMcPt");
//   TTree *signal_data = (TTree *) inputBackground->Get("ntp_signal");

   int const nOriginalSignalEntriesMCPt = hMcPt->Integral(hMcPt->FindBin(ptmin),hMcPt->FindBin(ptmax)); // Number of simulated D0/D0bar in this pT bin before efficiency
   int const nOriginalSignalEntries = nOriginalSignalEntriesMCPt; // Number of simulated D0/D0bar in this pT bin before efficiency


//    Float_t invYield=7.5*f1D0pp->Eval(pT)*0.61/42;//12: 0.4639 7: 0.4893//  10 is for nBin, 50 is just how it is saved...
//    NSig = 2*TMath::Pi()*pT*dpT*dy*2*invYield*nevt*BR*eff;
//    TString signalWeightExpression = TString::Format("7.5*weight*2.*3.14*D_pt*(%f)*2.*2.*0.0389)",
//                                                     ptmax-ptmin,
//                                                     (float)totalNumberOfEvents,
//                                                     (float)nOriginalSignalEntries,
//                                                     );

//   TString signalWeightExpression = TString::Format("1*weight*((%f/%f)*0.8*2.*3.14*D_pt*2*(%f)*2.*exp(-1.45-1.73*D_pt)*0.0389)", (float)totalNumberOfEvents, (float)nOriginalSignalEntries, ptmax-ptmin);
//   TString signalWeightExpression = TString::Format("((%f/%f)*0.8*2.*3.14*D_pt*2*(%f)*2.*exp(-1.45-1.73*D_pt))", (float)totalNumberOfEvents, (float)nOriginalSignalEntries, ptmax-ptmin);
   TString signalWeightExpression = "1";
//   dataloader->SetSignalWeightExpression(signalWeightExpression);
   TString backgroundWeightExpression = "1";


   // Apply additional cuts on the signal and background samples (can be different)
//   TCut mycuts = "D_pt<3 && D_pt>2 && k_pt>0.15 && pi1_pt>0.15 && k_dca>0.002 && pi1_dca>0.002 && cosTheta>0.5";

   TCut mycuts = Form("D_mass > 1. && D_mass < 3. && D_pt>%1.2f && D_pt<%1.2f && k_pt>%1.2f && pi1_pt>%1.2f && " //_mass > 1.76 && D_mass < 1.96
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && " //lets change to 0.1, was 0.2
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       ptmin, ptmax, tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

   TCut mycutb = mycuts;

   // check input pt distribution and yield
   TH1F *hPtSignal = new TH1F("hPtSignal", "hPtSignal", 100, 0, 10);
   TH1F *hPtBackgroundSameSign = new TH1F("hPtBackgroundSameSign", "hPtBackgroundSameSign", 100, 0, 10);
//   TH1F *hPtBackgroundSideBand = new TH1F("hPtBackgroundSideBand", "hPtBackgroundSideBand", 100, 0, 10);
   backgroundSameSign->Draw("D_pt>>hPtBackgroundSameSign", backgroundWeightExpression * mycuts, "e");
//   backgroundSideBand->Draw("D_pt>>hPtBackgroundSideBand", backgroundWeightExpression * mycuts, "e");
   TH1F *hPtBackground = (TH1F *) hPtBackgroundSameSign->Clone("hPtBackground");
   signal->Draw("D_pt>>hPtSignal", signalWeightExpression * mycuts, "e");
//   hPtBackground->Add(hPtBackgroundSideBand);

    dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

   // Cut optimisation

   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"]) {  // Adaptive Boost
       TString optionsTrees = Form("NTrees=%f:MaxDepth=%f", nTrees, treeDepth);
       TString option = optionsTrees + "!H:!V:MinNodeSize=2.5%:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1";

       //Brieman random forest
//       TString option = "!H:!V:NTrees=400:MaxDepth=7:UseRandomisedTrees:UseNvars=7:BaggedSampleFraction=0.6:nCuts=-1:CreateMVAPdfs";

       factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", option);
   }
   if (Use["BDTB"]) // Bagging
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod(dataloader, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );


   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;

   int nSignal = hPtSignal->GetEntries();
   int nBackground = hPtBackground->GetEntries();
   float const signalYield = hPtSignal->Integral();
   float const backgroundYield = hPtBackground->Integral();

   cout << "Nentries in hMcPt in the pT range: " << nOriginalSignalEntriesMCPt << endl;
//   cout << "Nentries in data signal: " << signal_data->GetEntries(mycuts) << endl;

   cout << endl << "signal total sample yield for sign.: " << signalYield << endl;
   cout << "background total sample yield for sign.: " << backgroundYield << endl << endl;

   cout << "------------------------------------------------------------------\n";
   cout << "Signal counts passed cuts: " << nSignal << "\n";
   cout << "Background counts passed cuts: " << nBackground << endl;
   cout << "Signal counts*weight passed cuts: " << signalYield << "\n";
   cout << "Background counts*weight passed cuts: " << backgroundYield << endl;
   cout << "------------------------------------------------------------------\n";

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
}
