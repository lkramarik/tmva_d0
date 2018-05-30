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
#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "" ) {
   float ptmin = 2;
   float ptmax = 3;

   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
   TString thisdir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   gROOT->SetMacroPath(thisdir + ":" + gROOT->GetMacroPath());
   gROOT->ProcessLine(".L TMVAGui.C");

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA_Dpm_pt2.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   factory->AddVariable("k_dca", "DCA (Kaon)" , 'F' );
   factory->AddVariable("pi1_dca","DCA (Pion1)" , 'F' );
   factory->AddVariable("dcaDaughters", "DCA_daughters" , 'F' );
   factory->AddVariable("cosTheta","cos(#theta)" , 'F' );
//   factory->AddVariable("D_decayL", "#lambda" , 'F' );
   factory->AddVariable("dcaD0ToPv", "dcaD0ToPv" , 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   TFile *inputBackground = new TFile("/home/lukas/work/dmesons/Dmaker_dAu/res_analyse/ntp/ntp_lukas_1704.root");
   TTree *backgroundSameSign = (TTree *) inputBackground->Get("ntp_background");
   TFile *inputBackgroundSide = new TFile("/home/lukas/work/tmva_d0/ntp_2401_sideband.root");
   TTree *backgroundSideBand = (TTree *) inputBackgroundSide->Get("ntp_sideband");
   Double_t backgroundWeight = 1;
   factory->AddBackgroundTree(backgroundSameSign, backgroundWeight);
   factory->AddBackgroundTree(backgroundSideBand, backgroundWeight);

   TFile *inputSignal = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.Large.root");
//    TFile *inputSignal = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.root");
   TTree *signal = (TTree *) inputSignal->Get("ntp_signal");
   std::cout << "--- TMVA D0       : Using input signal file: " << inputSignal->GetName() << std::endl;
   Double_t signalWeight = 1.;
   factory->AddSignalTree(signal, signalWeight);

   int nEntriesSignalTree = signal->GetEntries();
   int const totalNumberOfEvents = 120e6;
   TH1F*  hMcPt      = (TH1F*)inputSignal->Get("hMcPt");
   int const nOriginalSignalEntries = hMcPt->Integral(hMcPt->FindBin(ptmin),hMcPt->FindBin(ptmax)); // Number of simulated D0/D0bar in this pT bin before efficiency
   TString signalWeightExpression = TString::Format("1*weight*((%f/%f)*0.8*2.*3.14*D_pt*2*(%f)*2.*exp(-1.45-1.73*D_pt)*0.0389)", (float)totalNumberOfEvents, (float)nOriginalSignalEntries, ptmax-ptmin);
   factory->SetSignalWeightExpression(signalWeightExpression);

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "D_pt<3 && D_pt>2 && k_pt>0.15 && pi1_pt>0.15 && k_dca>0.002 && pi1_dca>0.002 && cosTheta>0.5"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = mycuts; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the factory how to use the training and testing events
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );


   // ---- Now you can tell the factory to train, test, and evaluate the MVAs
   // Train MVAs using the set of training events
   factory->TrainAllMethods();
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}