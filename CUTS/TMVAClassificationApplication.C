#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "fitting.C"
#include "tmvaCuts.h"

using tmvaCuts::PtBins;
using tmvaCuts::totalNumberOfEvents;
using namespace TMVA;

void makeDetailedAnalysis(TTree*, TTree*, TMVA::Reader*, TH1F*);
//Variables need to be defined here, in order to use them in different functions - they need to be somehow connected for data Ntuple and reader
Float_t k_pt, pi1_pt, k_dca, pi1_dca, dcaDaughters, cosTheta, D_decayL, dcaD0ToPv, D_pt, D_mass, ptmin, ptmax;

void TMVAClassificationApplication( TString myMethodList = "", int ptBin = 1, int pass = 1) {
    gROOT->LoadMacro("fitting.C++");

    ptmin = PtBins[ptBin];
    ptmax = PtBins[ptBin + 1];
    cout<<"I will try to apply TMVA you trained for D^{0} p_{T}: "<<ptmin<<" to "<<ptmax<<endl;

    bool makeDetailed = true;
    TMVA::Tools::Instance();
    std::map<std::string, int> Use;

    Use["Cuts"] = 1;
    Use["CutsD"] = 0;
    Use["CutsPCA"] = 0;
    Use["CutsGA"] = 0;
    Use["CutsSA"] = 0;

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
        std::vector <TString> mlist = gTools().SplitString(myMethodList, ',');
        for (UInt_t i = 0; i < mlist.size(); i++) {
            std::string regMethod(mlist[i]);
            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod
                          << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) {
                    std::cout << it->first << " ";
                }
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    // Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used. Their order needs to be same as in the training as well.
//    Float_t k_pt, pi1_pt, k_dca, pi1_dca, dcaDaughters, cosTheta, D_decayL, dcaD0ToPv, D_pt, D_mass;
    reader->AddVariable("dcaDaughters", &dcaDaughters);
    reader->AddVariable("dcaD0ToPv", &dcaD0ToPv);
    reader->AddVariable("k_dca", &k_dca);
//    reader->AddVariable("cosTheta", &cosTheta  );
//    reader->AddVariable("D_decayL", &D_decayL );
    reader->AddVariable("pi1_dca", &pi1_dca);
    reader->AddSpectator("D_pt", &D_pt);
    reader->AddSpectator("D_mass", &D_mass);

    // Book the MVA methods, load results of the training
    TString dir = "dataset/weights/";
    TString prefix = "tmvaD0";
    // Book method(s)
    for (std::map<std::string, int>::iterator it = Use.begin(); it != Use.end(); it++) {
        if (it->second) {
            TString methodName = TString(it->first) + TString(" method");
            TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            reader->BookMVA(methodName, weightfile);
        }
    }

    // Loading the root file with results from the training
    TFile *inputTraining(0);
    TString inputFname = prefix + ".root";
    if (!gSystem->AccessPathName(inputFname)) {
        inputTraining = TFile::Open(inputFname); // check if file in local directory exists
    }
    if (!inputTraining) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }

    TH1F *heffS;
    TH1F *heffB;
    if (Use["Cuts"]) {
        heffS = static_cast<TH1F *>(inputTraining->Get("dataset/Method_Cuts/Cuts/MVA_Cuts_effS"));
        heffB = static_cast<TH1F *>(inputTraining->Get("dataset/Method_Cuts/Cuts/MVA_Cuts_effB"));
    }

    // Select efficiency, that you want to use for getting cuts trained by the classifier
    // Efficiency setter for cut method
    Double_t effS = 0.5;

    if (Use["Cuts"]) {
        TMVA::MethodCuts *methodCuts = reader->FindCutsMVA("Cuts method");
        methodCuts->PrintCuts(effS);
        std::vector <Double_t> cutsMin;
        std::vector <Double_t> cutsMax;
        Double_t true_effS = methodCuts->GetCuts(effS, cutsMin, cutsMax);

        //This is just other type of printing, ready to be used for something else (put cuts in file? plotting significances?)
        //Ok, so it is mainly to crooscheck the order of cuts saved in the arrays
        cout<<"This needs to be same as above, just without fancy stuff:"<<endl;
        for (unsigned short i = 0; i < cutsMin.size(); ++i) {
            cout << cutsMin[i] << " " << cutsMax[i] << endl;
        }
    }

    if (Use["CutsGA"]) {
        TMVA::MethodCuts *mcuts = reader->FindCutsMVA("CutsGA method");
        if (mcuts) {
            std::vector <Double_t> cutsMin;
            std::vector <Double_t> cutsMax;
            mcuts->GetCuts(effS, cutsMin, cutsMax);
            for (unsigned short iGA = 0; iGA < cutsMin.size(); ++iGA) {
                cout << cutsMin[iGA] << " " << cutsMax[iGA] << endl;
            }
        }
    }

    // Prepare input tree (data, stuff, that you want to analyse (on what you want to apply your trained TMVA)
    TFile *input(0);
    TString fname = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.picoD0AnaMaker.0802.0415.root";
    if (!gSystem->AccessPathName(fname)) {
        input = TFile::Open(fname); // check if file in local directory exists
    }

    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVAClassificationApp  : Using input file: " << input->GetName() << std::endl;
    // Loop over the tree, on which you want to apply trained TMVA
    // Prepare the event tree
    // - Here the variable names have to corresponds to your tree
    // - You can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop

    //This makes something for one signal efficiency set of cuts
    if(!makeDetailed){
        TTree *theTree = (TTree *) input->Get("ntp_signal");
        theTree->SetBranchAddress("k_dca", &k_dca);
        theTree->SetBranchAddress("pi1_dca", &pi1_dca);
        theTree->SetBranchAddress("dcaDaughters", &dcaDaughters);
        theTree->SetBranchAddress("cosTheta", &cosTheta);
        theTree->SetBranchAddress("D_decayL", &D_decayL);
        theTree->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);

        std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
        TStopwatch sw;
        sw.Start();

        Int_t nSelCutsGA = 0;
        Int_t nGoodCuts = 0;
        for (Long64_t ievt = 0; ievt < theTree->GetEntries(); ievt++) {
            if (ievt % 10000000 == 0) std::cout << "Processing event: " << ievt << std::endl;
            theTree->GetEntry(ievt);
            if (Use["Cuts"]) {
                Bool_t goodEvent = reader->EvaluateMVA("Cuts method", effS); //cut evaluation: returns 1.0 if event passed, 0.0 otherwise
                if (goodEvent) nGoodCuts++;
            }

            if (Use["CutsGA"]) {
                Bool_t passed = reader->EvaluateMVA("CutsGA method", effS);
                if (passed) nSelCutsGA++;
            }
        }
        sw.Stop();
        std::cout << "End of event loop: "; sw.Print();

        if (Use["Cuts"]) std::cout << "Efficiency for Cuts method (ratio of 'good' to 'all' in analysed tree): " << double(nGoodCuts)/theTree->GetEntries()
                                   << " (for a required signal efficiency of " << effS << ")" << std::endl;

        if (Use["CutsGA"]) std::cout << "Efficiency for CutsGA method (ratio of 'good' to 'all' in analysed tree): " << double(nSelCutsGA)/theTree->GetEntries()
                                     << " (for a required signal efficiency of " << effS << ")" << std::endl;
    }

    //This does everything for all of the signal efficiencies...
    //Caution, this will currently crash for other methods than "Cuts". Everything that says "Cuts method" need to be overwritten.
    if(makeDetailed) makeDetailedAnalysis((TTree*)input->Get("ntp_signal"), (TTree*)input->Get("ntp_background"), reader, heffS);

    // Write histograms
    TFile *target  = new TFile("TMVApp.root","RECREATE");
    //SEM VYPLNIT CO CHCEME?
    heffS->Write();
    heffB->Write();

    target->Close();
    inputTraining->Close();
    input->Close();
    std::cout << "Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
    delete reader;
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

//___________________________________________________________________________________________________
int main( int argc, char** argv ) {
    TString methodList;
    for (int i=1; i<argc; i++) {
        TString regMethod(argv[i]);
        if(regMethod=="-b" || regMethod=="--batch") continue;
        if (!methodList.IsNull()) methodList += TString(",");
        methodList += regMethod;
    }
    TMVAClassificationApplication(methodList);
    return 0;
}

//___________________________________________________________________________________________________
void makeDetailedAnalysis(TTree* signal, TTree* background, TMVA::Reader* reader, TH1F* effS) {
    cout<<"--- Starting the detailed analysis"<<endl;
    TTree* ntp[2] = {signal, background};
    TString name[2] = {"signal", "background"};
    TMVA::MethodCuts *methodCuts = reader->FindCutsMVA("Cuts method");
    TFile *f = new TFile(Form("TMVA_cuts_pt_%.1f_%.1f.root", ptmin, ptmax),"RECREATE");

    const int nEffS = effS->GetNbinsX();
    Float_t effSAxisArray[nEffS];

    //Filling Ntuple with all of the cuts for all of the signall efficiencies
    TNtuple* ntpCuts = new TNtuple("ntpCuts","cutsTuple","effS:k_dca_min:k_dca_max:pi1_dca_min:pi1_dca_max:dcaDaughters_min:dcaDaughters_max:dcaD0ToPv_min:dcaD0ToPv_max");
    const int nNtVars = ntpCuts->GetNvar();
    Float_t ntVar[nNtVars];
    TMVA::MethodCuts *mcuts = reader->FindCutsMVA("Cuts method");

    for (int i = 1; i < effS->GetNbinsX()+1; ++i) {
        effSAxisArray[i-1] = effS->GetBinCenter(i);

        std::vector <Double_t> cutsMin;
        std::vector <Double_t> cutsMax;

        mcuts->GetCuts(effSAxisArray[i-1], cutsMin, cutsMax);
        int ii = 0;
        ntVar[ii++] = effSAxisArray[i-1];
        for (unsigned short iGA = 0; iGA < cutsMin.size(); ++iGA) {
            ntVar[ii++] = cutsMin[iGA];
            ntVar[ii++] = cutsMax[iGA];
        }
        ntpCuts->Fill(ntVar);
    }
    f->cd();
    ntpCuts->Write();

    //Initializing histograms for all of the signal efficiency bins (=> [nEffS]) and signal & background (=> [2])
    TH1F *his[nEffS][2];

    for (int jj = 0; jj < 2; ++jj) {
        for (int ii = 0; ii < nEffS; ii++) {
            his[ii][jj] = new TH1F(Form("D_mass_%.4f_", effSAxisArray[ii]) + name[jj], "D_mass;m [GeV]", 2000, 0.4, 2.4);
            his[ii][jj]->Sumw2();
        }
    }

    methodCuts->PrintCuts(effSAxisArray[30]); //just a test...that reader and everything  is correctly loaded

    //Projecting and evaluating cuts for signal and background  trees
    for (int k = 0; k < 2; ++k) {
        cout<<"--- Evaluating "<<name[k]<<endl;
        ntp[k]->SetBranchAddress("D_mass", &D_mass);
        ntp[k]->SetBranchAddress("D_pt", &D_pt);
        ntp[k]->SetBranchAddress("k_dca", &k_dca);
        ntp[k]->SetBranchAddress("pi1_dca", &pi1_dca);
        ntp[k]->SetBranchAddress("dcaDaughters", &dcaDaughters);
//        ntp[k]->SetBranchAddress("cosTheta", &cosTheta);
//        ntp[k]->SetBranchAddress("D_decayL", &D_decayL);
        ntp[k]->SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);

        for (Long64_t jentry = 0; jentry < ntp[k]->GetEntries(); jentry++) {
            ntp[k]->GetEntry(jentry);
            for (int bin = 0; bin < nEffS; ++bin) {
                Bool_t goodEvent = reader->EvaluateMVA("Cuts method", effSAxisArray[bin]); //cut evaluation: returns ok if event passed, 0.0 otherwise
                if (goodEvent) his[bin][k]->Fill(D_mass);
            }
        }

        f->cd();
        TList *listOut = new TList();
        for (int i = 0; i < nEffS; i++) {
            listOut->Add(his[i][k]);
        }
        listOut->Write("hists_"+name[k], 1, 0);
        delete listOut;
    }

    //Subtractring background from background, fitting and plotting of significances vs signal efficiency,....
    Double_t y[nEffS], x[nEffS], stat[nEffS];

    TH1F *hisSubtr[nEffS];
    TList *listOut = new TList();

    for (int l = 0; l < nEffS; ++l) {
        hisSubtr[l] = new TH1F(Form("D_mass_%.4f", effSAxisArray[l]), "D_mass;m [GeV]", 2000, 0.4, 2.4);
        hisSubtr[l] = (TH1F*)his[l][0]->Clone();
        stat[l] = his[l][0]->Integral(his[l][0] -> FindBin(1.7), his[l][0] -> FindBin(2));
        hisSubtr[l]->Add(his[l][1], -1);
        hisSubtr[l]->Rebin(10);
        listOut->Add(hisSubtr[l]);
        y[l] = fit(his[l][0], his[l][1], ptmin, ptmax, false, true, Form("%.4f", effSAxisArray[l]), effSAxisArray[l]);
        x[l] = effSAxisArray[l];
    }

    f->cd();
    listOut->Write("hists_D_mass", 1, 0);
    f->Close();

    //And let's start plotting stuff:
    gStyle->SetMarkerStyle(34);
    gStyle->SetOptFit(1);
    gStyle->SetStatY(0.899);
    gStyle->SetStatX(0.9);

    TCanvas *cSign = new TCanvas("cSign","cSign",1200,900);
    TGraph* gr = new TGraph(nEffS, x, y);
    gr -> SetMarkerColor(9);
    gr->SetTitle("");
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitle("Significance");
    gr->GetYaxis()->SetTitleOffset(1.1);
    gr->GetXaxis()->SetTitle("Signal efficiency");
    gr->SetMarkerSize(2.5);
    gr->Draw("ap");

    TFile *fSign = new TFile(Form("significance_pt_%.1f_%.1f.root", ptmin, ptmax),"RECREATE");
    gr->Write(Form("gr_sign_pt_%.1f_%.1f", ptmin, ptmax));
    fSign->Close();

    TLatex txR;
    txR.SetNDC();
    txR.SetTextSize(0.04);
    txR.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text1->SetTextSize(0.04);
    text1->SetShadowColor(0);
    text1->SetLineColor(0);
    text1->SetFillColor(0);
    text1->AddText("d+Au 200 GeV");
    text1->Draw("same");

    TCanvas *cStat = new TCanvas("cStat","cStat",1200,900);
    cStat -> SetLogy();
    TGraph* grStat = new TGraph(nEffS, x, stat);
    grStat->SetMarkerColor(46);
    grStat->SetTitle("");
    grStat->GetXaxis()->SetLabelSize(0.04);
    grStat->GetYaxis()->SetLabelSize(0.04);
    grStat->GetXaxis()->SetTitleSize(0.045);
    grStat->GetYaxis()->SetTitleSize(0.045);
    grStat->GetYaxis()->SetTitle("# pairs with mass [1.7, 2.0] GeV/c^{2}");
    grStat->GetYaxis()->SetTitleOffset(1.1);
    grStat->GetXaxis()->SetTitle("Signal efficiency");
    grStat->SetMarkerSize(2.5);
//    TGraphErrors* gr = new TGraph(n_bin, x, y);
    grStat -> Draw("ap");
    text1 -> Draw("same");
    txR.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

    //Fuck, it's just too long :(....
    cout<<"--- End the detailed analysis"<<endl;
}