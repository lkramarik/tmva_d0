#include "tmvaCuts.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"
#include "TLatex.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TCollection.h"
#include <iostream>
#include <iostream>
#include <fstream>
#include "analyse/FitD0Peak.hh"

using namespace std;
Double_t NSig = -1;
Long64_t Nbckg = -1;

void estimateNsigNbckg(Double_t, Double_t);
Double_t doSignificance(double, double, int, int);
void doRawYieldSim(Double_t, Double_t, int, int, Double_t);
void doRawYield(Double_t, Double_t, int, int, Double_t);
void project_bdt_oneCut(Double_t, Double_t, Double_t, Double_t, Double_t, TString);
void project_bdt_oneCut_SIM(Double_t, Double_t, Double_t, Double_t, Double_t, TString);
void plotTogether(Int_t, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *);
std::vector<TCut> mCuts;
TCut mPidCut;

Double_t doSignificance(double ptmin, double ptmax, int nTrees, int maxDepth) {
    TFile *inputBckg = new TFile(Form("pt_%.0f_%.0f/n%i_d%i/TMVA_bdt_d0_pt_%.1f_%.1f.root", ptmin, ptmax, nTrees, maxDepth, ptmin, ptmax));

    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.03);

    TH1D *hS = (TH1D*) inputBckg->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffS");
    hS->SetLineColor(4);
    hS->SetStats(0);
    hS->GetYaxis()->SetTitle("Efficiency");
    hS->GetYaxis()->SetTitleOffset(0.93);
    hS->GetXaxis()->SetTitle("BDT Response Cut");
    hS->GetXaxis()->SetTitleSize(0.04);
    hS->GetYaxis()->SetTitleSize(0.04);
    hS->SetTitle("");
    hS->SetLineWidth(2);

    TH1D *hB = (TH1D*) inputBckg->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffB");
    hB->SetLineColor(2);
    hB->SetStats(0);
    hB->SetTitle("");
    hB->SetLineWidth(2);

    //////////////////////////////////////////////////////
    TCanvas *cEff = new TCanvas("cEff", "cEff", 900, 1000);
    cEff->SetGrid();

    hS->Draw();
    hB->Draw("same");
    tx2.DrawLatex(0.1,0.92,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

    TLegend *legend = new TLegend(0.663, 0.784, 0.813, 0.89);
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(hS, "Signal", "pl");
    legend->AddEntry(hB, "Background", "pl");
    legend->Draw("same");

    cEff->SaveAs(Form("finalAnalysis/img/BDT_eff_%.1f_%.1f_n%i_d%i.png", ptmin, ptmax, nTrees, maxDepth));
    cEff->SaveAs(Form("finalAnalysis/img/BDT_eff_%.1f_%.1f_n%i_d%i.pdf", ptmin, ptmax, nTrees, maxDepth));

    //////////////////////////////////////////////////////

    TH1D *hSign = (TH1D*) hS->Clone("hSign");
    hSign->SetDirectory(0);
    hSign->SetLineColor(8);
    hSign->GetYaxis()->SetTitle("Significance");
    hSign->GetXaxis()->SetTitle("BDT Response Cut");
    hSign->SetTitle("");

    estimateNsigNbckg(ptmin, ptmax);

    if (NSig<0 || Nbckg<0) return -999;

    Double_t S, B;
    for (int i = 1; i < hS->GetNbinsX()+1; ++i) {
        S = hS->GetBinContent(i);
        B = hB->GetBinContent(i);
        hSign->SetBinContent(i, NSig*S/(sqrt(NSig*S + Nbckg*B)));
    }

    Double_t maxSign=hSign->GetBinCenter(hSign->GetMaximumBin());
    cout<<maxSign<<endl;
    TF1 *signFct = new TF1("signFct", "[0]+[1]*x-[2]*(x-[3])*(x-[3])", maxSign-0.05, maxSign+0.05);
    signFct->SetParameters(5, 0, 100, maxSign);

    ////////////////////////////////////////////////////////////
    TCanvas *cSign = new TCanvas("cSign", "cSign", 900, 1000);
    cSign->SetGrid();
    hSign->Fit("signFct","RLL");

    Double_t bdtResponse = signFct->GetMaximumX();
    cout<<bdtResponse<<endl;

    hSign->Draw();
    tx2.DrawLatex(0.1,0.92,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

    TPaveText *text5 = new TPaveText(0.114,0.816,0.625,0.89,"brNDC");
    text5->SetTextSize(0.025);
    text5->SetLineColor(0);
    text5->SetShadowColor(0);
    text5->SetFillStyle(0);
    text5->SetTextAlign(11);
    text5->AddText(Form("Maximum at BDT reponse cut: %3.4f", bdtResponse));
    text5->AddText(Form("For number of signal (background): %.0f (%lli)", NSig, Nbckg));
    text5->Draw("same");

    TLine *leftline2 = new TLine(bdtResponse, hSign->GetMinimum(), bdtResponse, 1.05*hSign->GetMaximum());
    leftline2->SetLineStyle(9);
    leftline2->SetLineWidth(2);
    leftline2->SetLineColor(2);
    leftline2->Draw("same");

    cSign->SaveAs(Form("finalAnalysis/img/BDT_significance_%.1f_%.1f_n%i_d%i.png", ptmin, ptmax, nTrees, maxDepth));
    cSign->SaveAs(Form("finalAnalysis/img/BDT_significance_%.1f_%.1f_n%i_d%i.pdf", ptmin, ptmax, nTrees, maxDepth));

    cSign->Close();
    cEff->Close();

    ////////////////////////////////////////////////////////////
    ofstream myfile;
    myfile.open(Form("finalAnalysis/BDT_cut_%.1f_%.1f_n%i_d%i.txt", ptmin, ptmax, nTrees, maxDepth));
    myfile << bdtResponse;
    myfile << "\n";
    myfile << signFct->Eval(bdtResponse);
    myfile.close();

    return bdtResponse;
}

//_______________________________________________________________________
void estimateNsigNbckg(Double_t ptmin = 3, Double_t ptmax = 5){
    Double_t dpT=ptmax-ptmin;
    Double_t pT=(ptmin+ptmax)/2;
    Double_t dy=2;
    Double_t invYield=0.0001; //D0 AuAu 50-80 %, Ncoll = cca 40 -> 0.0002, pt23
//    Double_t invYield=0.00015; //D0 AuAu 50-80 %, Ncoll = cca 40 -> 0.0002, pt12
//    Double_t invYield=0.00002; //D0 AuAu 50-80 %, Ncoll = cca 40 -> 0.0002, pt35
//    invYield=invYield*7/50;
    Double_t nevt = 140e6;
    Double_t BR=0.0395;

    Double_t invYieldPubl[] = {0.00678, 0.00329, 0.00216, 0.000217};
    Double_t invYieldPublErr[] = {0.00193, 0.00102, 0.00055, 0.000094};
    Double_t ptPubl[] = {0.3, 0.75, 1.25, 2.25};

    TGraphErrors* grInvYieldsPubl = new TGraphErrors(4, ptPubl, invYieldPubl, 0, invYieldPublErr);
    grInvYieldsPubl->SetMarkerColor(1);
    grInvYieldsPubl->SetMarkerStyle(21);
    grInvYieldsPubl->SetTitle("");
    grInvYieldsPubl->GetXaxis()->SetLabelSize(0.04);
    grInvYieldsPubl->GetYaxis()->SetLabelSize(0.04);
    grInvYieldsPubl->GetXaxis()->SetTitleSize(0.045);
    grInvYieldsPubl->GetXaxis()->SetLimits(0, 5);
    grInvYieldsPubl->GetYaxis()->SetLimits(0.00000001, 0.01);
    grInvYieldsPubl->GetYaxis()->SetTitleSize(0.045);
    grInvYieldsPubl->GetYaxis()->SetTitle("Inv. yield D^{0} d+Au 2004");
    grInvYieldsPubl->GetYaxis()->SetTitleOffset(1.6);
    grInvYieldsPubl->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    grInvYieldsPubl->SetMarkerSize(2.5);

//    TF1 *fD0 = new TF1("fD0", "[0]+x**[1]", 0.2, 3);
//    TF1 *fD0 = new TF1("fD0", "[0]+[1]*x+[2]*x*x+[3]*x*x*x*x", 0, 10);
    TF1 *fD0 = new TF1("fD0", "1/(2*pi)*[dndy]*([A]-1)*([A]-2)/([A]*[B]*([A]*[B]+1.864*([A]-2)))*pow(1+(sqrt(x*x+1.864*1.864)-1.864)/([A]*[B]),-[A])", 0.2, 6);
//    fD0->SetParLimits(1, 0.001, 0.5);
    fD0->SetParameter(1, 0.3);
    fD0->SetParLimits(2, 0.001, 0.05);
    fD0->SetParameter(2, 0.028);
    fD0->SetParameter(0,-1);
    fD0->SetParameter(1,-0.003);
//    TF1 *fD0 = new TF1("fD0", "[0]+TMath::Exp(-[1]*x)", 0, 5);
    grInvYieldsPubl->Fit(fD0, "R");

    cout<<fD0->Eval(1.5)<<" "<<fD0->Eval(2.5)<<" "<<fD0->Eval(4)<<endl;

    TCanvas *cInvYields = new TCanvas("invYielsPubl", "invYielsPubl", 900, 1000);
    gPad->SetLeftMargin(0.15);
    cInvYields->SetLogy();
    grInvYieldsPubl->Draw("ap");
    grInvYieldsPubl->Draw("X");

    cInvYields->SaveAs("finalAnalysis/grInvYieldsPubl.png");
    cInvYields->SaveAs("finalAnalysis/grInvYieldsPubl.pdf");
    cInvYields->Close();

    TCut mycuts = Form("D_mass > 1.815 && D_mass < 1.905 && D_pt>%1.2f && D_pt<%1.2f && k_pt>%1.2f && pi1_pt>%1.2f && "
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && "
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       ptmin, ptmax, tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

    TFile *inputSim = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMc.1605.root");
    TTree *signalTree = (TTree*) inputSim->Get("ntp_signal");
    Long64_t signal=signalTree->GetEntries(mycuts);
//    Long64_t all=signalTree->GetEntries();
    Long64_t all=signalTree->GetEntries("D_mass > 1.815 && D_mass < 1.905");

    cout<<(double)signal/(double)all<<endl;
//    Double_t eff=(double)signal/(double)all;

    TF1 *reco = new TF1("reco", "0.5*([0]+[1]*x+[2]*x*x+[3]*x*x*x)", 0, 10);
//    reco->SetParameters(0.0239696, -0.0503129, 0.0556806, -0.0041924);
    reco->SetParameters(0.00174319, -0.00388281, 0.0048162, -0.000424064);
    Double_t totaleff=reco->Eval(pT);
//    Double_t topo = 0.05; //from AuAu paper pt23
//    Double_t topo = 0.01; //from AuAu paper pt12
    Double_t topo = 0.1; //from AuAu paper pt35
    Double_t deteEff = totaleff/topo;
    Double_t eff = deteEff*(double)signal/(double)all;
//    eff=0.02;
    cout<<eff<<endl;
    NSig = 2*TMath::Pi()*pT*dpT*dy*2*invYield*nevt*BR*eff;

    cout<<"NSIG = "<<NSig<<endl;

    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.2506.root";
    TFile *inputBckg = new TFile(input);
    TTree *bckgTree = (TTree*) inputBckg->Get("ntp_background");

    Nbckg=bckgTree->GetEntries(mycuts);

    cout<<"Nbckg = "<<Nbckg<<endl;

    inputSim->Close();
    inputBckg->Close();

}

//_______________________________________________________________________
void doRawYield(Double_t ptmin, Double_t ptmax, int nTrees, int maxDepth, Double_t bdtCut) {
    TString inputFile=Form("/home/lukas/work/tmva_d0/BDT/finalAnalysis/out_local_pt_%.1f_%.1f_n%i_d%i.root", ptmin, ptmax, nTrees, maxDepth);
    gSystem->Exec(Form("cp /home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/out_local.root %s", ptmin, ptmax, nTrees, maxDepth, inputFile.Data()));
    project_bdt_oneCut(ptmin, ptmax, (Double_t)nTrees, (Double_t)maxDepth, bdtCut, inputFile);
}

//_______________________________________________________________________
void doRawYieldSIM(Double_t ptmin, Double_t ptmax, int nTrees, int maxDepth, Double_t bdtCut) {
    TString inputFile=Form("/home/lukas/work/tmva_d0/BDT/finalAnalysis/out_local_SIM_pt_%.1f_%.1f_n%i_d%i.root", ptmin, ptmax, nTrees, maxDepth);
    gSystem->Exec(Form("cp /home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/out_local_SIM.root %s", ptmin, ptmax, nTrees, maxDepth, inputFile.Data()));
    project_bdt_oneCut_SIM(ptmin, ptmax, (Double_t)nTrees, (Double_t)maxDepth, bdtCut, inputFile);
}

//_________________________________________________________________________________________________________________________________
void project_bdt_oneCut(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange, TString input) {
    bool mixed=false;
    TFile *f = new TFile(Form("finalAnalysis/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtRange, ptmin, ptmax, nTrees ,maxDepth),"RECREATE");

    TH1F *his[2];
    TString name[2] = {"signal", "background"};

    TFile* data = new TFile(input ,"r");
    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_"+name[0]), (TNtuple*)data -> Get("ntp_"+name[1])};

    float D_mass, D_pt, BDTresponse;

    TCut setCuts = "";
    for(unsigned int k = 0; k < mCuts.size(); ++k) {
        setCuts += mCuts[k];
    }
    setCuts += mPidCut;
    cout<<setCuts<<endl;

    for (int k = 0; k < 2; ++k) {
        his[k] = new TH1F(Form("D_mass_%.4f_", bdtRange) + name[k], "D_mass;m [GeV]", 2000, 0.4, 2.4);
        his[k]->Sumw2();
        ntp[k]->Project(his[k]->GetName(), "D_mass", setCuts);

        f->cd();
        TList *listOut = new TList();
        listOut->Add(his[k]);
        listOut->Write("hists_"+name[k], 1, 0);
        delete listOut;
    }

    if (mixed){
        TH1F *hisMxd = new TH1F(Form("D_mass_%.4f_ME", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
        hisMxd -> Sumw2();

        TFile* dataMxd = new TFile("../out_local_mix.root" ,"r");
        TNtuple* ntpMxd = (TNtuple*)dataMxd -> Get("ntp_ME");
        ntpMxd->Project(hisMxd->GetName(), "D_mass", setCuts);
        f -> cd();
        TList *listOut = new TList();
        listOut->Add(hisMxd);
        listOut->Write("hists_ME", 1, 0);
        delete listOut;
    }

    TList *listOut = new TList();

    TH1F *hisSubtr = new TH1F(Form("D_mass_%.4f", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
    hisSubtr = (TH1F*)his[0]->Clone();
    Double_t stat = his[0] -> Integral(his[0] -> FindBin(1.7), his[0] -> FindBin(2));
    hisSubtr->Add(his[1], -1);
    hisSubtr->Rebin(10);
    listOut->Add(hisSubtr);

    TString outFile=Form("finalAnalysis/signals_%.4fbdt_pt_%.1f_%.1f.root", bdtRange, ptmin, ptmax);
    FitD0Peak *fitmass = new FitD0Peak(his[0], his[1], ptmin, ptmax, outFile);
    fitmass->doStuff();

    Float_t ptBin=(ptmax+ptmin)/2;
//    TCanvas *cEff = new TCanvas("cEff", "cEff", 900, 1000);
//    cEff->SetGrid();
    TGraphErrors* grResults = new TGraphErrors();
    grResults->SetPoint(0, ptBin, fitmass->getRawYield());
    grResults->SetPointError(0, 0, fitmass->getRawYieldError());
    grResults->SetMarkerColor(1);
    grResults->SetMarkerStyle(21);
    grResults->SetTitle("");
    grResults->GetXaxis()->SetLabelSize(0.04);
    grResults->GetYaxis()->SetLabelSize(0.04);
    grResults->GetXaxis()->SetTitleSize(0.045);
    grResults->Draw("ap");

    delete fitmass;

    f->cd();
    grResults->Write("raw_yield");
//    grResults->Write("significance");
    listOut->Write("hists_D_mass", 1, 0);
    data->Close();
    f->Close();
}

//_________________________________________________________________________________________________________________________________
void project_bdt_oneCut_SIM(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange, TString input) {
    TFile *f = new TFile(Form("finalAnalysis/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtRange, ptmin, ptmax, nTrees ,maxDepth),"RECREATE");

    TFile* data = new TFile(input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get("ntp_signal");

    TH1F *his = new TH1F(Form("D_mass_%.4f_signal", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
    TH1F *hisEmpty = new TH1F(Form("D_mass_%.4f_background", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4); //fake background for FitD0Peak
    his -> Sumw2();

    TCut setCuts = "";
    for(unsigned int k = 0; k < mCuts.size(); ++k) {
        setCuts += mCuts[k];
    }

//    TCut=Form("",); //reconstructed D0

    Long64_t nAllSim=ntp->GetEntries(Form("D_pt>=%.1f && D_pt<%.1f", ptmin, ptmax)); //all from simu
    setCuts+="pid>0";
    Long64_t nAllPid=ntp->GetEntries(setCuts); //pid cuts
    setCuts+=mPidCut;

    cout<<nAllPid<<endl;

    cout<<setCuts<<endl;
    ntp->Project(his->GetName(), "D_mass", setCuts); //same as in the data
    Long64_t nDataReco=ntp->GetEntries(setCuts); //same as in the data
    cout<<nDataReco<<endl;
    f->cd();
    TList *listOut = new TList();
    listOut->Add(his);
    listOut->Write("hists_signal", 1, 0);
    delete listOut;

    TString outFile=Form("finalAnalysis/signals_SIM_%.4fbdt_pt_%.1f_%.1f.root", bdtRange, ptmin, ptmax);
    FitD0Peak *fitmass = new FitD0Peak(his, hisEmpty, ptmin, ptmax, outFile);
    fitmass->doStuff();

    Float_t ptBin=(ptmax+ptmin)/2;
    Float_t ptBinWidth=(ptmax-ptmin)/2;
    TGraphErrors* grResults = new TGraphErrors();
    grResults->SetPoint(0, ptBin, fitmass->getRawYield());
    grResults->SetPointError(0, 0, fitmass->getRawYieldError());
    grResults->SetMarkerColor(1);
    grResults->SetMarkerStyle(21);
    grResults->SetTitle("");
    grResults->GetXaxis()->SetLabelSize(0.04);
    grResults->GetYaxis()->SetLabelSize(0.04);
    grResults->GetXaxis()->SetTitleSize(0.045);

    TGraphErrors* grTpcAcc = new TGraphErrors();
    grTpcAcc->SetPoint(0, ptBin, (float)nAllPid/(float)nAllSim);
    grTpcAcc->SetPointError(0, ptBinWidth, 0);

    TGraphErrors* grTpcAccHftPidBDT = new TGraphErrors();
    grTpcAccHftPidBDT->SetPoint(0, ptBin, (float)nDataReco/(float)nAllSim);
    grTpcAccHftPidBDT->SetPointError(0, ptBinWidth, 0);

    delete fitmass;

    f->cd();
    grResults->Write("raw_yield");
    grTpcAcc->Write("grTpcAcc");
    grTpcAccHftPidBDT->Write("grTpcAccHftPidBDT");
//    grResults->Write("significance");
    data->Close();
    f->Close();
}

//_______________________________________________________________________
void correctYield(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange) {
    TFile *fData = new TFile(Form("finalAnalysis/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtRange, ptmin, ptmax, nTrees ,maxDepth),"READ");
    TFile *fSIM = new TFile(Form("finalAnalysis/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtRange, ptmin, ptmax, nTrees ,maxDepth),"READ");

    TGraphErrors *grRawYield = (TGraphErrors*) fData->Get("raw_yield");

    fData->Close();
    fSIM->Close();
}

//_______________________________________________________________________
void plotTogether(Int_t nBins, Double_t *ptmin, Double_t *ptmax, Double_t *nTrees, Double_t *maxDepth, Double_t *bdtResponse){
    TGraphErrors* gr = new TGraphErrors();
    Double_t rawYields[nBins], rawYieldsErr[nBins];
    Double_t pT[nBins], ptWidths[nBins], y[nBins], yE[nBins];

    TFile *fOut = new TFile("finalAnalysis/final_result.root", "RECREATE");

    for (int i = 0; i < nBins; ++i) {
        TFile *fData = new TFile(Form("finalAnalysis/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtResponse[i], ptmin[i], ptmax[i], nTrees[i], maxDepth[i]), "READ");
        gr = (TGraphErrors *) fData->Get("raw_yield");
        gr->GetPoint(0, pT[i], rawYields[i]);
        ptWidths[i] = gr->GetErrorX(0);
        rawYieldsErr[i] = gr->GetErrorY(0);
        fData->Close();
    }

    TGraphErrors* grYields = new TGraphErrors(nBins, pT, rawYields, ptWidths, rawYieldsErr);
    grYields->SetMarkerColor(1);
    grYields->SetMarkerStyle(21);
    grYields->SetTitle("");
    grYields->GetXaxis()->SetLabelSize(0.04);
    grYields->GetYaxis()->SetLabelSize(0.04);
    grYields->GetXaxis()->SetTitleSize(0.045);
    grYields->GetYaxis()->SetTitle("Raw yield");
    grYields->GetYaxis()->SetTitleOffset(1.6);
    grYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fOut->cd();
    grYields->Write("raw_yields_all");

    std::vector<TGraphErrors*> graphsSIM;
    Int_t colors[] = {1, 46, 9, 9, 40, 41, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    TString names[2] = {"grTpcAcc", "grTpcAccHftPidBDT"};

    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < nBins; ++i) {
            TFile *fData = new TFile(Form("finalAnalysis/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtResponse[i], ptmin[i], ptmax[i], nTrees[i], maxDepth[i]), "READ");
            gr = (TGraphErrors *) fData->Get(names[j]);
            gr->GetPoint(0, pT[i], y[i]);
            ptWidths[i] = gr->GetErrorX(0);
            yE[i] = gr->GetErrorY(0);
            fData->Close();
        }
        TGraphErrors* grP = new TGraphErrors(nBins, pT, y, ptWidths, yE);
        graphsSIM.push_back(grP);
        fOut->cd();
        grP->Write(names[j]);

    }
    TCanvas *out = new TCanvas("out", "out", 1200, 1000);

    TMultiGraph *mg = new TMultiGraph();
    mg->GetYaxis()->SetTitle("#varepsilon");
    mg->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mg->GetXaxis()->SetTitleSize(0.96);
    mg->GetYaxis()->SetLimits(0, 1);
    mg->SetTitle("");
    TLegend *legend = new TLegend(0.126, 0.71, 0.277, 0.85);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    for (unsigned short j = 0; j < graphsSIM.size(); ++j) {
        graphsSIM[j]->SetMarkerColor(colors[j]);
        graphsSIM[j]->SetMarkerStyle(markers[j]);
        graphsSIM[j]->SetMarkerSize(1.7);
        graphsSIM[j]->SetLineColor(colors[j]);
        graphsSIM[j]->GetYaxis()->SetLimits(0,1.2);
        graphsSIM[j]->SetName(Form("%i",j));
        legend -> AddEntry(graphsSIM[j], names[j], "p");
        mg->Add(graphsSIM[j]);
    }
    mg->Draw("ap");
    TPaveText *text5 = new TPaveText(0.232,0.86,0.30,0.88,"brNDC");
    text5->SetTextSize(0.035);
    text5->SetLineColor(0);
    text5->SetShadowColor(0);
    text5->SetFillColor(0);
    text5->SetTextFont(42);
    text5->AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
//    mg->GetYaxis()->SetRangeUser(0.,20.);
    text5->Draw("same");
    legend->Draw("same");

    fOut->Close();

}

//_______________________________________________________________________
void BDTCutEstimate() {
    Double_t ptMin[]={1, 2};
    Double_t ptMax[]={2, 3};
    Double_t nTrees[]={300, 300};
    Double_t depth[]={3, 2};
    const int nBins=sizeof(ptMin)/ sizeof(Double_t);

    Double_t bdtResponse[nBins];

    gROOT->ProcessLine(".L analyse/FitD0Peak.cpp++");

    for (int i = 0; i < nBins; ++i) {
        bdtResponse[i] = doSignificance(ptMin[i], ptMax[i], nTrees[i], depth[i]); //bdt response from previous measurements
//    mCuts.push_back(Form("BDTresponse>=%.3f", bdtResponse[i]));
        mCuts.push_back(Form("D_pt>=%.3f && D_pt<%.3f", ptMin[i], ptMax[i]));
        mPidCut=Form("BDTresponse>=%.3f", bdtResponse[i]);

        doRawYield(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]); //raw yield with ideal BDT response and given bdt training - take out_local.root from the correct folder, make cuts and plot it
        doRawYieldSIM(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]); //raw yield with ideal BDT response and given bdt training - take out_local.root from the correct folder, make cuts and plot it
//    doEfficiency(2, 3, 300, 2, bdtResponse[i]);
        correctYield(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]);
        mCuts.clear();
        mCuts.shrink_to_fit();
    }

    plotTogether(nBins, ptMin, ptMax, nTrees, depth, bdtResponse);

}



