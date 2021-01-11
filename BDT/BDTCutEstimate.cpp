///////////////////////////////////
/*
 * To set before running on your data:
    In doSignificance():
        inputBckg - root file with results from TMVA training (generated in TMVA training)
        cEff->SaveAs(..) - where you want to save
        inputBckg -

    In estimateNsigNbckg():
            TFile *inputSim = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMC.0910.fullEff.root"); // simulated sig input
            TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.0110.root"; - real data background (+ after that,  getting NTUple)



    Main function is BDTCutEstimate() at the end of this file.
    There you set:
        D meson pt bins
        BDT setups
        bdtResponse[] - cuts on BDT



*/
/////////////////////////////////////
#include "tmvaCuts.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
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
#include "TGaxis.h"
#include "TDatime.h"
#include <iostream>
#include <iostream>
#include <fstream>
#include "/home/lukas/work/tmva_d0/BDT/analyse/FitD0Peak.hh"

using namespace std;
Double_t NSig = -1;
Long64_t Nbckg = -1;

void estimateNsigNbckg(Double_t, Double_t);
Double_t doSignificance(double, double, int, int);
void doRawYieldSim(Double_t, Double_t, int, int, Double_t);
void doRawYield(Double_t, Double_t, int, int, Double_t);
void project_bdt_oneCut(Double_t, Double_t, Double_t, Double_t, Double_t, TString, TString, TString);
void project_bdt_oneCut_SIM(Double_t, Double_t, Double_t, Double_t, Double_t, TString, TString,TString);
void correctYield(Double_t, Double_t, Double_t, Double_t, Double_t);
void correctYieldTopoOnly(Double_t, Double_t, Double_t, Double_t, Double_t);
void plotTogether(Int_t, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *);
void plotTogetherSIM(Int_t, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *);
void drawSignificanceData(Int_t, Double_t *, Double_t *, Double_t *, Double_t *, Double_t *);
float ratioError(float, float);

void setCurrentFolder(TString);

TGraphErrors* makeGraphOnePoint(float, float, float, float);


TString folderDate;
TString simulationFileName = "out_local_SIM.root";
TString dataFileName = "out_local.root";

std::vector<TCut> mCuts;
TCut mBDTCut;
TCut precutsTMVA;

Double_t doSignificance(double ptmin, double ptmax, int nTrees, int maxDepth) {
    //root file resulting from TMVA training - with all of the information about training....
    TFile *inputBckg = new TFile(Form("pt_%.0f_%.0f/n%i_d%i/TMVA_bdt_d0_pt_%.1f_%.1f.root", ptmin, ptmax, nTrees, maxDepth, ptmin, ptmax));

    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.03);

    Float_t maxYaxis=1.2;

    TH1D *hS = (TH1D*) inputBckg->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffS"); //signal efficiency histo
    hS->SetLineColor(4);
    hS->SetStats(0);
    hS->GetYaxis()->SetTitle("Efficiency");
    hS->GetYaxis()->SetTitleOffset(0.93);
    hS->GetXaxis()->SetTitle("BDT Response Cut");
    hS->GetXaxis()->SetLabelFont(62);
    hS->GetXaxis()->SetTitleFont(62);
    hS->GetYaxis()->SetLabelFont(62);
    hS->GetYaxis()->SetTitleFont(62);
    hS->GetXaxis()->SetTitleSize(0.04);
    hS->GetYaxis()->SetTitleSize(0.04);
    hS->GetYaxis()->SetRangeUser(0, maxYaxis);
    hS->SetTitle("");
    hS->SetLineWidth(3);

    TH1D *hB = (TH1D*) inputBckg->Get("dataset/Method_BDT/BDT/MVA_BDT_trainingEffB"); //background efficiency histo
    hB->SetLineColor(2);
    hB->SetStats(0);
    hB->SetTitle("");
    hB->SetLineWidth(3);

    //////////////////////////////////////////////////////
    TCanvas *cEff = new TCanvas("cEff", "cEff", 1100, 800);
//    cEff->SetGrid();
    hS->Draw();
    hB->Draw("same");
    tx2.DrawLatex(0.1,0.92,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));
    cEff->Update();
    TLegend *legend = new TLegend(0.115665, 0.78, 0.26, 0.894);
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(hS, "Signal", "pl");
    legend->AddEntry(hB, "Background", "pl");
    legend->Draw("same");

    //////////////////////////////////////////////////////
    TH1D *hSign = (TH1D*) hS->Clone("hSign");
    hSign->SetDirectory(0);
    hSign->SetLineColor(8);
    hSign->SetTitle("");

    TH1D *hPurity = (TH1D*) hS->Clone("hSign");
    hPurity->SetDirectory(0);
    hPurity->SetLineColor(4);
    hPurity->SetLineStyle(9);
    hPurity->SetTitle("");

    estimateNsigNbckg(ptmin, ptmax);
    if (NSig<0 || Nbckg<0) return -999;

    Double_t S, B, signCheck;
    for (int i = 1; i < hS->GetNbinsX()+1; ++i) {
        S = hS->GetBinContent(i);
        B = hB->GetBinContent(i);
        signCheck=sqrt(NSig*S + (Double_t)Nbckg*B);
        if (signCheck==0.){
            hSign->SetBinContent(i,0);
            hPurity->SetBinContent(i,0);
        }
        else  {
            hSign->SetBinContent(i, NSig*S/(sqrt(NSig*S + (Double_t)Nbckg*B)));
            hPurity->SetBinContent(i, NSig*S/(NSig*S + (Double_t)Nbckg*B));
        }

        if (hSign->GetMaximum()<hSign->GetBinContent(i)) hSign->SetMaximum(hSign->GetBinContent(i));

    }

    hPurity->Draw("same");

    Float_t rightmax = 1.1*hSign->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hSign->Scale(scale);

    // draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(8);
    axis->SetTextColor(8);
    axis->SetTextSize(0.04);
    axis->SetLabelColor(8);
    axis->SetTitle("Significance");
    axis->SetTitleOffset(1.1);
    axis->Draw("same");

    hSign->Draw("same");
    legend->AddEntry(hSign, "S/#sqrt{S+B}", "pl");

    Double_t maxSign=hSign->GetBinCenter(hSign->GetMaximumBin());
    cout<<"Maximum significance from the historgram is at: "<<maxSign<<endl;
    TF1 *signFct = new TF1("signFct", "[0]+[1]*x-[2]*(x-[3])*(x-[3])", maxSign-0.16, maxSign+0.16);
    signFct->SetParameters(5/rightmax, 0, 100/rightmax, maxSign);
    signFct->SetLineColor(28);
    hSign->Fit("signFct","RLL");
    cEff->cd();

    Double_t bdtResponse = signFct->GetMaximumX();
    if (abs(maxSign-bdtResponse)>0.05){
        hSign->Fit("signFct","LL", "", maxSign-0.1, maxSign+0.1);
        bdtResponse = signFct->GetMaximumX();
    }

    if (abs(maxSign-bdtResponse)>0.05){
        hSign->Fit("signFct","LL", "", maxSign-0.05, maxSign+0.05);
        bdtResponse = signFct->GetMaximumX();
    }

    TPaveText *text5 = new TPaveText(0.48725,0.9034,0.908,0.977,"brNDC");
    text5->SetTextSize(0.025);
    text5->SetLineColor(0);
    text5->SetShadowColor(0);
    text5->SetFillStyle(0);
    text5->SetTextAlign(11);
    text5->AddText(Form("Maximum significance at BDT reponse cut: %3.4f", bdtResponse));
    text5->AddText(Form("For number of signal (background): %.0f (%lli)", NSig, Nbckg));
    text5->Draw("same");

    TLine *leftline2 = new TLine(bdtResponse, hS->GetMinimum(), bdtResponse, 1.1*hSign->GetMaximum());
    leftline2->SetLineStyle(9);
    leftline2->SetLineWidth(2);
    leftline2->SetLineColor(28);
    leftline2->Draw("same");

    legend->Draw("same");

    cEff->SaveAs(Form("finalAnalysis/%s/img/BDT_significance_%.1f_%.1f_n%i_d%i.png", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
    cEff->SaveAs(Form("finalAnalysis/%s/img/BDT_significance_%.1f_%.1f_n%i_d%i.pdf", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
    cEff->Close();

    ////////////////////////////////////////////////////////////
    TString nameData=Form("pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f",(float)ptmin, (float)ptmax, (float)nTrees, (float)maxDepth);
    TFile *fSignificanceData = TFile::Open(Form("pt_%.0f_%.0f/n%i_d%i/analyse/significance_%s.root", ptmin, ptmax, nTrees, maxDepth, nameData.Data()));
    if (fSignificanceData) {
        TCanvas *cSignData = new TCanvas("cSignData", "cSignData", 1100, 800);
        TGraphErrors *grSignData = (TGraphErrors *) fSignificanceData->Get(Form("gr_sign_%s", nameData.Data()));
        Double_t x,y,maximum=0;
        for (int i = 1; i < grSignData->GetN(); ++i) {
            grSignData->GetPoint(i,x,y);
            if (y>maximum) maximum=y;
            if (abs(x-bdtResponse)<0.02) maximum=y;
        }
        maximum=1.2*maximum;
        grSignData->Draw("ap");
        grSignData->GetYaxis()->SetRangeUser(0,maximum);
        grSignData->GetXaxis()->SetLimits(0,1);
        grSignData->SetMarkerSize(1.2);

        leftline2->SetY1(0);
        leftline2->SetY2(maximum);
        leftline2->Draw("same");
        tx2.DrawLatex(0.1,0.92,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

        cSignData->SaveAs(Form("finalAnalysis/%s/img/root/BDT_significance_DATA_%.1f_%.1f_n%i_d%i.root", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
        cSignData->SaveAs(Form("finalAnalysis/%s/img/BDT_significance_DATA_%.1f_%.1f_n%i_d%i.pdf", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
        cSignData->SaveAs(Form("finalAnalysis/%s/img/BDT_significance_DATA_%.1f_%.1f_n%i_d%i.png", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
        cSignData->Close();
    }
    fSignificanceData->Close();
    ////////////////////////////////////////////////////////////

    ofstream myfile;
    myfile.open(Form("finalAnalysis/%s/BDT_cut_%.1f_%.1f_n%i_d%i.txt", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth));
    myfile << bdtResponse;
    myfile << "\n";
    myfile << signFct->Eval(bdtResponse);
    myfile.close();

    inputBckg->Close();
    return bdtResponse;
}

//___________________________________________________________________________________
void estimateNsigNbckg(Double_t ptmin, Double_t ptmax){
    Double_t dpT=ptmax-ptmin;
    Double_t pT=(ptmin+ptmax)/2;
    Double_t dy=2;
    Double_t nevt = 93e6;
    Double_t BR=0.0395;

    auto *fppSys = new TFile("out_ppsys.root","READ");
    TF1 *f1D0pp = (TF1*) fppSys->Get("Levynew_pp");
    cout<<f1D0pp->Eval(1.5)<<" "<<f1D0pp->Eval(2.5)<<" "<<f1D0pp->Eval(4)<<endl;
    fppSys->Close();

    Double_t massMin=1.76;
    Double_t massMax=1.96;

    if (ptmin>2.9) {
        massMin=1.;
        massMax=3.;
    }

    TCut massCut = Form("D_mass >= %1.2f && D_mass < %1.2f", massMin, massMax);
    TCut mycuts = Form("D_pt>=%1.2f && D_pt<%1.2f && k_pt>%1.2f && pi1_pt>%1.2f && "
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.1 && "
                       "pi1_dca>%f && pi1_dca<0.1 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       ptmin, ptmax, tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

    TFile *inputSim = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMC.0910.fullEff.root");
    TTree *signalTree = (TTree*) inputSim->Get("ntp_signal");
    TCut pid = "pid>0.5";

    cout<<mycuts<<endl;

    Long64_t signal=signalTree->GetEntries(mycuts+pid);
    Long64_t all=signalTree->GetEntries(Form("D_pt>=%1.2f && D_pt<%1.2f", ptmin, ptmax));

//    TF1 *reco = new TF1("reco", "0.5*([0]+[1]*x+[2]*x*x+[3]*x*x*x)", 0, 10);
// //   reco->SetParameters(0.0239696, -0.0503129, 0.0556806, -0.0041924);
//    reco->SetParameters(0.00174319, -0.00388281, 0.0048162, -0.000424064);
//    Double_t totaleff=reco->Eval(pT);

    Double_t eff = (double)signal/(double)all;
    cout<<"Estimated efficiency is: "<<eff<<endl;

//    Float_t invYield=fD0->Eval(pT);
    Float_t invYield=10*f1D0pp->Eval(pT)/50;//12: 0.4639 7: 0.4893//  10 is for nBin, 50 is just how it is saved...

    NSig = 2*TMath::Pi()*pT*dpT*dy*2*invYield*nevt*BR*eff*0.7;

    cout<<"NSIG = "<<NSig<<endl;

    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.0110.root";
    TFile *inputBckg = new TFile(input);
    TTree *bckgTree = (TTree*) inputBckg->Get("ntp_background");
    Nbckg=bckgTree->GetEntries(mycuts+massCut); //mascut helped to move sign max to the left
    cout<<"Nbckg = "<<Nbckg<<endl;

    inputSim->Close();
    inputBckg->Close();
}

//_______________________________________________________________________
void doRawYield(Double_t ptmin, Double_t ptmax, int nTrees, int maxDepth, Double_t bdtCut) {
    TString inputFile=Form("/home/lukas/work/tmva_d0/BDT/finalAnalysis/%s/out_local_pt_%.1f_%.1f_n%i_d%i.root", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth);
    gSystem->Exec(Form("cp /home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/%s %s", ptmin, ptmax, nTrees, maxDepth, dataFileName.Data(), inputFile.Data()));
    project_bdt_oneCut(ptmin, ptmax, (Double_t)nTrees, (Double_t)maxDepth, bdtCut, inputFile, "signal", "background");
}

//_______________________________________________________________________
void doRawYieldSIM(Double_t ptmin, Double_t ptmax, int nTrees, int maxDepth, Double_t bdtCut) {
    TString inputFile=Form("/home/lukas/work/tmva_d0/BDT/finalAnalysis/%s/out_local_SIM_pt_%.1f_%.1f_n%i_d%i.root", folderDate.Data(), ptmin, ptmax, nTrees, maxDepth);
    gSystem->Exec(Form("cp /home/lukas/work/tmva_d0/BDT/pt_%.0f_%.0f/n%i_d%i/%s %s", ptmin, ptmax, nTrees, maxDepth, simulationFileName.Data(), inputFile.Data()));
    project_bdt_oneCut_SIM(ptmin, ptmax, (Double_t)nTrees, (Double_t)maxDepth, bdtCut, inputFile, "ntp_signal", "");
}

//_________________________________________________________________________________________________________________________________
void project_bdt_oneCut(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange, TString input, TString signalTupleName, TString backgTupleName) {
    bool mixed=false;
    TFile *f = new TFile(Form("finalAnalysis/%s/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"RECREATE");
    TFile* data = new TFile(input ,"r");

    TH1F *his[2];
    TString name[2] = {signalTupleName, backgTupleName};
    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_"+name[0]), (TNtuple*)data -> Get("ntp_"+name[1])};

    float D_mass, D_pt;

    TCut setCuts = "";
    TCut pTCut = Form("D_pt>=%.3f && D_pt<%.3f", ptmin, ptmax);
    setCuts+=pTCut;

    for(unsigned int k = 0; k < mCuts.size(); ++k) { //for now, mCuts has probably only pT, take care about it
        setCuts += mCuts[k];
    }
    setCuts += mBDTCut;
    setCuts += precutsTMVA;
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

    TString outFile=Form("finalAnalysis/%s/signals_%.4fbdt_pt_%.1f_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax);
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

    TGraphErrors* grResultsFit = new TGraphErrors();
    grResultsFit->SetPoint(0, ptBin, fitmass->getRawYieldFit());
    grResultsFit->SetPointError(0, 0, fitmass->getRawYieldFitError());
    grResultsFit->SetMarkerColor(1);
    grResultsFit->SetMarkerStyle(21);
    grResultsFit->SetTitle("");
    grResultsFit->GetXaxis()->SetLabelSize(0.04);
    grResultsFit->GetYaxis()->SetLabelSize(0.04);
    grResultsFit->GetXaxis()->SetTitleSize(0.045);

    TCut cutMass = Form("D_mass>%f && D_mass<%f", fitmass->getMean()-3*fitmass->getSigma(), fitmass->getMean()+3*fitmass->getSigma());
    setCuts+=cutMass;

    Float_t nSignalTuple = ntp[0]->GetEntries(setCuts);
    Float_t nBackgroundTuple = ntp[1]->GetEntries(setCuts);
    nSignalTuple-=nBackgroundTuple;
//    Float_t significanceTuple = (Float_t)nSignalTuple/(sqrt((Float_t nSignalTuple))
    Float_t significanceTuple = nSignalTuple/(sqrt(nSignalTuple+2*nBackgroundTuple));
    cout<<"---------------Significance from tuple projection------------------"<<endl;
    cout<<setCuts<<endl;
    cout<<"---------------background is: "<<nBackgroundTuple<<endl;
    cout<<"---------------signal is: "<<nSignalTuple<<endl;
    cout<<"---------------significance is "<<significanceTuple<<endl;

    delete fitmass;

    f->cd();
    grResultsFit->Write("raw_yield_fit");
    grResults->Write("raw_yield");
//    grResults->Write("significance");
    listOut->Write("hists_D_mass", 1, 0);
    data->Close();
    f->Close();
}

//_________________________________________________________________________________________________________________________________
void project_bdt_oneCut_SIM(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange, TString input, TString tupleName, TString weight) {
    TFile *f = new TFile(Form("finalAnalysis/%s/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"RECREATE");

    TFile* data = new TFile(input ,"r");
    TNtuple* ntp = (TNtuple*)data -> Get(tupleName);

    TH1F *his = new TH1F(Form("D_mass_%.4f_signal", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
    TH1F *hisEmpty = new TH1F(Form("D_mass_%.4f_background", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4); //fake background for FitD0Peak
    his -> Sumw2();

    TCut setCuts;
    TCut hijingCompareCuts;

    for(unsigned int k = 0; k < mCuts.size(); ++k) { //for now, mCuts has probably only pT, take care about it
        setCuts += mCuts[k];
    }

//    TCut pTCut = Form("D_pt>=%.3f && D_pt<%.3f", ptmin, ptmax);
    TCut pTCut = Form("D_ptSIM>=%.3f && D_ptSIM<%.3f", ptmin, ptmax);
    setCuts+=pTCut;

    Long64_t nAllSim=ntp->GetEntries(setCuts*weight); //all from simu

//    setCuts+="k_pt>0.15 && pi1_pt>0.15 && etas>0";
//    Long64_t nAllSim=ntp->GetEntries(setCuts*weight); //all from simu

    setCuts+="k_pt>0.15 && pi1_pt>0.15";
    Long64_t nAcc=ntp->GetEntries(setCuts*weight);

    setCuts+="tpc>0";
    Long64_t nTpcAcc=ntp->GetEntries(setCuts*weight);

    setCuts+="hft>0";
    hijingCompareCuts=setCuts;
    Long64_t nTpcAccHft=ntp->GetEntries(setCuts*weight);

    setCuts+="pid>0";
    Long64_t nTpcAccHftPid=ntp->GetEntries(setCuts*weight);

    setCuts+=precutsTMVA;
    hijingCompareCuts+=precutsTMVA;
    Long64_t nTpcAccHftPreCuts=ntp->GetEntries(hijingCompareCuts*weight); //pid cuts
    Long64_t nTpcAccHftPidPreCuts=ntp->GetEntries(setCuts*weight); //pid cuts

    setCuts+=mBDTCut;
    hijingCompareCuts+=mBDTCut;
    Long64_t nTpcAccHftPreCutsBDT=ntp->GetEntries(hijingCompareCuts*weight); //pid cuts
    Long64_t nTpcAccHftPidPreCutsBDT=ntp->GetEntries(setCuts*weight); //pid cuts

    cout<<"Cuts set for final eff est. "<<setCuts<<endl;
    ntp->Project(his->GetName(), "D_mass", setCuts+mBDTCut); //same as in the data

    f->cd();
    TList *listOut = new TList();
    listOut->Add(his);
    listOut->Write("hists_signal", 1, 0);
    delete listOut;

    Float_t ptBin=(ptmax+ptmin)/2;
    Float_t ptBinWidth=(ptmax-ptmin)/2;

    TString outFile=Form("finalAnalysis/%s/signals_SIM_%.4fbdt_pt_%.1f_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax);

//    FitD0Peak *fitmass = new FitD0Peak(his, hisEmpty, ptmin, ptmax, outFile);
//    fitmass->doStuff();
//
//    TGraphErrors* grResults = new TGraphErrors();
//    grResults->SetPoint(0, ptBin, fitmass->getRawYield());
//    grResults->SetPointError(0, 0, fitmass->getRawYieldError());
//    grResults->SetMarkerColor(1);
//    grResults->SetMarkerStyle(21);
//    grResults->SetTitle("");
//    grResults->GetXaxis()->SetLabelSize(0.04);
//    grResults->GetYaxis()->SetLabelSize(0.04);
//    grResults->GetXaxis()->SetTitleSize(0.045);
//
//    delete fitmass;

    f->cd();
//    grResults->Write("raw_yield");

    TGraphErrors* grAcc = makeGraphOnePoint(ptBin, (float)nAcc/(float)nAllSim, ptBinWidth, ratioError((float)nAcc,(float)nAllSim) );
    grAcc->Write("grAcc");

    TGraphErrors* grTpcAcc = makeGraphOnePoint(ptBin, (float)nTpcAcc/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAcc,(float)nAllSim));
    grTpcAcc->Write("grTpcAcc");

    TGraphErrors* grTpcAccHft = makeGraphOnePoint(ptBin, (float)nTpcAccHft/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHft,(float)nAllSim));
    grTpcAccHft->Write("grTpcAccHft");

    TGraphErrors* grTpcAccHftPid = makeGraphOnePoint(ptBin, (float)nTpcAccHftPid/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHftPid,(float)nAllSim));
    grTpcAccHftPid->Write("grTpcAccHftPid");

    TGraphErrors* grTpcAccHftPidPreCuts = makeGraphOnePoint(ptBin, (float)nTpcAccHftPidPreCuts/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHftPidPreCuts,(float)nAllSim));
    grTpcAccHftPidPreCuts->Write("grTpcAccHftPidPreCuts");

    TGraphErrors* grTpcAccHftPidPreCutsBDT = makeGraphOnePoint(ptBin, (float)nTpcAccHftPidPreCutsBDT/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHftPidPreCutsBDT,(float)nAllSim));
    grTpcAccHftPidPreCutsBDT->Write("grTpcAccHftPidPreCutsBDT");

    TGraphErrors* grTpcAccHftPreCuts = makeGraphOnePoint(ptBin, (float)nTpcAccHftPreCuts/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHftPreCuts,(float)nAllSim));
    grTpcAccHftPreCuts->Write("grTpcAccHftPreCuts");

    TGraphErrors* grTpcAccHftPreCutsBDT = makeGraphOnePoint(ptBin, (float)nTpcAccHftPreCutsBDT/(float)nAllSim, ptBinWidth, ratioError((float)nTpcAccHftPreCutsBDT,(float)nAllSim));
    grTpcAccHftPreCutsBDT->Write("grTpcAccHftPreCutsBDT");

//    grResults->Write("significance");
    data->Close();
    f->Close();
}

//_______________________________________________________________________
float ratioError(float ratioUp, float ratioDown){
    float result = ratioUp*ratioUp*ratioDown + ratioDown*ratioDown*ratioUp;
    result=sqrt(result);
    result=result/(ratioDown*ratioDown);
    return result;
}

//_______________________________________________________________________
TGraphErrors* makeGraphOnePoint(float x, float y, float xE, float yE){
    TGraphErrors* gr = new TGraphErrors();
    gr->SetPoint(0, x, y);
    gr->SetPointError(0, xE, yE);
    return gr;
}

//_______________________________________________________________________
void correctYield(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange) {
    TFile *fData = new TFile(Form("finalAnalysis/%s/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"UPDATE");
    TFile *fSIM = new TFile(Form("finalAnalysis/%s/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"READ");

    Double_t rawYieldErr, rawYield, pT, effError, eff;

    TGraphErrors *grRawYield = (TGraphErrors*) fData->Get("raw_yield");
    grRawYield->GetPoint(0, pT, rawYield);
    rawYieldErr = grRawYield->GetErrorY(0);

    TGraphErrors *grEff = (TGraphErrors*) fSIM->Get("grTpcAccHftPidPreCutsBDT");
    grEff->GetPoint(0, pT, eff);
    effError = grEff->GetErrorY(0);

    Double_t dpT=ptmax-ptmin;
//    pT=(ptmin+ptmax)/2;
    Double_t dy=2;
    Double_t nevt = 90e6;
    Double_t BR=0.0395;  //(3.950+-0.031)%
    Double_t BRerror = 0.00031;
    Double_t cons = 2*TMath::Pi()*pT*dpT*dy*2*nevt;

    Double_t invYield = rawYield/(cons*BR*eff);
    Double_t invYieldError = rawYieldErr/(cons*BR*eff);// + pow(rawYield*effError/(cons*eff*eff), 2);
//    invYieldError = sqrt(invYieldError);

    TGraphErrors* grResults = new TGraphErrors();
    grResults->SetPoint(0, pT, invYield);
    grResults->SetPointError(0, dpT, invYieldError);
    grResults->SetMarkerColor(1);
    grResults->SetMarkerStyle(21);
    grResults->SetTitle("");
    grResults->GetXaxis()->SetLabelSize(0.04);
    grResults->GetYaxis()->SetLabelSize(0.04);
    grResults->GetXaxis()->SetTitleSize(0.045);

    auto *fppSys = new TFile("out_ppsys.root","READ");
    TF1 *f1D0pp = (TF1*) fppSys->Get("Levynew_pp");
    cout<<f1D0pp->Eval(1.5)<<" "<<f1D0pp->Eval(2.5)<<" "<<f1D0pp->Eval(4)<<endl;
    fppSys->Close();

    Float_t invYieldPP=f1D0pp->Eval(pT)/50;//12: 0.4639 7: 0.4893//  10 is for nBin, 50 is just how it is saved...
    Float_t raa=invYield/invYieldPP;

    TGraphErrors* grResultsRAA = new TGraphErrors();
    grResultsRAA->SetPoint(0, pT, raa);
    grResultsRAA->SetPointError(0, dpT, invYieldError*raa/invYield);
    grResultsRAA->SetMarkerColor(1);
    grResultsRAA->SetMarkerStyle(21);
    grResultsRAA->SetTitle("");
    grResultsRAA->GetXaxis()->SetLabelSize(0.04);
    grResultsRAA->GetYaxis()->SetLabelSize(0.04);
    grResultsRAA->GetXaxis()->SetTitleSize(0.045);

    fData->cd();
    grResults->Write("grInvYield", 1, 0);
    grResultsRAA->Write("grResultsRAA", 1, 0);

    cout<<"Raw yield:"<<endl;
    cout<<rawYield<<" "<<rawYieldErr<<endl;
    cout<<"Invariant yield:"<<endl;
    cout<<invYield<<" "<<invYieldError<<endl;

    fData->Close();
    fSIM->Close();
}

//_______________________________________________________________________
void correctYieldTopoOnly(Double_t ptmin, Double_t ptmax, Double_t nTrees, Double_t maxDepth, Double_t bdtRange) {
    TFile *fData = new TFile(Form("finalAnalysis/%s/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"UPDATE");
    TFile *fSIM = new TFile(Form("finalAnalysis/%s/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtRange, ptmin, ptmax, nTrees ,maxDepth),"READ");

    Double_t rawYieldErr, rawYield, pT, effError, eff;
    Double_t dpT=ptmax-ptmin;

    TGraphErrors *grRawYield = (TGraphErrors*) fData->Get("raw_yield");
    grRawYield->GetPoint(0, pT, rawYield);
    rawYieldErr = grRawYield->GetErrorY(0);

    TGraphErrors *grEff = (TGraphErrors*) fSIM->Get("grPreCutsBDT");
    grEff->GetPoint(0, pT, eff);
    effError = grEff->GetErrorY(0);

    Double_t invYield = rawYield/eff;
    Double_t invYieldError = rawYieldErr/eff;// + pow(rawYield*effError/(cons*eff*eff), 2);

    TGraphErrors* grResults = new TGraphErrors();
    grResults->SetPoint(0, pT, invYield);
    grResults->SetPointError(0, dpT/2, invYieldError);
    grResults->SetMarkerColor(1);
    grResults->SetMarkerStyle(21);
    grResults->SetTitle("");
    grResults->GetXaxis()->SetLabelSize(0.04);
    grResults->GetYaxis()->SetLabelSize(0.04);
    grResults->GetXaxis()->SetTitleSize(0.045);

    fData->cd();
    grResults->Write("grRawYieldTopoCorrected", 1, 0);

    fData->Close();
    fSIM->Close();
}


//_______________________________________________________________________
void plotTogether(Int_t nBins, Double_t *ptmin, Double_t *ptmax, Double_t *nTrees, Double_t *maxDepth, Double_t *bdtResponse) {
    TGraphErrors *gr = new TGraphErrors();
    Double_t rawYields[nBins], rawYieldsErr[nBins];
    Double_t pT[nBins], ptWidths[nBins], y[nBins], yE[nBins];

    TFile *fOut = new TFile(Form("finalAnalysis/%s/final_result.root", folderDate.Data()), "RECREATE");
    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

//    TString names1[] = {"grResultsRAA","grInvYield", "raw_yield", "grTpcAcc", "grTpcAccHft", "grTpcAccHftPid", "grTpcAccHftPidPreCuts", "grTpcAccHftPidPreCutsBDT"};
//    TString axisName[] = {"R_{AA} = dAu/pp","Invariant yield", "Raw yield", "grTpcAcc", "grTpcAccHft", "grTpcAccHftPid", "grTpcAccHftPidPreCuts", "grTpcAccHftPidPreCutsBDT"};

//    TString names1[] = {"grRawYieldTopoCorrected"};
//    TString axisName[] = {"Raw yield/#varepsilon_{BDT}"};

//    TString names1[] = {"grInvYield", "grResultsRAA", "raw_yield"};
    TString names1[] = {"grResultsRAA", "grInvYield", "raw_yield", "raw_yield_fit"};
    TString axisName[] = {"R_{AA} = dAu/pp", "Invariant yield", "Raw yield"};

    const int ngraphsData = sizeof(names1) / sizeof(TString);

    for (int k = 0; k < ngraphsData; ++k) {
        cout << "plotting: " << names1[k] << endl;
        for (int i = 0; i < nBins; ++i) {
            TFile *fData = new TFile(Form("finalAnalysis/%s/D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtResponse[i], ptmin[i], ptmax[i], nTrees[i], maxDepth[i]), "READ");
            gr = (TGraphErrors *) fData->Get(names1[k]);
            gr->GetPoint(0, pT[i], rawYields[i]);
            cout << rawYields[i] << endl;
            ptWidths[i] = gr->GetErrorX(0);
            rawYieldsErr[i] = gr->GetErrorY(0);
            fData->Close();
        }

        TGraphErrors *grYields = new TGraphErrors(nBins, pT, rawYields, 0, rawYieldsErr);
        grYields->SetMarkerColor(1);
        grYields->SetMarkerStyle(21);
        grYields->SetTitle("");
        grYields->GetXaxis()->SetLabelSize(0.04);
        grYields->GetYaxis()->SetLabelSize(0.04);
        grYields->GetXaxis()->SetTitleSize(0.045);
        grYields->GetYaxis()->SetTitle(axisName[k]);
        grYields->GetYaxis()->SetTitleOffset(1.6);
        grYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fOut->cd();
        grYields->Write(names1[k] + "_all");
    }

    fOut->Close();
    return;
}

//_______________________________________________________________________
void plotTogetherSIM(Int_t nBins, Double_t *ptmin, Double_t *ptmax, Double_t *nTrees, Double_t *maxDepth, Double_t *bdtResponse) {
    TFile *fOut = new TFile(Form("finalAnalysis/%s/final_result_SIM.root", folderDate.Data()), "RECREATE");
    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    TGraphErrors *gr = new TGraphErrors();
    Double_t pT[nBins], ptWidths[nBins], y[nBins], yE[nBins];

    std::vector<TGraphErrors*> graphsSIM;

    TString names[] = {"grAcc",
                       "grTpcAcc",
                       "grTpcAccHft",
                       "grTpcAccHftPid",
                       "grTpcAccHftPidPreCuts",
                       "grTpcAccHftPidPreCutsBDT",
                       "grTpcAccHftPreCuts",
                       "grTpcAccHftPreCutsBDT"};
    TString legendNames[] = {"Acc",
                             "TPC #times Acc",
                             "TPC #times Acc #times HFT",
                             "TPC #times Acc #times HFT #times PID",
                             "TPC #times Acc #times HFT #times PID #times Topo pre-cuts",
                             "TPC #times Acc #times HFT #times PID #times Topo pre-cuts #times BDT",
                             "TPC #times Acc #times HFT #times Topo pre-cuts",
                             "TPC #times Acc #times HFT #times Topo pre-cuts #times BDT"};

    const int ngraphs=sizeof(names)/ sizeof(TString);

    for (int j = 0; j < ngraphs; ++j) {
        for (int i = 0; i < nBins; ++i) {
            TFile *fData = new TFile(Form("finalAnalysis/%s/D0_SIM_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", folderDate.Data(), bdtResponse[i], ptmin[i], ptmax[i], nTrees[i], maxDepth[i]), "READ");
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
    out->SetLogy();
    TMultiGraph *mg = new TMultiGraph();
//    mg->GetYaxis()->SetRangeUser(0.000001,2);
    mg->SetMaximum(5.);
    mg->GetYaxis()->SetTitle("#varepsilon");
    mg->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
    mg->GetXaxis()->SetTitleSize(0.96);
    mg->GetYaxis()->SetLimits(0, 1);
    mg->SetTitle("");
    TLegend *legend = new TLegend(0.12354, 0.7497, 0.275, 0.8898);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.03);
    for (unsigned short j = 1; j < graphsSIM.size()-2; ++j) {
        graphsSIM[j]->SetMarkerColor(colors[j]);
        graphsSIM[j]->SetMarkerStyle(2);
        graphsSIM[j]->SetMarkerSize(1.7);
        graphsSIM[j]->SetLineColor(colors[j]);
        graphsSIM[j]->SetLineWidth(2);
        graphsSIM[j]->GetYaxis()->SetLimits(0,1.2);
        graphsSIM[j]->SetName(Form("%i",j));
        legend -> AddEntry(graphsSIM[j], legendNames[j], "pl");
        mg->Add(graphsSIM[j]);
    }
    mg->Draw("ap");
    TPaveText *text5 = new TPaveText(0.724,0.925,0.793,0.945,"brNDC");
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

    return;
}

//_______________________________________________________________________
void drawSignificanceData(Int_t nBins, Double_t *ptmin, Double_t *ptmax, Double_t *nTrees, Double_t *maxDepth, Double_t *bdtResponse) {
    TString nameData;
    TCanvas *cSignData = new TCanvas("cSignData", "cSignData", 1100, 800);
    Double_t maximum=0;
    Double_t maxSign[nBins];
    TGraphErrors *grSignData[nBins];
    TLine *line[nBins];
    Int_t colors[] = {8, 46, 9, 9, 40, 41, 42, 28, 2};
    Int_t markers[] = {21, 20, 34, 9, 40, 41, 42, 28, 2};

    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);

    for (int i = 0; i < nBins; ++i) {
        nameData = Form("pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f", (float) ptmin[i], (float) ptmax[i], (float) nTrees[i], (float) maxDepth[i]);
        TFile *fSignificanceData = TFile::Open(Form("pt_%.0f_%.0f/n%i_d%i/analyse/significance_%s.root", ptmin[i], ptmax[i], (int)nTrees[i], (int)maxDepth[i], nameData.Data()));

        if (fSignificanceData) {
            grSignData[i] = (TGraphErrors *) fSignificanceData->Get(Form("gr_sign_%s", nameData.Data()));
            grSignData[i]->SetMarkerColor(colors[i]);
            grSignData[i]->SetMarkerStyle(markers[i]);
            Double_t x, y;
            for (int k = 1; k < grSignData[i]->GetN(); ++k) {
                grSignData[i]->GetPoint(k, x, y);
                if (abs(x - bdtResponse[i]) < 0.0075 && y > maximum) maximum = y;
                if (abs(x - bdtResponse[i]) < 0.0075) maxSign[i]=y;
            }
        }
        fSignificanceData->Close();
    }

    grSignData[0]->GetYaxis()->SetLabelSize(0.045);
    grSignData[0]->GetXaxis()->SetLabelSize(0.045);
    grSignData[0]->GetXaxis()->SetTitleSize(0.055);
    grSignData[0]->GetYaxis()->SetTitleSize(0.055);
    grSignData[0]->GetYaxis()->SetTitleOffset(0.8);
    grSignData[0]->GetXaxis()->CenterTitle();
    grSignData[0]->GetYaxis()->CenterTitle();

    TLegend *legend = new TLegend(0.126, 0.71, 0.277, 0.85);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);

    maximum=1.2*maximum;
    for (int j = 0; j < nBins; ++j) {
        grSignData[j]->GetYaxis()->SetRangeUser(0,9);
        grSignData[j]->GetXaxis()->SetLimits(0.08,0.88);
        grSignData[j]->SetMarkerSize(1.2);

        legend->AddEntry(grSignData[j], Form("%0.f<p_{T}<%0.f GeV/c", ptmin[j], ptmax[j]), "p");

        if (j==0) grSignData[j]->Draw("ap");
        else grSignData[j]->Draw("p same");
        line[j] = new TLine();
        line[j]->SetLineStyle(2);
        line[j]->SetLineWidth(2);
        line[j]->SetLineColor(colors[j]);
        line[j]->SetX1(bdtResponse[j]);
        line[j]->SetX2(bdtResponse[j]);
        line[j]->SetY1(0);
        line[j]->SetY2(9);
//        line[j]->SetY2(maxSign[j]);
        line[j]->Draw("same");
    }

    TPaveText *text5 = new TPaveText(0.232,0.86,0.30,0.88,"brNDC");
    text5->SetTextSize(0.035);
    text5->SetLineColor(0);
    text5->SetShadowColor(0);
    text5->SetFillColor(0);
    text5->SetTextFont(42);
    text5->AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
    text5->Draw("same");
    legend->Draw("same");

}
//_______________________________________________________________________
void setCurrentFolder(TString name) {
    folderDate=name;
}
//_______________________________________________________________________
void BDTCutEstimate() {
//    Double_t ptMin[]={1, 2, 3};
//    Double_t ptMax[]={2, 3, 5};
//

    folderDate="";

    if(folderDate=="") {
        auto *time = new TDatime();
        Int_t day = time->GetDay();
        Int_t month = time->GetMonth();
        Int_t hour = time->GetHour();
        Int_t minute = time->GetMinute();
        folderDate = Form("%i%i_%i%i", month, day, hour, minute);
        cout << folderDate << endl;
        gSystem->Exec(Form("mkdir -p finalAnalysis/%s", folderDate.Data()));
        gSystem->Exec(Form("mkdir finalAnalysis/%s/img", folderDate.Data()));
        gSystem->Exec(Form("mkdir finalAnalysis/%s/img/root", folderDate.Data()));
    }

    //GOOD SET
//    Double_t ptMin[]={1,2,3};
//    Double_t ptMax[]={2,3,5};
//
//    Double_t nTrees[]={100,150,400};
//    Double_t depth[]={3,3,3};
//
//    Double_t bdtResponse[]={0.7552, 0.64516, 0.53154};


    Double_t ptMin[16];
    Double_t ptMax[16];
    Double_t nTrees[16];
    Double_t depth[16];
    Double_t bdtResponse[16];

    for (int j = 0; j < 16; ++j) {
        ptMin[j] = 1+j*0.25;
        ptMax[j] = 1+(j+1)*0.25;
    }
    const int nBins=sizeof(ptMin)/ sizeof(Double_t);

    TString inputSim[nBins];

    for (int k = 0; k < nBins; ++k) {
        if (ptMin[k]>=1 && ptMax[k]<=2){
            nTrees[k]=100;
            depth[k]=3;
            bdtResponse[k]=0.7552;
            inputSim[k]="/home/lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root";
        }

        if (ptMin[k]>=2 && ptMax[k]<=3){
            nTrees[k]=150;
            depth[k]=3;
            bdtResponse[k]=0.64516;
            inputSim[k]="/home/lukas/work/tmva_d0/BDT/pt_2_3/n150_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root";

        }

        if (ptMin[k]>=3 && ptMax[k]<=5){
            nTrees[k]=400;
            depth[k]=3;
            bdtResponse[k]=0.53154;
            inputSim[k]="/home/lukas/work/tmva_d0/BDT/pt_3_5/n400_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root";

        }
    }

    for (int l = 0; l < nBins; ++l) {
        cout<<ptMin[l]<<" "<<ptMax[l]<<" "<<nTrees[l]<<" "<<depth[l]<<" "<<bdtResponse[l]<<endl;
    }

    gROOT->ProcessLine(".L analyse/FitD0Peak.cpp++");

    precutsTMVA = Form("k_pt>%1.2f && pi1_pt>%1.2f && "
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
//    TString input = "/hom//lukas/work/tmva_d0/BDT/pt_1_2/n100_d3/out_local_SIM_ntpTMVA_full_D0.toyMc.0303.root";
//    TString input = "/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMC.0910.fullEff.root";

    for (int i = 0; i < nBins; ++i) {
////    for (int i = 0; i < 1; ++i) {
//        bdtResponse[i] = doSignificance(ptMin[i], ptMax[i], nTrees[i], depth[i]); //bdt response from previous measurements
//        mCuts.push_back(Form("D_pt>=%.3f && D_pt<%.3f", ptMin[i], ptMax[i]));
//        mCuts.push_back("refMult>10");
        mBDTCut=Form("BDTresponse>=%.3f", bdtResponse[i]);

//        project_bdt_oneCut_SIM(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i], inputSim[i], "signal", "background");
//

//        doRawYield(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]); //raw yield with ideal BDT response and given bdt training - take out_local.root from the correct folder, make cuts and plot it
//        doRawYieldSIM(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]); //raw yield with ideal BDT response and given bdt training - take out_local.root from the correct folder, make cuts and plot it
//        correctYield(ptMin[i], ptMax[i], nTrees[i], depth[i], bdtResponse[i]);
        mCuts.clear();
        mCuts.shrink_to_fit();
    }

    plotTogetherSIM(nBins, ptMin, ptMax, nTrees, depth, bdtResponse);
//    drawSignificanceData(nBins, ptMin, ptMax, nTrees, depth, bdtResponse);
    cout<<"Work is done. Everything is in folder "<<folderDate<<endl;

}



