#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TEventList.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "FitD0Peak.hh"

using namespace std;
using namespace TMath;
ClassImp(FitD0Peak)

FitD0Peak::FitD0Peak(TH1F* sigInput, TH1F* bckgInput, Float_t ptMinInput, Float_t ptMaxInput, TString mOutputFileName) : TObject(),
                   mSigma(999), mSigmaE(999), mMean(999), mMeanE(999), mHeight(999), mHeightE(999),
                   significanceBins(999), SoverB(999), rawYieldError(999), rawYield(999),
                   mFitRMin(1.7), mFitRMax(2.0),
                   peakMin(1.82), peakMax(1.88),
                   nsigma(3), lines(false) {
    sigOrig=(TH1F*)sigInput->Clone(sigInput->GetName());
    bckgOrig=(TH1F*)bckgInput->Clone(bckgInput->GetName());
    bckgToWork=(TH1F*)bckgOrig->Clone(bckgOrig->GetName());

    isMxdEv=false;
    scale=false;

    ptMin=ptMinInput;
    ptMax=ptMaxInput;

    sigOrig->Rebin(10);
    bckgOrig->Rebin(10);
    bckgToWork->Rebin(10);

    binSize=sigOrig->GetBinWidth(sigOrig->FindBin(1.864));

    fOut = new TFile(mOutputFileName,"recreate");
    text1 = new TPaveText(0.719,0.9229,0.998,0.980,"brNDC");
    text1->SetTextSize(0.04);
    text1->SetShadowColor(0);
    text1->SetLineColor(0);
    text1->SetFillColor(0);
    text1->AddText("d+Au 200 GeV");

    tx2.SetNDC();
    tx2.SetTextSize(0.04);

    setHistoStyle(sigOrig, 46, 20);
    setHistoStyle(bckgOrig, 1, 25);
    setHistoStyle(bckgToWork, 2, 25);

    sigSubtracted=(TH1F*)sigOrig->Clone(sigOrig->GetName());
    setHistoStyle(sigSubtracted, 9, 20);

    mTupleSignal=new TNtuple();
    mTupleBackground=new TNtuple();

//    doStuff();


}

//____________________________________________________________________________________________________________________________
void FitD0Peak::doStuff(){
        Float_t intLowLow = 1.7;
        Float_t intLowUp = 1.8;
        Float_t intUpLow = 1.92;
        Float_t intUpUp = 2.;

        if (isMxdEv) bckgToWork->Add(mxdOrig, 1);

        Float_t hBackIntegral = bckgToWork -> Integral(bckgToWork -> FindBin(intLowLow), bckgToWork->FindBin(intLowUp), "") + bckgToWork -> Integral(bckgToWork->FindBin(intUpLow), bckgToWork->FindBin(intUpUp), ""); //integrals of original histos
        Float_t hSignIntegral = sigOrig -> Integral(sigOrig -> FindBin(intLowLow), sigOrig -> FindBin(intLowUp), "") + sigOrig -> Integral(sigOrig->FindBin(intUpLow), sigOrig->FindBin(intUpUp), ""); //integrals of original histos

        std::cout<<"Backgroung Integral: "<<hBackIntegral<<endl;
        std::cout<<"Signal Integral :"<<hSignIntegral<<endl;

        if(scale) bckgToWork->Scale(hSignIntegral/hBackIntegral);
        sigSubtracted->Add(bckgToWork,-1);


        Double_t nentriesSig = sigOrig->Integral(sigOrig->FindBin(1.7), sigOrig->FindBin(2),"");
        fitComeOn();
        fitFunction();

        plotPub();
        plotOrig();
        plotPubWithResidual();
//        ofstream yieldsF;
//        yieldsF.open("raw_yields_fct.txt", std::ios::app);
//        yieldsF<<ptMin<<" "<<ptMax<<" "<<integral_function_yield/binSize<<" "<<integral_function_yield_error/binSize<<endl;
//        yieldsF.close();

//        ofstream yields;
//        yields.open("raw_yields.txt", std::ios::app);
//        yields<<ptMin<<" "<<ptMax<<" "<<rawYield<<" "<<rawYieldError<<endl;
//        yields.close();
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::setHistoStyle(TH1F* histo, Int_t color, Int_t marker) {
    histo->Sumw2();

    histo->SetTitle("");
    histo->SetMarkerStyle(marker);
    histo->SetLineColor(color);
    histo->SetMarkerColor(color);
    histo->SetMarkerSize(2);

    histo->GetYaxis()->SetTitle("Counts");
    histo->GetYaxis()->SetTitleOffset(1.);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetTitleFont(42);
    histo->GetYaxis()->CenterTitle(kTRUE);

    histo->GetXaxis()->SetTitle("Invariant mass, m_{K#pi} [GeV/c^{2}]");
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetXaxis()->SetTitleFont(42);
    histo->GetXaxis()->SetRangeUser(0.65,2.4);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->CenterTitle(kTRUE);

    histo->GetXaxis()->SetRangeUser(mFitRMin,mFitRMax);
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::plotPub(){
    TCanvas *cP = new TCanvas("cPublish","cPublish",1200,1000);
    gStyle->SetOptStat(0);
    gPad->SetMargin(0.1,0.05,0.13,0.08);

    sigSubtracted->Draw();
//    fun0->Draw("same");
    gStyle->SetOptFit(111);

    TPaveStats *st = (TPaveStats*)sigSubtracted->FindObject("stats");
    st->SetX1NDC(0.7); //new x start position
    st->SetY1NDC(0.7); //new x end position
    st->SetX2NDC(0.949833); //new x start position
    st->SetY2NDC(0.920598); //new x end position

    TPaveText *textPub1 = new TPaveText(0.13,0.717,0.205,0.823,"brNDC");
    textPub1->SetTextSize(0.04);
    textPub1->SetLineColor(0);
    textPub1->SetShadowColor(0);
    textPub1->SetFillColor(0);
    textPub1->SetTextFont(42);
    textPub1->AddText(Form("Significance: %0.3g, S/B: %0.3g", significanceBins, SoverB));
    textPub1->AddText(Form("Raw yield: %0.2g #pm %0.2g", rawYield, rawYieldError));
    textPub1->SetTextAlign(12);
    textPub1->Draw("same");

    TLegend *legendPub = new TLegend(0.127,0.817, 0.287, 0.907,"","brNDC");
    legendPub->AddEntry(sigSubtracted, "US - LS", "pl");
//    legendPub->AddEntry(fun0, "Gaussian + linear fit", "pl");
    legendPub->AddEntry(sigSubtracted->GetFunction("fun0"), "Gaussian + linear fit", "pl");
    legendPub->SetFillStyle(0);
    legendPub->SetLineColor(0);
    legendPub->SetTextSize(0.04);
    legendPub->Draw("same");

    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptMin, ptMax));
    text1->Draw("same");

    fOut->cd();
    cP->Write(Form("%.1f_%.1f_publ",ptMin, ptMax));
    cP->Close();
    delete cP;
    delete textPub1;
    delete legendPub;
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::plotPubWithResidual(){
    TCanvas *cP = new TCanvas("cPublish1","cPublish1",1200,1000);
    gPad->SetMargin(0.1,0.05,0.13,0.08);
    gStyle->SetOptStat(0);

    sigSubtractedResidualBckg->SetStats(0);
    sigSubtractedResidualBckg->Draw();

    sigSubtracted->Draw("same");
//    gStyle->SetOptFit(111);

    TPaveStats *st = (TPaveStats*)sigSubtracted->FindObject("stats");
    st->SetOptFit(111);
    st->SetX1NDC(0.7); //new x start position
    st->SetY1NDC(0.7); //new x end position
    st->SetX2NDC(0.949833); //new x start position
    st->SetY2NDC(0.920598); //new x end position

    TPaveText *textPub1 = new TPaveText(0.13,0.63,0.20,0.737,"brNDC");
    textPub1->SetTextSize(0.03);
    textPub1->SetLineColor(0);
    textPub1->SetShadowColor(0);
    textPub1->SetFillColor(0);
    textPub1->SetTextFont(42);
    textPub1->AddText(Form("Significance: %0.3g, S/B: %0.3g", significanceBins, SoverB));
    textPub1->AddText(Form("Raw yield: %0.2g #pm %0.2g", rawYield, rawYieldError));
    textPub1->SetTextAlign(12);
    textPub1->Draw("same");

    TF1 *resGaus = new TF1("resGaus","gaus(0)",mFitRMin,mFitRMax);
    resGaus->SetParameters(fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resGaus->SetLineColor(15);
//    resGaus->SetLineStyle(9);
    resGaus->Draw("same");

    TF1 *resLinear = new TF1("resLinear","pol1",mFitRMin,mFitRMax);
    resLinear->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
    resLinear->SetLineStyle(3);
    resLinear->SetLineWidth(4);
    resLinear->SetLineColor(1);
    resLinear->Draw("same");


    TLegend *legendPub = new TLegend(0.127,0.753, 0.287, 0.907,"","brNDC");
    legendPub->AddEntry(sigSubtracted, "US - LS", "pl");
    legendPub->AddEntry(sigSubtractedResidualBckg, "US - LS - residual", "pl");
    legendPub->AddEntry(sigSubtracted->GetFunction("fun0"), "Gaussian + linear fit", "pl");
    legendPub->AddEntry(resLinear, "Residual background (linear)", "pl");
    legendPub->AddEntry(resGaus, "Gaussian peak", "pl");
    legendPub->SetFillStyle(0);
    legendPub->SetLineColor(0);
    legendPub->SetTextSize(0.03);
    legendPub->Draw("same");

    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptMin, ptMax));
    text1->Draw("same");

    fOut->cd();
    cP->Write(Form("%.1f_%.1f_publ_peak",ptMin, ptMax));
    cP->Close();
    delete cP;
    delete textPub1;
    delete legendPub;
    delete resGaus;
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::plotOrig(){
    TCanvas *c5 = new TCanvas("c5","c5",1200,1000);
    gPad->SetMargin(0.1,0.05,0.13,0.08);
    gStyle->SetOptStat(0);

    //    FINAL PUBLICATION QM18
//    TPaveText *text5 = new TPaveText(0.74,0.71,0.81,0.94,"brNDC");
//    text5->SetTextSize(0.05);
//    text5->SetLineColor(0);
//    text5->SetShadowColor(0);
//    text5->SetFillColor(0);
//    text5->SetTextFont(42);
//    text5->AddText("STAR Preliminary");
//    text5->AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
//    text5->AddText(Form("%3.1f < p_{T} < %3.1f GeV/c", ptMin, ptMax));
//    TText *t = text5->GetLineWith("STAR");
//    t->SetTextColor(kRed);

    sigOrig->SetStats(0);
    bckgOrig->SetStats(0);
    if(isMxdEv) mxdOrig->SetStats(0);


    sigOrig->Draw("");
    bckgOrig->Draw("same");
    sigOrig->Draw("same");
    if(isMxdEv) mxdOrig->Draw("same");
    if(isMxdEv) bckgToWork->Draw("same");

    if (lines) {
        leftLine->Draw("same");
        rightLine->Draw("same");
    }

    TLegend *legend3 = new TLegend(0.1137,0.812, 0.351, 0.913,"","brNDC");
    legend3->AddEntry(bckgOrig, "Like-sign (LS) background", "pl");
    legend3->AddEntry(sigOrig, "Unlike-sign (US) signal", "pl");
    legend3->SetFillStyle(0);
    legend3->SetFillColorAlpha(kBlue, 0.0);
    legend3->SetLineColor(0);
    legend3->SetTextSize(0.04);
    legend3->Draw("same");

    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptMin, ptMax));
    text1->Draw("same");
//    text5 -> Draw("same");
    fOut->cd();
    c5->Write(Form("%.1f_%.1f_mass_zoom", ptMin, ptMax));
    c5->Close();

    delete c5;
    delete legend3;
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::fitFunction() {
    TH1F* histo=(TH1F*)sigOrig->Clone(sigOrig->GetName());

    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", mFitRMin,mFitRMax);
    funLS->SetParameters(1.,-5.,1.,1.84,0.01);
//    funLS->SetLineStyle(7);
    funLS->SetLineStyle(1);
    funLS->SetParName(0,"intercept");
    funLS->SetParName(1,"slope");
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    funLS->SetParLimits(3,peakMin,peakMax);
    funLS->SetLineColor(9);

    histo->Fit(funLS,"RLI");

    TF1 *fgausLS = new TF1("fgausLS","gaus",mFitRMin,mFitRMax);
    fgausLS->SetParameters(funLS->GetParameter(2),funLS->GetParameter(3),funLS->GetParameter(4));

    TF1 *flinLS = new TF1("flinLS","pol1",mFitRMin,mFitRMax);
    flinLS->SetParameters(funLS->GetParameter(0),funLS->GetParameter(1));
    flinLS->SetLineStyle(7);
    flinLS->SetLineWidth(1);
    flinLS->SetLineColor(46);

//    Double_t SLSError=fgausLS->IntegralError(funLS->GetParameter(3)-nsigma*funLS->GetParameter(4), funLS->GetParameter(3)+nsigma*funLS->GetParameter(4));
    Double_t SLS=fgausLS->Integral(funLS->GetParameter(3)-nsigma*funLS->GetParameter(4), funLS->GetParameter(3)+nsigma*funLS->GetParameter(4))/binSize;
    Double_t SLSError=0;

    Double_t BLS=flinLS->Integral(funLS->GetParameter(3)-nsigma*funLS->GetParameter(4), funLS->GetParameter(3)+nsigma*funLS->GetParameter(4))/binSize;
    std::cout<<"Fitting without backround signal, background: "<<SLS<<" "<<BLS<<endl;
    Double_t signLS=abs(SLS)/sqrt(abs(SLS)+abs(2*BLS));
    std::cout<<"Significance: "<<signLS<<endl;

    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    gPad->SetMargin(0.1,0.05,0.13,0.08);
    gStyle->SetOptStat(0);
//    gStyle->SetOptFit(111);
    gStyle->SetStatX(0.949833);
    gStyle->SetStatY(0.920598);
    histo->Draw();
    gStyle->SetOptFit(111);
    gStyle->SetStatFontSize(6);

    histo->GetFunction("funLS")->SetLineColor(46);

    TPaveStats *st = (TPaveStats*)histo->FindObject("stats");
    st->SetX1NDC(0.7); //new x start position
    st->SetY1NDC(0.7); //new x end position
    st->SetX2NDC(0.949833); //new x start position
    st->SetY2NDC(0.920598); //new x end position

    TPaveText *textPub1 = new TPaveText(0.13,0.717,0.205,0.823,"brNDC");
    textPub1->SetTextSize(0.04);
    textPub1->SetLineColor(0);
    textPub1->SetShadowColor(0);
    textPub1->SetFillColor(0);
    textPub1->SetTextFont(42);
    textPub1->AddText(Form("Significance: %0.3g, S/B: %0.3g", signLS, SLS/BLS));
//    textPub1->AddText(Form("Raw yield: %0.2g #pm %0.2g", SLS, SLSError));
    textPub1->AddText(Form("Raw yield: %0.2g", SLS));
    textPub1->SetTextAlign(12);
    textPub1->Draw("same");

    TLegend *legendPub = new TLegend(0.127,0.817, 0.287, 0.907,"","brNDC");
    legendPub->AddEntry(histo, "Unlike-sign (US) signal", "pl");
    legendPub->AddEntry(histo->GetFunction("funLS"), "Gaussian + linear fit", "pl");
    legendPub->SetFillStyle(0);
    legendPub->SetLineColor(0);
    legendPub->SetTextSize(0.04);
    legendPub->Draw("same");


    flinLS->Draw("same");
//    fgausLS->Draw("same");
    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptMin, ptMax));
    text1->Draw("same");

    fOut->cd();
    c1->Write(Form("%.1f_%.1f_fit_fct", ptMin, ptMax));
    c1->Close();

    delete c1;
    delete fgausLS;
    delete flinLS;
    delete textPub1;
    delete histo;
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::fitComeOn() {
    TCanvas *c1 = new TCanvas("c1","c1",1200,900);

    fun0 = new TF1("fun0","pol1(0)+gaus(2)", mFitRMin, mFitRMax);
    fun0->SetParameters(1.,1.,1.,1.84,0.01);
    fun0->SetLineColor(2);
    fun0->SetLineStyle(7);
    fun0->SetParName(0,"intercept");
    fun0->SetParName(1,"slope");
    fun0->SetParName(2,"height");
    fun0->SetParName(3,"mean");
    fun0->SetParName(4,"sigma");
    fun0->SetParLimits(3, peakMin, peakMax);
    fun0->SetLineColor(9);

    sigSubtracted->Fit(fun0, "RLI");

    //getting results from the fit to function and parameters
    TF1 *resGaus = new TF1("resGaus","gaus(0)",mFitRMin,mFitRMax);
    resGaus->SetParameters(fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resGaus->SetLineColor(4);
    resGaus->SetLineStyle(9);

    TF1 *resLinear = new TF1("resLinear","pol1",mFitRMin,mFitRMax);
    resLinear->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
    resLinear->SetLineStyle(7);
    resLinear->SetLineWidth(1);
    resLinear->SetLineColor(46);

    mMean = abs(fun0->GetParameter(3));
    mMeanE = abs(fun0->GetParError(3));
    mSigma = abs(fun0->GetParameter(4));
    mSigmaE = abs(fun0->GetParError(4));
    mHeight = abs(fun0->GetParameter(2));
    mHeightE = abs(fun0->GetParError(2));

    //lines around peak, sigma same as in the raw yield
    leftLine=new TLine(mMean - nsigma*mSigma, sigSubtracted->GetMaximum(), mMean - nsigma*mSigma, sigSubtracted->GetMinimum());
    leftLine->SetLineColor(46);
    rightLine=new TLine(mMean + nsigma*mSigma, sigSubtracted->GetMaximum(), mMean + nsigma*mSigma, sigSubtracted->GetMinimum());
    rightLine->SetLineColor(46);

    //residual backg subtraction
    sigSubtractedResidualBckg=(TH1F*)sigSubtracted->Clone(sigSubtracted->GetName());
    setHistoStyle(sigSubtractedResidualBckg, 15, 20);
    sigSubtractedResidualBckg->Add(resLinear, -1);

    bckgAddedResidualBckg=(TH1F*)bckgToWork->Clone(bckgToWork->GetName());
    setHistoStyle(bckgAddedResidualBckg, 15, 5);
    bckgAddedResidualBckg->Add(resLinear, 1);

    //Evaluation of significances and stuff:
    //Peak range, bins
    Int_t intLow = sigSubtractedResidualBckg->FindBin(mMean-nsigma*mSigma);
    Int_t intUp = sigSubtractedResidualBckg->FindBin(mMean+nsigma*mSigma);

    std::cout<<"Integral range: "<<intLow<<" "<<intUp<<endl;
    std::cout<<"Integral real range from centre of bins: "<<sigSubtracted->GetBinCenter(intLow)<<" "<<sigSubtracted->GetBinCenter(intUp)<<endl;

    rawYield = sigSubtractedResidualBckg->IntegralAndError(intLow, intUp, rawYieldError, "");

    Double_t B = bckgAddedResidualBckg->Integral(intLow, intUp,"");
    SoverB = rawYield/B;
    significanceBins = rawYield/TMath::Sqrt(rawYield+2*B);

    std::cout<<"raw yield with error from histo (IntegralAndError): "<<rawYield<<" "<<rawYieldError<<endl;
    std::cout<<"S: "<<rawYield<<endl;
    std::cout<<"B: "<<B<<endl;
    std::cout<<"S/B: "<<SoverB<<endl;
    std::cout<<"start fit end"<<endl;
    c1->Close();
    delete c1;
    delete resGaus;
    delete resLinear;
}

//____________________________________________________________________________________________________________________________
void FitD0Peak::makeTupleCalculations(){
    if (mTupleSignal==NULL || mTupleBackground==NULL) {
        std::cout<<"tuple not specified. end."<<endl;
        return;
    }



}

//____________________________________________________________________________________________________________________________
void FitD0Peak::addMixedEventBckg(TH1F* histoMxdEv) {
    //mxd event background neet to be scaled to the signal
    isMxdEv=true;
    scale=true;
    mxdOrig=(TH1F*)histoMxdEv->Clone(histoMxdEv->GetName());
    mxdOrig->Rebin(10);
    setHistoStyle(mxdOrig, 8, 26);

    Float_t intLowLow = 1.7;
    Float_t intLowUp = 1.8;
    Float_t intUpLow = 1.92;
    Float_t intUpUp = 2.;

    Double_t mixIntegral = mxdOrig->Integral(mxdOrig->FindBin(intLowLow), mxdOrig->FindBin(intLowUp), "") + mxdOrig->Integral(mxdOrig->FindBin(intUpLow), mxdOrig->FindBin(intUpUp), ""); //integrals of original histos
    Double_t sigIntegral = sigOrig->Integral(sigOrig->FindBin(intLowLow), sigOrig->FindBin(intLowUp), "") + sigOrig->Integral(sigOrig->FindBin(intUpLow), sigOrig->FindBin(intUpUp), ""); //integrals of original histos

    mxdOrig->Scale(sigIntegral/mixIntegral);
    cout<<"mxd added"<<endl;
}

//____________________________________________________________________________________________________________________________
FitD0Peak::~FitD0Peak() {
    // destructor
    delete bckgOrig;
    delete sigOrig;
    delete sigSubtractedResidualBckg;
    delete sigSubtracted;
    delete rightLine;
    delete leftLine;
    delete fun0;
    delete text1;

    fOut->Close();
    delete fOut;
}