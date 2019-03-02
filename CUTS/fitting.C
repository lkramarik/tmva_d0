#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
using namespace std;
using namespace TMath;

Double_t fit(TH1F* , TH1F* , Double_t, Double_t, bool , bool, TString, float);

void fitting(TString input = "D0_bdt_cuts.root", bool scale = false, bool subtract = true,  Double_t ptminText = 2, Double_t ptmaxText = 5, float bdtcut = 0.1) {
    cout << input << endl;
    TFile *data = new TFile(input, "r");
    TList *listS = (TList *) data->Get("hists_signal");
    TList *listB = (TList *) data->Get("hists_background");

    //just number of pairs:
    cout << Form("hB_%.1f_%.1f", ptminText, ptmaxText) << endl;

    TH1F *hInvMassBack = (TH1F*) listB->FindObject(Form("D_mass_%.4f_background", bdtcut));
    TH1F *hInvMassSign = (TH1F*) listS->FindObject(Form("D_mass_%.4f_signal", bdtcut));

    fit(hInvMassSign, hInvMassBack, ptminText,  ptmaxText, scale, subtract, input, bdtcut);
}


Double_t fit(TH1F* hInvMassSign, TH1F* hInvMassBack, Double_t ptminText, Double_t ptmaxText, bool scale, bool subtract, TString input, float bdtcut){
//fit(TH1F* hInvMassSign, TH1F* hInvMassBack, Double_t ptminText, Double_t ptmaxText, bool scale, bool subtract, TString input, float bdtcut){
    int rebin = 10; //this is good
    Double_t nsigma = 3;
    Float_t fitRMin = 1.7;
    Float_t fitRMax = 2.;
    TString intLowLow = "1.7";
    TString intLowUp = "1.8";
    TString intUpLow = "1.92";
    TString intUpUp = "2.2";
    const float rotwthmin = 1.84; // peak mean fitting range
    const float rotwthmax = 1.88; //peak mean fitting range

    const float MKSize    = 1.;

    const int Nbins = 2000/rebin;
    Double_t binSize = (2.4-0.4)*rebin/2000;

    hInvMassBack -> Sumw2();
    hInvMassSign -> Sumw2();
    hInvMassSign -> Rebin(rebin);
    hInvMassBack -> Rebin(rebin);

    hInvMassSign -> SetMarkerColor(46);
    hInvMassSign -> SetLineColor(46);

    Double_t hBackIntegral = hInvMassBack -> Integral(hInvMassBack -> FindBin(intLowLow.Atof()), hInvMassBack -> FindBin(intLowUp.Atof()), "") + hInvMassBack -> Integral(hInvMassBack -> FindBin(intUpLow.Atof()), hInvMassBack -> FindBin(intUpUp.Atof()), "");
    Double_t hSignIntegral = hInvMassSign -> Integral(hInvMassSign -> FindBin(intLowLow.Atof()), hInvMassSign -> FindBin(intLowUp.Atof()), "") + hInvMassSign -> Integral(hInvMassSign -> FindBin(intUpLow.Atof()), hInvMassSign -> FindBin(intUpUp.Atof()), "");

    cout<<hInvMassBack -> GetNbinsX()<<endl;
    cout<<hInvMassSign -> GetNbinsX()<< endl;

    cout<<"Backgroung Integral: "<<hBackIntegral<<endl;
    cout<<"Signal Integral :"<<hSignIntegral<<endl;

    TF1 *funLS = new TF1("funLS","pol1(0)+gaus(2)", fitRMin,fitRMax);
    funLS->SetParameters(1.,1.,1.,1.84,0.01);
    funLS->SetLineColor(2);
//    funLS->SetLineStyle(7);
    funLS->SetLineStyle(1);
    funLS->SetParName(0,"intercept");
    funLS->SetParName(1,"slope");
    funLS->SetParName(2,"height");
    funLS->SetParName(3,"mean");
    funLS->SetParName(4,"sigma");
    funLS->SetParLimits(3,rotwthmin,rotwthmax);
    funLS->SetLineColor(9);

    double xLS[Nbins],yLS[Nbins],yLSe[Nbins], xWS[Nbins], yWS[Nbins], yWSe[Nbins];
    for (int ibi=0; ibi<Nbins; ibi++) {
        xLS[ibi]  = hInvMassSign -> GetBinCenter(ibi);
        yLS[ibi]  = hInvMassSign -> GetBinContent(ibi);
        yLSe[ibi] = hInvMassSign -> GetBinError(ibi);
        xWS[ibi] = hInvMassBack -> GetBinCenter(ibi);
        yWS[ibi] = hInvMassBack -> GetBinContent(ibi);
        yWSe[ibi] =  hInvMassBack -> GetBinError(ibi);
    }

    TGraphErrors *gBackground = new TGraphErrors(Nbins,xWS,yWS,0,yWSe);
    gBackground->SetMarkerStyle(20);
    gBackground->SetMarkerSize(0.9);
    gBackground->SetMarkerColor(kBlack);
    gBackground->SetLineColor(kBlack);
    gBackground->GetYaxis()->SetTitle("Raw Counts");
    gBackground->GetYaxis()->SetTitleOffset(1.1);
    gBackground->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    gBackground->SetTitle("");
    gBackground->PaintStats(funLS);
    gBackground->GetXaxis()->SetRangeUser(fitRMin,fitRMax);
//    gBackground->GetYaxis()->SetRangeUser(-70.,500.0);

    TGraphErrors *gLS = new TGraphErrors(Nbins,xLS,yLS,0,yLSe);
    gLS->SetMarkerStyle(20);
    gLS->SetMarkerSize(0.9);
    gLS->SetMarkerColor(kBlack);
    gLS->SetLineColor(kBlack);
    gLS->GetYaxis()->SetTitle("Raw Counts");
    gLS->GetYaxis()->SetTitleOffset(1.1);
    gLS->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    gLS->SetTitle("");
    gLS->PaintStats(funLS);
    gLS->GetXaxis()->SetRangeUser(fitRMin,fitRMax);
//    gLS->GetYaxis()->SetRangeUser(-70.,500.0);
    gLS->Fit(funLS,"OR");

    TF1 *fgausLS    = new TF1("fgausLS","gaus",fitRMin,fitRMax);
    fgausLS->SetParameters(funLS->GetParameter(2),funLS->GetParameter(3),funLS->GetParameter(4));

    TF1 *flinLS   = new TF1("flinLS","pol1",fitRMin,fitRMax);
    flinLS->SetParameters(funLS->GetParameter(0),funLS->GetParameter(1));
    flinLS->SetLineStyle(7);
    flinLS->SetLineWidth(1);
    flinLS->SetLineColor(46);

    Double_t SLS=fgausLS->Integral(funLS->GetParameter(3)-nsigma*funLS->GetParameter(4), funLS->GetParameter(3)+nsigma*funLS->GetParameter(4))/binSize;
    Double_t BLS=flinLS->Integral(funLS->GetParameter(3)-nsigma*funLS->GetParameter(4), funLS->GetParameter(3)+nsigma*funLS->GetParameter(4))/binSize;
    cout<<"Fitting without backround signal, background: "<<SLS<<" "<<BLS<<endl;
    Double_t signLS=abs(SLS)/sqrt(abs(SLS)+abs(BLS));
    cout<<"Significance: "<<signLS<<endl;

    if(scale)   hInvMassBack -> Scale(hSignIntegral/hBackIntegral);
    float err[Nbins], errS[Nbins];
    for (int j=0; j<Nbins; j++) {
        err[j] = hInvMassBack -> GetBinError(j);
        errS[j] = hInvMassSign -> GetBinError(j);
    }

//    for (j=0; j<Nbins; j++) {
//        err[j] = err[j]*hSignIntegral/hBackIntegral;
//        hInvMassBack -> SetBinError(j, err[j]);
//    }   //scaling error

    TH1F *hInvMassSignOrig = (TH1F*)hInvMassSign->Clone("signal_orig");
    Double_t nentriesSig = hInvMassSignOrig->Integral(hInvMassSignOrig->FindBin(1.7), hInvMassSignOrig->FindBin(2),"");

    if(subtract)   hInvMassSign -> Add(hInvMassBack,-1);

    TF1 *fun0 = new TF1("fun0","pol1(0)+gaus(2)", fitRMin,fitRMax);
    fun0->SetParameters(1.,1.,1.,1.84,0.01);
    fun0->SetLineColor(2);
    fun0->SetLineStyle(7);
    fun0->SetParName(0,"intercept");
    fun0->SetParName(1,"slope");
    fun0->SetParName(2,"height");
    fun0->SetParName(3,"mean");
    fun0->SetParName(4,"sigma");
    fun0->SetParLimits(3,rotwthmin,rotwthmax);
    fun0->SetLineColor(9);

    double mm[Nbins],ym[Nbins],yme[Nbins],ym1[Nbins], mme[Nbins];
    for (int ib=0; ib<Nbins; ib++) {
        mme[ib] = binSize/2;
        mm[ib]  = hInvMassSign -> GetBinCenter(ib);
        ym[ib]  = hInvMassSign -> GetBinContent(ib);
        yme[ib] = hInvMassSign -> GetBinError(ib);
    }
    TGraphErrors *gm = new TGraphErrors(Nbins,mm,ym,mme  ,yme);
    gm->SetMarkerStyle(20);
    gm->SetMarkerSize(2.);
    gm->SetMarkerColor(9);
    gm->SetLineColor(9);
    gm->Fit(fun0,"OR");
    gm->GetYaxis()->SetTitle("Counts");
    gm->GetYaxis()->SetTitleOffset(0.7);
    gm->GetYaxis()->SetTitleSize(0.045);
    gm->GetXaxis()->SetTitleSize(0.045);
    gm->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    gm->GetXaxis()->CenterTitle(kTRUE);
    gm->SetTitle("");
    TF1 *resGaus = new TF1("resGaus","gaus(0)",fitRMin,fitRMax);
    resGaus->SetParameters(fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resGaus->SetLineColor(4);
    resGaus->SetLineStyle(9);
    gm->PaintStats(resGaus);
    gm->GetXaxis()->SetRangeUser(fitRMin,fitRMax);
//    gm->GetYaxis()->SetRangeUser(-8.,30.0);

    Double_t mean = abs(fun0->GetParameter(3));
    Double_t sigma = abs(fun0->GetParameter(4));
    cout<<"mean from fit: "<<mean<<endl;
    cout<<"sigma from fit: "<<sigma<<endl;

    TF1 *resfunm = new TF1("resfunm","pol1",fitRMin,fitRMax);
    resfunm->SetParameters(fun0->GetParameter(0),fun0->GetParameter(1));
    resfunm->SetLineStyle(7);
    resfunm->SetLineWidth(1);
    resfunm->SetLineColor(46);
    //fun0->Draw("same");

    TF1 *resfunm1 = new TF1("resfunm1","gaus(0)",fitRMin,fitRMax);
    resfunm1->SetParameters(fun0->GetParameter(2),fun0->GetParameter(3),fun0->GetParameter(4));
    resfunm1->SetLineColor(4);
    resfunm1->SetLineStyle(9);

    double mSig[Nbins],ySig[Nbins],ySige[Nbins],ySig1[Nbins];
    for (int h = 0; h < Nbins; ++h) {
        mSig[h]=mm[h];
        ySig[h]=ym[h]-resfunm->Eval(mm[h]);
        ySige[h]=yme[h];//-resfunm1->Eval(mm[h]);
    }
    TGraphErrors *gSignal = new TGraphErrors(Nbins,mSig,ySig,0,ySige);
    gSignal->SetLineColor(kBlack);
    gSignal->Fit(resfunm1,"OR");

    Double_t integral_function_yield = resfunm1->Integral(mean-nsigma*sigma, mean+nsigma*sigma)/binSize;
    Double_t integral_function_yield_error = resfunm1->IntegralError(mean-nsigma*sigma, mean+nsigma*sigma);

    //Subtraction of residual background and fit by gauss
//    auto mg = new TMultiGraph();
//    mg->SetTitle(" ;p_{T} [GeV/c];Raw counts");
//    mg->Add(gm);
//    mg->Add(gSignal);
//    mg->Draw("ap");
//    mg->GetXaxis()->SetLimits(1.72,2.0);

//    resfunm->Draw("same");
//    resfunm1->Draw("same");


    //Evaluation of significances and stuff
    Double_t intLow = hInvMassSignOrig -> FindBin(mean-nsigma*sigma);
    Double_t intUp = hInvMassSignOrig -> FindBin(mean+nsigma*sigma);
    cout<<"Integral range: "<<intLow<<" "<<intUp<<endl;
    cout<<"Integral range: "<<hInvMassSign->FindBin(mean-nsigma*sigma)<<" "<<hInvMassSign->FindBin(mean+nsigma*sigma)<<endl;
    cout<<"Integral real range of centre of bins: "<<hInvMassSignOrig->GetBinCenter(intLow)<<" "<<hInvMassSignOrig->GetBinCenter(intUp)<<endl;
    cout<<"Integral raw yield from fct.: "<<integral_function_yield<<endl;
    cout<<"Integral raw yield from fct. error: "<<integral_function_yield_error/binSize<<endl;

    Int_t nbinsInPeak = intUp-intLow;
    Double_t sumBins=0, sumBinsErr=0;
    for (int i = intLow; i < intUp+1; ++i) {
        sumBins += hInvMassSign->GetBinContent(i);
        sumBinsErr += (hInvMassSign->GetBinError(i))*(hInvMassSign->GetBinError(i));
    }
    cout<<"Raw yield manual counting: "<<sumBins<<" "<<sqrt(sumBinsErr)<<endl;

    //Tgraph of signal with residual back subtracted
    for(int j1=0; j1<Nbins; j1++) ym1[j1] = ym[j1] - resfunm->Eval(mm[j1]);
    TGraphErrors *gm1 = new TGraphErrors(Nbins,mm,ym1,0,yme);
    gm1->SetMarkerStyle(24);
    gm1->SetMarkerSize(MKSize-0.1);
    gm1->SetMarkerColor(2);
    gm1->SetLineColor(2);
//    gm1->Draw("same,AP");

    TH1F* hInvMassSignClean = new TH1F("hInvMassSignClean", "hInvMassSignClean", 2000, 0.4, 2.4);
    TH1F* hInvMassBackClean = new TH1F("hInvMassBackClean", "hInvMassBackClean", 2000, 0.4, 2.4);
    hInvMassSignClean->Rebin(rebin);
    hInvMassBackClean->Rebin(rebin);

    //Histos after residual background subtr.
    for (Int_t bin = 1; bin < hInvMassSignClean->GetNbinsX(); ++bin) {
        hInvMassSignClean->SetBinContent(bin, hInvMassSign->GetBinContent(bin)- resfunm->Eval(hInvMassSign->GetBinCenter(bin)));
        hInvMassSignClean->SetBinError(bin, hInvMassSign->GetBinError(bin));
        hInvMassBackClean->SetBinContent(bin, hInvMassBack->GetBinContent(bin)+ resfunm->Eval(hInvMassBack->GetBinCenter(bin)));
        hInvMassBackClean->SetBinError(bin, hInvMassBack->GetBinError(bin));
    }
    Double_t rawYieldError;
    Double_t rawYield = hInvMassSignClean -> IntegralAndError(hInvMassSignClean->FindBin(mean-nsigma*sigma), hInvMassSignClean->FindBin(mean+nsigma*sigma) , rawYieldError, "");
    Double_t S = rawYield;
    Double_t B = hInvMassBackClean -> Integral(hInvMassBackClean->FindBin(mean-nsigma*sigma), hInvMassBackClean->FindBin(mean+nsigma*sigma) ,"");
    Double_t SoverB = S/B;
    Double_t significance_bins = S/TMath::Sqrt(S+2*B);

    cout<<"raw yield with error from histo (IntegralAndError): "<<rawYield<<" "<<rawYieldError<<endl;
    cout<<"S: "<<S<<endl;
    cout<<"B: "<<B<<endl;
    cout<<"S/B: "<<SoverB<<endl;
    Double_t significance_fit = integral_function_yield/integral_function_yield_error ;
    cout<<"Significance from fit: "<<abs(significance_fit)<<endl;
    cout<<"Significance from number of entries in bins: "<<significance_bins<<endl;

    TPaveText *text1 = new TPaveText(0.719,0.9229,0.998,0.980,"brNDC");
    text1 -> SetTextSize(0.04);
    text1 -> SetShadowColor(0);
    text1 -> SetLineColor(0);
    text1 -> SetFillColor(0);
    text1 -> AddText("d+Au 200 GeV");

////// graphs signal, bck nofct
    TCanvas *c3 = new TCanvas("c3","c3",1200,1000);
//    gStyle->SetOptFit(1);
//    gStyle->SetStatY(0.899);
//    gStyle->SetStatY(0.899);
//    gStyle->SetStatX(0.9);

    gm->Draw("ap");
//    rightline->Draw("");
//    leftline->Draw("same");

    TPaveText *text = new TPaveText(0.172,0.735,0.443,0.8300,"brNDC");
    text -> SetTextSize(0.025);
    text->SetTextColor(39);
    text -> SetLineColor(0);
    text -> SetShadowColor(0);
    text -> SetFillColor(0);
//    text -> AddText(input);
//    text -> AddText(Form("N_entries in sig. 1.7-2.0 GeV/c^{2}: %g", nentriesSig));
    text -> AddText(Form("Bin size: %g GeV/c^{2}", binSize));
    TString paveSc;
    if (scale) {TString paveSc = "Scaling integral: " + intLowLow + "-" + intLowUp +", " + intUpLow + "-" + intUpUp + " GeV/c^{2}";}
    else TString paveSc = "No scaling";
    text -> AddText(paveSc);
//    text -> AddText(Form("Significance from fit: %0.3g", significance_fit));
//    text -> AddText(Form("Significance: %0.3g", significance_bins));
    text -> Draw("same");

    TPaveText *text3 = new TPaveText(0.0,0.01,0.27,0.0350,"brNDC");
    text3 -> SetTextSize(0.02);
    text3->SetTextColor(39);
    text3 -> SetLineColor(0);
    text3 -> SetShadowColor(0);
    text3 -> SetFillColor(0);
    text3 -> AddText(input);
    text3 -> Draw("same");

    TPaveText *text4 = new TPaveText(0.29,0.848,0.229,0.873,"brNDC");
    text4 -> SetTextSize(0.03);
    text4 -> SetLineColor(0);
    text4 -> SetShadowColor(0);
    text4 -> SetFillColor(0);
    text4 -> AddText(Form("Significance: %0.3g, S/B: %0.4g", significance_bins, S/B));
    text4 -> Draw("same");

    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.04);
    tx2.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptminText, ptmaxText));

    text1 -> Draw("same");
    gStyle->SetOptFit();

////// graphsLS signal
    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
//    gStyle->SetOptFit(1);
//    gStyle->SetStatY(0.899);
//    gStyle->SetStatX(0.9);

    gLS->Draw("ap");

    TPaveText *textLS = new TPaveText(0.29,0.848,0.229,0.873,"brNDC");
    textLS -> SetTextSize(0.03);
    textLS -> SetLineColor(0);
    textLS -> SetShadowColor(0);
    textLS -> SetFillColor(0);
    textLS -> AddText(Form("Significance: %0.3g, S/B: %0.4g", signLS, signLS));
    textLS -> Draw("same");

    TLatex txR;
    txR.SetNDC();
    txR.SetTextSize(0.04);
    txR.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptminText, ptmaxText));

    flinLS->Draw("same");
    fgausLS->Draw("same");
    text1 -> Draw("same");
//         gStyle->SetOptFit();


/////////// PUBLSISH
    TCanvas *cP = new TCanvas("cPublish","cPublish",1200,1000);
    gPad->SetMargin(0.1,0.05,0.13,0.08);

    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.949833);
    gStyle->SetStatY(0.920598);

    gm->SetTitle(" ;Invariant mass, m_{K#pi} [GeV/c^{2}];Counts");
    gm->GetXaxis()->SetLabelSize(0.045);
    gm->GetYaxis()->SetLabelSize(0.045);
    gm->GetXaxis()->SetTitleSize(0.05);
    gm->GetYaxis()->SetTitleSize(0.05);
    gm->GetXaxis()->CenterTitle(kTRUE);
    gm->GetXaxis()->SetTitleOffset(1.0);
    gm->GetYaxis()->SetTitleOffset(1.);
    gm->GetYaxis()->CenterTitle(kTRUE);

    gm->Draw("apz");
    TPaveText *textPub1 = new TPaveText(0.13,0.717,0.205,0.823,"brNDC");
    textPub1 -> SetTextSize(0.04);
    textPub1 -> SetLineColor(0);
    textPub1 -> SetShadowColor(0);
    textPub1 -> SetFillColor(0);
    textPub1 -> SetTextFont(42);
    textPub1 -> AddText(Form("Significance: %0.3g, S/B: %0.3g", significance_bins, SoverB));
    textPub1 -> AddText(Form("Raw yield: %0.2g #pm %0.2g", S, rawYieldError));
    textPub1 -> SetTextAlign(12);
//    textPub -> Draw("same");
    textPub1 -> Draw("same");
    TLegend *legendPub = new TLegend(0.127,0.817, 0.287, 0.907,"","brNDC");
    legendPub -> AddEntry(gm, "US - LS", "pl");
    legendPub -> AddEntry(fun0, "Gaussian + linear fit", "pl");
    legendPub -> SetFillStyle(0);
    legendPub -> SetLineColor(0);
    legendPub -> SetTextSize(0.04);
    legendPub -> Draw("same");
//    0.719,0.9229,0.998,0.980
    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptminText, ptmaxText));
    text1 -> Draw("same");
///////////////PLOT
    TCanvas *c4 = new TCanvas("c4","c4",1200,1000);
    gPad->SetLeftMargin(0.15);

    hInvMassBack->SetTitle("");
    hInvMassBack->SetStats(0);
    hInvMassBack->SetLineColor(46);
    hInvMassBack->SetMarkerColor(46);
    hInvMassBack->SetMarkerStyle(2);
    hInvMassSign->SetTitle("");
    hInvMassSign->SetStats(0);
    hInvMassSign->SetLineColor(9);
    hInvMassSign->SetMarkerStyle(20);
    hInvMassSign->SetMarkerColor(9);
    hInvMassSign->GetYaxis()->SetTitle("Raw Counts");
    hInvMassSign->GetYaxis()->SetTitleOffset(1.7);
    hInvMassSign->GetXaxis()->SetTitle("Mass_{K#pi} (GeV/c^{2})");
    hInvMassSign->GetXaxis()->SetLabelFont(42);
    hInvMassSign->GetXaxis()->SetTitleFont(42);
    hInvMassSign->GetYaxis()->SetLabelFont(42);
    hInvMassSign->GetYaxis()->SetTitleFont(42);
    hInvMassSign -> GetXaxis()->SetRangeUser(0.65,2.4);
    hInvMassSign -> Draw("");

    TLine *leftline1 = new TLine(mean-3*sigma,hInvMassSign->GetMinimum(),mean-3*sigma,hInvMassSign->GetMaximum());
    leftline1->SetLineStyle(9);
    leftline1->SetLineColor(28);
    leftline1->Draw("same");

    TLine *leftline2 = new TLine(mean+3*sigma,hInvMassSign->GetMinimum(),mean+3*sigma,hInvMassSign->GetMaximum());
    leftline2->SetLineStyle(9);
    leftline2->SetLineColor(28);
    leftline2->Draw("same");

//    legend -> AddEntry(hInvMassBack, Form("scaled background (%g-%g & %g-%g GeV/c^{2})", intLowLow.Atof(), intLowUp.Atof(), intUpLow.Atof(), intUpUp.Atof()), "pl");
    TLatex tx1;
    tx1.SetNDC();
    tx1.SetTextSize(0.04);
    tx1.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptminText, ptmaxText));
    text1 -> Draw("same");

    TLegend *legend = new TLegend(0.40,0.845, 0.560, 0.897,"","brNDC");
    legend -> AddEntry(hInvMassSignOrig, "Unlike-sign - like-sign pairs", "pl");
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    legend -> Draw("same");

/////////////SIGNAL AND BACKGROUND
    TCanvas *c5 = new TCanvas("c5","c5",1200,1000);
    gPad->SetMargin(0.1,0.05,0.13,0.08);
//    gPad->SetLeftMargin(0.15);
//    gPad->SetMargin(0.12,0.1,0.18,0.05);
//    gStyle->SetOptFit(0);
//    gStyle->SetStatY(0.899);
//    gStyle->SetStatX(0.9);
//    gStyle->SetTitleFontSize(.045);
//    gStyle->SetLabelSize(.04);

    hInvMassBack -> SetTitle("");
    hInvMassBack -> SetStats(0);
    hInvMassBack -> SetLineColor(kRed);
    hInvMassBack -> SetMarkerColor(kRed);
    hInvMassBack -> SetMarkerStyle(25);
    hInvMassBack -> SetMarkerSize(2);

    TPaveText *text5 = new TPaveText(0.74,0.71,0.81,0.94,"brNDC");
    text5 -> SetTextSize(0.05);
    text5 -> SetLineColor(0);
    text5 -> SetShadowColor(0);
    text5 -> SetFillColor(0);
    text5 -> SetTextFont(42);
    text5 -> AddText("STAR Preliminary");
    text5 -> AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
    text5 -> AddText(Form("%3.1f < p_{T} < %3.1f GeV/c", ptminText, ptmaxText));
    TText *t = text5->GetLineWith("STAR");
    t->SetTextColor(kRed);

    hInvMassSignOrig->SetTitle("");
    hInvMassSignOrig->SetStats(0);
    hInvMassSignOrig->SetLineColor(kBlack);
    hInvMassSignOrig->SetMarkerStyle(20);
    hInvMassSignOrig->SetMarkerSize(2.);
    hInvMassSignOrig->SetMarkerColor(kBlack);
    hInvMassSignOrig->GetYaxis()->SetTitle("Counts");
    hInvMassSignOrig->GetYaxis()->SetTitleOffset(1.0);

    hInvMassSignOrig->GetXaxis()->SetTitleOffset(1.);
    hInvMassSignOrig->GetYaxis()->CenterTitle(kTRUE);
    hInvMassSignOrig->GetXaxis()->SetTitle("Invariant mass, m_{K#pi} [GeV/c^{2}]");
    hInvMassSignOrig->GetXaxis()->SetLabelFont(42);
//    hInvMassSignOrig->GetXaxis()->SetLabelOffset(0.015);
    hInvMassSignOrig->GetXaxis()->SetTitleFont(42);
    hInvMassSignOrig->GetXaxis()->CenterTitle(kTRUE);
    hInvMassSignOrig->GetYaxis()->SetLabelFont(42);
    hInvMassSignOrig->GetYaxis()->SetTitleFont(42);
    hInvMassSignOrig->GetYaxis()->SetTitleSize(0.05);
    hInvMassSignOrig->GetXaxis()->SetTitleSize(0.05);
    hInvMassSignOrig->GetYaxis()->SetLabelSize(0.045);
    hInvMassSignOrig->GetXaxis()->SetLabelSize(0.045);

    hInvMassSignOrig->GetXaxis()->SetRangeUser(fitRMin,fitRMax);
    hInvMassSignOrig->Draw("");
    hInvMassBack->Draw("same");
    hInvMassSignOrig->Draw("same");

    TLegend *legend3 = new TLegend(0.1137,0.812, 0.351, 0.913,"","brNDC");
    legend3 -> AddEntry(hInvMassBack, "Like-sign (LS) background", "pl");
    legend3 -> AddEntry(hInvMassSignOrig, "Unlike-sign (US) signal", "pl");
    legend3 -> SetFillStyle(0);
    legend3->SetFillColorAlpha(kBlue, 0.0);
    legend3 -> SetLineColor(0);
    legend3 -> SetTextSize(0.04);
    legend3 -> Draw("same");
    tx2.DrawLatex(0.1,0.940,Form("p_{T}: %3.1f-%3.1f GeV/c", ptminText, ptmaxText));
    text1 -> Draw("same");
//    text5 -> Draw("same");
///////////////////////////////////////////////////////////////////////////////////////

    TFile *fOut = new TFile(Form("signals_%.3f", bdtcut)+input,"recreate");
    c3 -> Write(Form("%.1f_%.1f_signal",ptminText, ptmaxText));
    c4 -> Write(Form("%.1f_%.1f_mass", ptminText, ptmaxText));
    c5 -> Write(Form("%.1f_%.1f_mass_zoom", ptminText, ptmaxText));
    c1 -> Write(Form("%.1f_%.1f_signal_fct", ptminText, ptmaxText));
    cP -> Write(Form("%.1f_%.1f_publ",ptminText, ptmaxText));
    fOut->Close();



//    cP->SaveAs(Form("%.1f_%.1f_bdt_%.3f_.png",ptminText, ptmaxText, bdtcut));

    ofstream yieldsF;
    yieldsF.open("raw_yields_fct.txt", std::ios::app);
    yieldsF<<ptminText<<" "<<ptmaxText<<" "<<integral_function_yield/binSize<<" "<<integral_function_yield_error/binSize<<endl;
    yieldsF.close();

    ofstream yields;
    yields.open("raw_yields.txt", std::ios::app);
    yields<<ptminText<<" "<<ptmaxText<<" "<<rawYield<<" "<<rawYieldError<<endl;
    yields.close();
    c1->Close();
    c3->Close();
    c4->Close();
    c5->Close();
//    cP->Close();

    return significance_bins;
}


//void plotDCA() {
//    TString folder = "";
////    TString folder = "res_ntp/p17id/";
//    TString input = "res_plots_signal.root";
//    TFile* data = new TFile(folder + input ,"r");
//
//    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);
//
////    hEvents = (TH1F*)data -> Get("hcosTheta");
//    hEvents = (TH1F*)data -> Get("hk_dca");
//    hEvents1 = (TH1F*)data -> Get("hpi1_dca");
//    hEvents ->GetYaxis()->SetTitleOffset(1.5);
////    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
//    hEvents -> SetTitle("");
//    hEvents -> SetStats(0);
////    hEvents -> SetFillColor(18);
////    hEvents1 -> SetFillColor(18);
//    hEvents->SetLineWidth(3); //2
//    hEvents1->SetLineWidth(3); //2
//
//    hEvents -> SetLineColor(46);
//    hEvents1 -> SetLineColor(9);
//    hEvents ->GetXaxis()->SetRangeUser(0.,0.04);
////    hEvents ->GetXaxis()-> SetTitle("cos(pointing angle diff.) [-]");
//    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
//    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
////    hEvents -> SetFillStyle(3344);
//
////    gPad->SetLogy();
//    gPad->SetLeftMargin(0.15);
//    TLegend *legend = new TLegend(0.397,0.82, 0.55, 0.87,"","brNDC");
//    legend -> SetFillStyle(0);
//    legend -> SetLineColor(0);
//    legend -> SetTextSize(0.035);
//    legend -> AddEntry(hEvents, "kaons", "pl");
//    legend -> AddEntry(hEvents1, "pions", "pl");
//
//    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
//    leftline2->SetLineWidth(3);
//    leftline2->SetLineStyle(9);
//    leftline2->SetLineColor(28);
//
////    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
////    text1 -> SetTextSize(0.04);
////    text1 -> SetShadowColor(0);
////    text1 -> SetLineColor(0);
////    text1 -> SetFillColor(0);
////    text1 -> AddText("d+Au 200 GeV");
//
//    TPaveText *text2 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
//    text2 -> SetTextSize(0.03);
//    text2 -> SetShadowColor(0);
//    text2 -> SetLineColor(0);
//    text2 -> SetFillColor(0);
//    text2->SetTextColor(28);
//    text2 -> AddText("cuts I., II.");
//
////    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
////    leftline2->SetLineWidth(2);
////    leftline2->SetLineStyle(9);
////    leftline2->SetLineColor(8);
//
//    hEvents -> Draw();
//    hEvents1 -> Draw("same");
//    legend -> Draw("same");
//    leftline2->Draw("same");
////    text1->Draw("same");
//    text2->Draw("same");
//}

//void plot() {
//    TString folder = "";
////    TString folder = "res_ntp/p17id/";
//    TString input = "res_plots_signal.root";
//    TFile* data = new TFile(folder + input ,"r");
//
//    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);
//
//    hEvents = (TH1F*)data -> Get("hcosTheta");
////    hEvents = (TH1F*)data -> Get("hk_dca");
////    hEvents1 = (TH1F*)data -> Get("hpi1_dca");
//    hEvents ->GetYaxis()->SetTitleOffset(1.9);
////    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
//    hEvents -> SetTitle("");
//    hEvents -> SetStats(0);
//    hEvents -> SetFillColor(18);
//    hEvents->SetLineWidth(3); //2
//
//    hEvents -> SetLineColor(46);
//    hEvents ->GetXaxis()->SetRangeUser(0.75,1);
//    hEvents ->GetXaxis()-> SetTitle("cos(pointing angle diff.) [-]");
////    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
//    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
////    hEvents -> SetFillStyle(3344);
//
////    gPad->SetLogy();
//    gPad->SetLeftMargin(0.15);
//    TLegend *legend = new TLegend(0.397,0.82, 0.55, 0.87,"","brNDC");
//    legend -> SetFillStyle(0);
//    legend -> SetLineColor(0);
//    legend -> SetTextSize(0.035);
////    legend -> AddEntry(hEvents, "kaons", "pl");
////    legend -> AddEntry(hEvents1, "pions", "pl");
//
//    TLine *leftline2 = new TLine(0.8,hEvents->GetMinimum()-0.05*hEvents->GetMinimum(),0.8,hEvents->GetMaximum()+0.04*hEvents->GetMaximum());
//    leftline2->SetLineWidth(3);
//    leftline2->SetLineStyle(9);
//    leftline2->SetLineColor(28);
//
//    TLine *leftline3 = new TLine(0.95,hEvents->GetMinimum()-0.05*hEvents->GetMinimum(),0.95,hEvents->GetMaximum()+0.04*hEvents->GetMaximum());
//    leftline3->SetLineWidth(3);
//    leftline3->SetLineStyle(9);
//    leftline3->SetLineColor(9);
//
////    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
////    text1 -> SetTextSize(0.04);
////    text1 -> SetShadowColor(0);
////    text1 -> SetLineColor(0);
////    text1 -> SetFillColor(0);
////    text1 -> AddText("d+Au 200 GeV");
//
//    TPaveText *text2 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
//    text2 -> SetTextSize(0.03);
//    text2 -> SetShadowColor(0);
//    text2 -> SetLineColor(0);
//    text2 -> SetFillColor(0);
//    text2->SetTextColor(28);
//    text2 -> AddText("cuts I.");
//
//    TPaveText *text3 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
//    text3 -> SetTextSize(0.03);
//    text3 -> SetShadowColor(0);
//    text3 -> SetLineColor(0);
//    text3 -> SetFillColor(0);
//    text3->SetTextColor(9);
//    text3 -> AddText("cuts II.");
////    TLine *leftline2 = new TLine(0.003,hEvents->GetMinimum(),0.003,hEvents->GetMaximum()+0.05*hEvents->GetMaximum());
////    leftline2->SetLineWidth(2);
////    leftline2->SetLineStyle(9);
////    leftline2->SetLineColor(8);
//
//    hEvents -> Draw();
//    legend -> Draw("same");
//    leftline2->Draw("same");
//    leftline3->Draw("same");
////    text1->Draw("same");
//    text2->Draw("same");
//    text3->Draw("same");
//}



//void plotEv() {
//
//    TString folder = "res_ntp/p17id/";
//    TString input = "res_ntp_dau_p17id_wide.root";
//
//    TFile* data = new TFile(folder + input ,"r");
//
//    TCanvas *c3 = new TCanvas("c3","c3",1000,1200);
//
//    hEvents = (TH1F*)data -> Get("hEventStat1");
////    hEvents ->GetYaxis()->SetTitleOffset(1.5);
////    hInvMassBackPlus = (TH1F*)data -> Get("background plus");
//    hEvents -> SetTitle("");
//    hEvents -> SetStats(0);
//    hEvents -> SetFillColor(46);
//    hEvents->SetLineWidth(2);
//    hEvents -> SetLineColor(46);
//    hEvents ->GetYaxis()-> SetTitleOffset(1.4);
////    hEvents ->GetXaxis()-> SetTitle("DCA daughters [cm]");
////    hEvents ->GetYaxis()-> SetTitle("Number of pairs");
////    hEvents -> SetFillStyle(3344);
//
//    gPad->SetLeftMargin(0.15);
//    hEvents -> Draw();
//
//}

