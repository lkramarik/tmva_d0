//
// Created by lukas on 10. 10. 2019.
//
TF1* f1D0pp;
TF1* fD0scaledAuAu;
TF1 *fD0dAu;
TF1 *fD0scaledpp;

double ratioFct(double *x, double *par){
//    return fD0scaledAuAu->Eval(x[0])/fD0scaledpp->Eval(x[0]);
    return fD0dAu->Eval(x[0])/fD0scaledpp->Eval(x[0]);
}

double scaleFct(double *x, double *par){
//    return fD0dAu->Eval(x[0])/fD0scaledAuAu->Eval(x[0]);
    return 8*f1D0pp->Eval(x[0])/50;
}


void D0invYield() {
    Double_t nColl = (10.48+16.11+24.59+36.13)/4;
    cout<<"n coll is "<<nColl;
    nColl = 21.37; //from paper
    Double_t pT[] = {0.242, 0.728, 1.22, 1.71, 2.21, 2.71, 3.36, 4.36, 5.37, 6.59};
    Double_t invYield[] = {3.81e-3, 2.72e-03, 1.10e-03, 5.58e-04, 2.52e-04, 1.16e-04, 2.98e-05, 5.28e-06, 2.07e-06, 3.38e-07};
    const int nBins = sizeof(pT) / sizeof(Double_t);

    TGraph *gr = new TGraph(nBins, pT, invYield);
//    TGraph *gr = new TGraph(nBins, pT, y);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(9);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitle("Inv. yield D^{0} 2016 60-80%");
    gr->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    TF1 *fD0 = new TF1("fD0", "(1+x)/(2*pi)*[dndy]*([A]-1)*([A]-2)/([A]*[B]*([A]*[B]+1.864*([A]-2)))*pow(1+(sqrt(x*x+1.864*1.864)-1.864)/([A]*[B]),-[A])", 0.2, 7);
    fD0->SetParLimits(0, 0, 10);
    fD0->SetParLimits(1, 0, 1);
    fD0->SetParLimits(2, 0, 0.1);
    fD0->SetParameter(0, 6.5);
    fD0->SetParameter(1, 0.17);
    fD0->SetParameter(2, 0.009);

    gr->Fit("fD0","RI");


//    gr->Draw("ap");

    TFile *fData = new TFile("D0spectraAuAu.root","RECREATE");
    fD0->Write("f1D0spectraAuAu6080");
    gr->Write("grD0spectraAuAu6080");

    fD0scaledAuAu = new TF1("fD0scaledAuAu", "fD0(x)*10/21.37", 0.2, 7);
//    fD0scaledAuAu->Draw("same");
    fD0scaledAuAu->Write("f1D0scaledAuAu6080");
    fData->Close();

    //dAu
    Double_t invYieldPubl[] = {0.00678, 0.00329, 0.00216, 0.000217};
    Double_t invYieldPublErr[] = {0.00193, 0.00102, 0.00055, 0.000094};
    Double_t ptPubl[] = {0.3, 0.75, 1.25, 2.25};
    Double_t ptPublE[] = {0.3, 0.75, 1.25, 2.25};

    TGraphErrors* grInvYieldsPubl = new TGraphErrors(4, ptPubl, invYieldPubl, 0, invYieldPublErr);
    grInvYieldsPubl->SetMarkerStyle(20);
    grInvYieldsPubl->GetXaxis()->SetLimits(0,6);

    Double_t xbins[] = {0, 0.1, 0.5, 1, 1.5, 3, 4.5, 6};
    TH1D* hdAuPubl=new TH1D("hdAuPubl", "hdAuPubl", 10, xbins );
    hdAuPubl->SetStats(0);
    hdAuPubl->GetXaxis()->SetRangeUser(0., 6.0);
    hdAuPubl->GetYaxis()->SetRangeUser(0., 0.01);

    for (int j = 2; j < 6; ++j) {
        hdAuPubl->SetBinContent(j, invYieldPubl[j-1]);
        hdAuPubl->SetBinError(j, invYieldPublErr[j-1]);
    }

    fD0dAu = new TF1("fD0dAu", "1/(2*pi)*[dndy]*([A]-1)*([A]-2)/([A]*[B]*([A]*[B]+1.864*([A]-2)))*pow(1+(sqrt(x*x+1.864*1.864)-1.864)/([A]*[B]),-[A])", 0.25, 6);
    fD0dAu->SetParLimits(0, -100000, -5000);
    fD0dAu->SetParLimits(1, -1, 1);
    fD0dAu->SetParameter(1, 0.3);
    fD0dAu->SetParLimits(2, 0.001, 0.05);
    fD0dAu->SetParameter(2, 0.028);

    grInvYieldsPubl->Fit(fD0dAu, "", "", 0.25, 3);
//    hdAuPubl->Fit(fD0dAu, "", "", 0.25, 6);

    TH1D *hint = new TH1D("hint",
                          "Fitted gaussian with .95 conf.band", 100, 0, 6);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
    //Now the "hint" histogram has the fitted function values as the
    //bin contents and the confidence intervals as bin errors
    hint->SetStats(kFALSE);
    hint->SetFillColor(17);
    hint->SetTitle(0);
    hint->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hint->GetYaxis()->SetTitle("D^{0} dAu 2004");

    TFile *fDataDau = new TFile("D0spectra_dAu.root","RECREATE");
    fD0dAu->Write("f1D0spectra_dAu");
    grInvYieldsPubl->Write("grD0spectra_dAu");
    fDataDau->Close();


    TCanvas *cSign = new TCanvas("cSign","cSign",900,900);
    gPad->SetLogy();

    hint->Draw("e3 same");
    grInvYieldsPubl->Draw("p same");
    cSign->SaveAs("fitdAuLevy.png");
    cSign->SaveAs("fitdAuLevy.pdf");

    Double_t ratioFit[4], ratioFitError[4];
    for (int j = 0; j < 4; ++j) {
        ratioFit[j]=invYieldPubl[j]/fD0dAu->Eval(ptPubl[j]);
        ratioFitError[j]=invYieldPublErr[j]/invYieldPubl[j];
    }

    TCanvas *cRatioFit = new TCanvas("cRatioFit","cRatioFit",900,900);
    TGraphErrors* grInvYieldsPublRatioFit = new TGraphErrors(4, ptPubl, ratioFit, 0, invYieldPublErr);
    grInvYieldsPublRatioFit->SetTitle("");
    grInvYieldsPublRatioFit->SetMarkerStyle(20);
    grInvYieldsPublRatioFit->Draw("ap");


    //pp
    TFile *fppSys = new TFile("errorEstimate/out_ppsys.root","READ");
    f1D0pp = (TF1*) fppSys->Get("Levynew_pp");
//    cout<<f1D0pp->GetName()<<endl;
    f1D0pp->SetName("f1D0pp");
    fD0scaledpp = new TF1("fD0scaledpp", scaleFct, 0.2, 7);


    cout<<fD0dAu->Eval(3)/f1D0pp->Eval(3)<<endl;
    cout<<fD0scaledAuAu->Eval(3)<<endl;
    cout<<f1D0pp->Eval(0.908)<<endl;
    cout<<f1D0pp->Eval(1.57)<<endl;
    cout<<f1D0pp->Eval(2.45)<<endl;

    cout<<fD0scaledpp->Eval(0.908)<<endl;
    cout<<fD0scaledpp->Eval(1.57)<<endl;
    cout<<fD0scaledpp->Eval(2.45)<<endl;

    Double_t raaDau[4], raaDauError[4];
    for (int i = 0; i < 4; ++i) {
        raaDau[i]=invYieldPubl[i]/fD0scaledpp->Eval(ptPubl[i]);
        raaDauError[i]=raaDau[i]*invYieldPublErr[i]/invYieldPubl[i];
    }

    TCanvas *c2 = new TCanvas("c2","c2",900,900);

    TGraphErrors* grRdAu = new TGraphErrors(4, ptPubl, raaDau, 0, raaDauError);
    grRdAu->SetTitle("");
    grRdAu->GetXaxis()->SetTitle("p_{T} GeV/c");
    grRdAu->GetYaxis()->SetTitle("R_{dAu}");
    grRdAu->GetXaxis()->SetLimits(0,6);
    grRdAu->GetYaxis()->SetRangeUser(0,2);
    grRdAu->SetMarkerStyle(20);

    TLegend *legend = new TLegend(0.492, 0.747, 0.643, 0.89);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    legend -> AddEntry(grRdAu, "D^{0} dAu 2004 / 8x pp Levy fit", "p");

    TF1* fRatio = new TF1("fRatio", ratioFct, 0.15, 7);
    fRatio->SetTitle("");
//    legend -> AddEntry(fRatio, "dAu Levy fit / pp Levy fit", "l");

    grRdAu->Draw("ap");
//    fRatio->Draw("same");
    legend->Draw("same");
    c2->SaveAs("rdAu_oldData.png");
    c2->SaveAs("rdAu_oldData.pdf");



}
