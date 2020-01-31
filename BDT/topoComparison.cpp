#include "topoComparison.h"


#include "TCut.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPad.h"
#include "TFile.h"
#include <iostream>

ClassImp(topoComparison)



topoComparison::topoComparison() :
                   mCuts(), mVars(), mTuples(){
    cout<<"ahoj"<<endl;
    mOutputFileName="out.root";

    auto* time=new TDatime();
    Int_t day=time->GetDay();
    Int_t month=time->GetMonth();
    Int_t hour=time->GetHour();
    Int_t minute=time->GetMinute();
    TString folderDate=Form("%i%i_%i%i", month, day, hour, minute);
    cout<<folderDate<<endl;
    mFolder=Form("topoProjections/%s/", folderDate.Data());

    gSystem->Exec(Form("mkdir topoProjections/%s", folderDate.Data()));
//    gSystem->Exec(Form("mkdir finalAnalysis/%s/img", folderDate.Data()));
//    gSystem->Exec(Form("mkdir finalAnalysis/%s/img/root", folderDate.Data()));
}

//____________________________________________________________________________________________
void topoComparison::project() {
    TCut cuts = connectCuts();
//    cout<<"Cuts you set are: "<<cuts<<endl;

//    std::vector<TH1F*> mHists;

    TFile* outF = new TFile(mFolder+mOutputFileName ,"RECREATE");

    const int nTuples = mTuples.size();
    TString nameHisto;
    TH1F* mHists[nTuples];

    if (mVars.size()==0 || mTuples.size()==0){
        cout<<"no vars or tuples to project"<<endl;
        return;
    }

    for(unsigned int iVar = 0; iVar < mVars.size(); ++iVar) {
        auto *c = new TCanvas("c","c",900,900);
        auto *legend = new TLegend(0.6, 0.81, 0.75, 0.89);
        gPad->SetLeftMargin(0.15);

        for(unsigned int iNt = 0; iNt < mTuples.size(); ++iNt) {
            nameHisto = mVars[iVar] + "_" + mTuples[iNt]->GetName();
            cout << nameHisto << endl;
            mHists[iNt] = new TH1F(nameHisto, nameHisto, mnBins[iVar], mLimsMin[iVar], mLimsMax[iVar]);
            mHists[iNt]->SetNameTitle(nameHisto, "");
            mHists[iNt]->GetXaxis()->SetTitle(mTitleX[iVar]);
            mHists[iNt]->GetYaxis()->SetTitle("1/N_{entries}");
            mTuples[iNt]->Project(nameHisto, mVars[iVar], cuts);
            setHisto(mHists[iNt], iNt);
            mHists[iNt]->Write();
            mHists[iNt]->Scale(1 / mHists[iNt]->GetEntries());
            legend->AddEntry(mHists[iNt], mTuples[iNt]->GetName(), "pl");

        }
//        float maximumY = findMaximumHistos(*mHists, mTuples.size());

        float max = 0;
        for (unsigned int i1 = 0; i1 < mTuples.size(); ++i1) {
            if (mHists[i1]->GetMaximum()>max) max = mHists[i1]->GetMaximum();
        }

        for (unsigned int i = 0; i < mTuples.size(); ++i) {
            mHists[i]->GetYaxis()->SetRangeUser(0, max*1.1);
            if (i==0) mHists[i]->Draw(); else mHists[i]->Draw("same");
        }
        legend->SetFillStyle(0);
        legend->SetLineColor(0);
        legend->SetTextSize(0.035);
        legend->Draw("same");
        c->SaveAs(mFolder+mVars[iVar]+".png");
        c->Clear();
        c->Close();
        delete c;
        delete legend;
    }
    outF->Close();
}

//____________________________________________________________________________________________
void topoComparison::setHisto(TH1F* histo, int iterator){
    Int_t colors[] = {1, 46, 8, 9};
    histo->SetMarkerStyle(20);
    histo->Sumw2();
    histo->SetStats(0);
    histo->SetMarkerColor(colors[iterator]);
    histo->SetLineColor(colors[iterator]);
}

//___________________________________________________________________________________________
float topoComparison::findMaximumHistos(TH1F* hists, int nhists){
//    int nHists = sizeof(hists) / sizeof(TH1F);
    float max = 0;
    for (int i = 0; i < nhists; ++i) {
        if (hists[i].GetMaximum()>max) max = hists[i].GetMaximum();
    }
    return max;
}

//____________________________________________________________________________________________
TCut topoComparison::connectCuts() {
    TCut setCuts = "";
    for(unsigned int k = 0; k < mCuts.size(); ++k) {
        setCuts += mCuts[k];
    }
    return setCuts;
}