#include "/home/lukas/work/tmva_d0/BDT/efficiencyScripts/efficiencyCompare.cpp"

void plotTogetherEfficiencyCompareImagesMoreFiles();
void plotTogetherEfficiencyCompareImages(TString fileName1, TString fileName2, TString folder, TString legend1, TString legend2);

//__________________________________________________________________________________________________________________________________________
void plotTogetherEfficiencyCompareImages(TString fileName1, TString fileName2, TString folder, TString legend1, TString legend2) {
    folderDate=folder;

    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    TString inputPlot[]={fileName1, fileName2};
    TString legendNamesPlot[]={legend1, legend2};

    TString grNames[] = {"grTpcAccHftPidPreCutsBDT",
                         "grTpcAcc",
                         "grTpcAccHftPid",
                         "grTpcAccHft",
                         "grTpcAccHftPidPreCuts",
                         "grTpcAccHftPreCuts",
                         "grTpcAccHftPreCutsBDT",
                         "grAcc"};

    TString grNamesAxis[] = {"Acc+TPC+HFT+PID+BDT efficiency",
                             "Acc+TPC efficiency",
                             "Acc+TPC+HFT+PID efficiency",
                             "Acc+TPC+HFT efficiency",
                             "Acc+TPC+HFT+PID+precuts efficiency",
                             "Acc+TPC+HFT+precuts efficiency",
                             "Acc+TPC+HFT+BDT efficiency",
                             "Acc efficiency"};
    const int nGrNames = sizeof(grNames) / sizeof(TString);

    plot(inputPlot, 2, legendNamesPlot, grNames, nGrNames, grNamesAxis);
    TString folderOut = Form("finalAnalysis/efficiencyRatio/%s/%s_to_%s", folder.Data(), fileName2.Data(), fileName1.Data() );
    gSystem->Exec(Form("mkdir -p %s", folderOut.Data()));
    gSystem->Exec(Form("mv finalAnalysis/efficiencyRatio/*.png %s/.", folderOut.Data()));
    gSystem->Exec(Form("mv finalAnalysis/efficiencyRatio/*.eps %s/.", folderOut.Data()));
    gSystem->Exec(Form("mv finalAnalysis/efficiencyRatio/final_eff_ratio_SIM.root %s/.", folderOut.Data()));
}

//_________________________________________________________________________________________________________________________________________
void plotTogetherEfficiencyCompareImagesMoreFiles() {
    TString folder="narrowBins";
    folderDate=folder;

    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

//    TString inputPlot[] = {"ntp_FS_hijing_nonPrimary_DCA1cm_recoEvent_D0weightsHJ",
//                                "ntp_fullEvent_full_production.vtx.3M.1308.AllD0"};

    ///toto chcem:
//    TString inputPlot[] = {"ntp.FS.1cm.primaries.newE.1509.goodPid",
//                           "ntp.FS.1cm.all.newE.1509.goodPid",
//                           "ntp.FS.data.hft2tof2.global.3009.goodPid",
//                           "ntp.FS.hft2.global.1multEdge.0810.goodPid",
//                           "ntp.FS.hft2.global.4multEdge.0810.goodPid",
//                           "ntp_FS_data_primaries_HS_dca1_0302.goodPid_dca1"};

    TString inputPlot[] = {"ntp_FS_data_global_HS_dca1_0302.goodPid_dca1",
                           "ntp_FS_data_primaries_HS_dca1_0302.goodPid_dca1"};

    TString legendNamesPlot[] = {"global",
                                 "primaries"};

//    TString legendNamesPlot[] = {"primary, DCA < 1 cm",
//                                 "global, DCA < 1 cm",
//                                 "global, evts with N_{HFT}>2, N_{TOF}>2",
//                                 "global, evts with N_{HFT}>2, 1 multiplicity class in FS",
//                                 "global, evts with N_{HFT}>2, 4 multiplicity classes in FS",
//                                 "nn"};

    int nFiles= sizeof(inputPlot) / sizeof(TString);
    cout<<"number of input files: "<<nFiles<<endl;


    TString grNames[] = {"grTpcAccHftPidPreCutsBDT",
                         "grTpcAcc",
                         "grTpcAccHftPid",
                         "grTpcAccHft",
                         "grTpcAccHftPidPreCuts",
                         "grTpcAccHftPreCuts",
                         "grTpcAccHftPreCutsBDT",
                         "grAcc"};

    TString grNamesAxis[] = {"Acc+TPC+HFT+PID+BDT efficiency",
                             "Acc+TPC efficiency",
                             "Acc+TPC+HFT+PID efficiency",
                             "Acc+TPC+HFT efficiency",
                             "Acc+TPC+HFT+PID+precuts efficiency",
                             "Acc+TPC+HFT+precuts efficiency",
                             "Acc+TPC+HFT+BDT efficiency",
                             "Acc efficiency"};

    const int nGrNames = sizeof(grNames) / sizeof(TString);


    auto *time = new TDatime();
    Int_t day = time->GetDay();
    Int_t month = time->GetMonth();
    Int_t hour = time->GetHour();
    Int_t minute = time->GetMinute();
    TString folderDate = Form("%i%i_%i%i", month, day, hour, minute);
    cout << folderDate << endl;

    TString folderOut = Form("finalAnalysis/efficiencyRatio/%s/moreFiles/%s", folder.Data(), folderDate.Data());

    for (int i = 0; i < 1; ++i) {
        cout<<inputPlot[i]<<" "<<legendNamesPlot[i]<<endl;
    }
    plot(inputPlot, nFiles, legendNamesPlot, grNames, nGrNames, grNamesAxis);
    gSystem->Exec(Form("mkdir -p %s/eps", folderOut.Data()));
    gSystem->Exec(Form("mv finalAnalysis/efficiencyRatio/*.png %s/.", folderOut.Data()));
    gSystem->Exec(Form("mv finalAnalysis/efficiencyRatio/*.eps %s/eps/.", folderOut.Data()));
}