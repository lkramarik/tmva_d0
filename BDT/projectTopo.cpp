#include "topoComparison.h"
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLine.h"
#include "TList.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TSystem.h"
#include "TROOT.h"


void projectTopo() {
    gROOT->ProcessLine(".L topoComparison.cpp++");

    auto* dataF = new TFile("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.0110.root" ,"r");
    auto* simF = new TFile("/home/lukas/work/tmva_d0/sim/ntpTMVA_D0.toyMC.0910.fullEff.root" ,"r");
    auto* ntpData = (TNtuple*)dataF -> Get("ntp_signal");
    auto* ntpSim = (TNtuple*)simF -> Get("ntp_signal");

    topoComparison topo;
    topo.setOutFileName("topologyProjections.root");
    topo.addNtpToCompare(ntpData, "data");
    topo.addNtpToCompare(ntpSim, "sim");
    topo.addCut("D_mass<3");
    topo.addCut("D_mass>1.7");
    topo.addVarToCompare("k_dca", 0, 0.5, "DCA_{K} [cm]");
    topo.addVarToCompare("pi1_dca", 0, 0.3, "DCA_{#pi} [cm]");

    topo.project();

    dataF->Close();
    simF->Close();

}
