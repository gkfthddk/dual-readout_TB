#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "TBevt.h"
#include "TButility.h"
#include "functions.cc"

#include "TROOT.h"
#include "TStyle.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

// This macro will draw integrated ADC of module and pre-shower detector before and after the DWC & PS cut
// This macro will store 1D histograms in .root format

// How to execute
// On local or batch, run the following command :
// ./TBintADC.exe <run number> <max # of events to process>

int main(int argc, char** argv) {
    gStyle->SetOptFit(1);

    int runNum = std::stoi(argv[1]);
    int maxEntry = std::stoi(argv[2]);

    // Load ntuples using TChain. Need to change "getNtupleChain()" in function.cc to set proper path to your ntuples
    TChain* evtChain = getNtupleChain(runNum);
    TBevt<TBwaveform>* anEvt = new TBevt<TBwaveform>(); 
    evtChain->SetBranchAddress("TBevt", &anEvt);

    // Check total # of events to process
    int totalEntry = evtChain->GetEntries();
    std::cout << "Total entries : " << totalEntry << std::endl;
    if ( (maxEntry > 0) && (maxEntry < totalEntry) ) {
        totalEntry = maxEntry;
        std::cout << "Will process maximum " << std::to_string(totalEntry) << " events" << std::endl;
    }

    // Load mapping & pedestal info
    TButility utility = TButility();
    utility.loading("/gatbawi/dream/mapping/mapping_Aug2022TB.root");
    utility.loadped( ("/gatbawi/dream/ped/mean/Run" + std::to_string(runNum) + "_pedestalHist_mean.root").c_str() );
    
    // Get PID info
    TFile* pidFile = TFile::Open(("./box/auxPID_Run_" + std::to_string(runNum) + ".root").c_str());
    
    // Get DWC offset for DWC corr cut
    std::vector<TBcid> dwc1_cid;
    dwc1_cid.push_back(TBcid(1,17)); // R
    dwc1_cid.push_back(TBcid(1,19)); // L
    dwc1_cid.push_back(TBcid(1,21)); // U
    dwc1_cid.push_back(TBcid(1,23)); // D
    std::vector<TBcid> dwc2_cid;
    dwc2_cid.push_back(TBcid(1,25)); // R
    dwc2_cid.push_back(TBcid(1,27)); // L
    dwc2_cid.push_back(TBcid(1,29)); // U
    dwc2_cid.push_back(TBcid(1,31)); // D
    TH2D* dwc1_pos   = (TH2D*) pidFile->Get("dwc1_pos");
    TH2D* dwc2_pos   = (TH2D*) pidFile->Get("dwc2_pos");
    std::vector<float> dwc1_offset = getDWCoffset(dwc1_pos);
    std::vector<float> dwc2_offset = getDWCoffset(dwc2_pos);

    // Preparing histograms
    TH1F* psIntHist = new TH1F("psIntADC", "psIntADC;IntADC;evt", 1000, -5000., 200000.);
    TH1F* psIntHist_PID = new TH1F("psIntADC_PID", "psIntADC_PID;IntADC;evt", 1000, -5000., 200000.);
    // M1 C
    TH1F* M1T1C_Hist = new TH1F("M1T1C", "M1T1C;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T1C_Hist_PID = new TH1F("M1T1C_PID", "M1T1C_PID;IntADC;evt", 1000, -10000, 50000);
    // M1 S
    TH1F* M1T1S_Hist = new TH1F("M1T1S", "M1T1S;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T1S_Hist_PID = new TH1F("M1T1S_PID", "M1T1S_PID;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T1C_Hist_bad = new TH1F("M1T1C_bad", "M1T1C_bad;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T1S_Hist_bad = new TH1F("M1T1S_bad", "M1T1S_bad;IntADC;evt", 1000, -10000, 50000);
    // M1 C
    TH1F* M1T2C_Hist = new TH1F("M1T2C", "M1T2C;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T2C_Hist_PID = new TH1F("M1T2C_PID", "M1T2C_PID;IntADC;evt", 1000, -10000, 50000);
    // M1 S
    TH1F* M1T2S_Hist = new TH1F("M1T2S", "M1T2S;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T2S_Hist_PID = new TH1F("M1T2S_PID", "M1T2S_PID;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T2C_Hist_bad = new TH1F("M1T2C_bad", "M1T2C_bad;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T2S_Hist_bad = new TH1F("M1T2S_bad", "M1T2S_bad;IntADC;evt", 1000, -10000, 50000);
    // M1 C
    TH1F* M1T3C_Hist = new TH1F("M1T3C", "M1T3C;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T3C_Hist_PID = new TH1F("M1T3C_PID", "M1T3C_PID;IntADC;evt", 1000, -10000, 50000);
    // M1 S
    TH1F* M1T3S_Hist = new TH1F("M1T3S", "M1T3S;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T3S_Hist_PID = new TH1F("M1T3S_PID", "M1T3S_PID;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T3C_Hist_bad = new TH1F("M1T3C_bad", "M1T3C_bad;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T3S_Hist_bad = new TH1F("M1T3S_bad", "M1T3S_bad;IntADC;evt", 1000, -10000, 50000);
    // M1 C
    TH1F* M1T4C_Hist = new TH1F("M1T4C", "M1T4C;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T4C_Hist_PID = new TH1F("M1T4C_PID", "M1T4C_PID;IntADC;evt", 1000, -10000, 50000);
    // M1 S
    TH1F* M1T4S_Hist = new TH1F("M1T4S", "M1T4S;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T4S_Hist_PID = new TH1F("M1T4S_PID", "M1T4S_PID;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T4C_Hist_bad = new TH1F("M1T4C_bad", "M1T4C_bad;IntADC;evt", 1000, -10000, 50000);
    TH1F* M1T4S_Hist_bad = new TH1F("M1T4S_bad", "M1T4S_bad;IntADC;evt", 1000, -10000, 50000);

    TH1F* psIntHist_bad = new TH1F("psIntADC_bad", "psIntADC_bad;IntADC;evt", 1000, -5000., 200000.);

    TH2D* dwc1pos_bad = new TH2D("dwc1pos_bad","dwc1pos_bad",100,-50,50,100,-50,50);
    TH2D* dwc2pos_bad = new TH2D("dwc2pos_bad","dwc2pos_bad",100,-50,50,100,-50,50);
    TH2D* dwcdx_bad = new TH2D("dwcdx_bad","dwcdx_bad",100,-50,50,100,-50,50);

    // Get channel ID of pre-shower
    TBcid pscid = utility.getcid(TBdetector::detid::preshower);
    // Get pedestal of pre-shower
    float psPed = utility.retrievePed(pscid);
    // Get PS PID cut using (mean value of PS 1-mip fit) * 2.5
    TF1* psFitFunc = (TF1*) pidFile->Get("psFitFunc");
    float psPIDcut = psFitFunc->GetParameter(1); // This is not the proper way to set cut!! Details will be presented in workshop

    // Exercise 1 : Get cid and pedestal of module 1 tower 1 C ans S channel
    TBcid M1T1C_cid = utility.getcid(TBdetector::detid::PMT,1,4,true); // Your answer here
    TBcid M1T1S_cid = utility.getcid(TBdetector::detid::PMT,1,4,false); // Your answer here
    float M1T1C_ped = utility.retrievePed(M1T1C_cid); // Your answer here
    float M1T1S_ped = utility.retrievePed(M1T1S_cid); // Your answer here
    TBcid M1T2C_cid = utility.getcid(TBdetector::detid::PMT,1,4,true); // Your answer here
    TBcid M1T2S_cid = utility.getcid(TBdetector::detid::PMT,1,4,false); // Your answer here
    float M1T2C_ped = utility.retrievePed(M1T2C_cid); // Your answer here
    float M1T2S_ped = utility.retrievePed(M1T2S_cid); // Your answer here
    TBcid M1T3C_cid = utility.getcid(TBdetector::detid::PMT,1,4,true); // Your answer here
    TBcid M1T3S_cid = utility.getcid(TBdetector::detid::PMT,1,4,false); // Your answer here
    float M1T3C_ped = utility.retrievePed(M1T3C_cid); // Your answer here
    float M1T3S_ped = utility.retrievePed(M1T3S_cid); // Your answer here
    TBcid M1T4C_cid = utility.getcid(TBdetector::detid::PMT,1,4,true); // Your answer here
    TBcid M1T4S_cid = utility.getcid(TBdetector::detid::PMT,1,4,false); // Your answer here
    float M1T4C_ped = utility.retrievePed(M1T4C_cid); // Your answer here
    float M1T4S_ped = utility.retrievePed(M1T4S_cid); // Your answer here
    TBcid mu_cid = utility.getcid(TBdetector::detid::muon);
    float mu_ped = utility.retrievePed(mu_cid);
    
    std::string outFile = "./box/Tree_Run_" + std::to_string(runNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    TTree * event = new TTree("event","event info");
    float psIntADC = 0.f;
    float muIntADC = 0.f;
        float M1T1C_IntADC = 0.f;
        float M1T2C_IntADC = 0.f;
        float M1T3C_IntADC = 0.f;
        float M1T4C_IntADC = 0.f;
        float M1T1S_IntADC = 0.f;
        float M1T2S_IntADC = 0.f;
        float M1T3S_IntADC = 0.f;
        float M1T4S_IntADC = 0.f;

    float dwc1_H=0.;
    float dwc1_V=0.;
    float dwc2_H=0.;
    float dwc2_V=0.;
    bool dwcPID=false;
    bool PID=false;
    event->Branch("psIntADC",&psIntADC,"psIntADC/F");
    event->Branch("psPIDcut",&psPIDcut,"psPIDcut/F");
    event->Branch("muIntADC",&muIntADC,"muIntADC/F");
    event->Branch("M1T1C_IntADC",&M1T1C_IntADC,"M1T1C_IntADC/F");
    event->Branch("M1T1S_IntADC",&M1T1S_IntADC,"M1T1S_IntADC/F");
    event->Branch("M1T2C_IntADC",&M1T2C_IntADC,"M1T2C_IntADC/F");
    event->Branch("M1T2S_IntADC",&M1T2S_IntADC,"M1T2S_IntADC/F");
    event->Branch("M1T3C_IntADC",&M1T3C_IntADC,"M1T3C_IntADC/F");
    event->Branch("M1T3S_IntADC",&M1T3S_IntADC,"M1T3S_IntADC/F");
    event->Branch("M1T4C_IntADC",&M1T4C_IntADC,"M1T4C_IntADC/F");
    event->Branch("M1T4S_IntADC",&M1T4S_IntADC,"M1T4S_IntADC/F");
    event->Branch("dwc1_H",&dwc1_H,"dwc1_H/F");
    event->Branch("dwc1_V",&dwc1_V,"dwc1_V/F");
    event->Branch("dwc2_H",&dwc2_H,"dwc2_H/F");
    event->Branch("dwc2_V",&dwc2_V,"dwc2_V/F");
    event->Branch("dwcPID",&dwcPID,"dwcPID/B");
    event->Branch("PID",&PID,"PID/B");

    // Evt Loop
    for (int iEvt = 0; iEvt < totalEntry; iEvt++) {
        printProgress(iEvt + 1, totalEntry);

        evtChain->GetEntry(iEvt);

        // Aux data & waveform
        TBwaveform psData = anEvt->data(pscid);
        std::vector<short> psWaveform = psData.waveform();

        // DWC data & waveform
        std::vector<TBwaveform> dwc1_data;
        for (TBcid cid : dwc1_cid) {
            dwc1_data.push_back(anEvt->data(cid));
        }
        std::vector<TBwaveform> dwc2_data;
        for (TBcid cid : dwc2_cid) {
            dwc2_data.push_back(anEvt->data(cid));
        }
        std::vector< std::vector<short> > dwc1_waveform;
        for (TBwaveform data : dwc1_data) {
            dwc1_waveform.push_back(data.waveform());
        }
        std::vector< std::vector<short> > dwc2_waveform;
        for (TBwaveform data : dwc2_data) {
            dwc2_waveform.push_back(data.waveform());
        }
        std::vector<float> dwc1_peakTime;
        for (std::vector<short> waveform : dwc1_waveform) {
            dwc1_peakTime.push_back(getPeakTime(waveform));
        }
        std::vector<float> dwc2_peakTime;
        for (std::vector<short> waveform : dwc2_waveform) {
            dwc2_peakTime.push_back(getPeakTime(waveform));
        }
        std::vector<float> dwc1_correctedPosition = getDWC1position(dwc1_peakTime, dwc1_offset.at(0), dwc1_offset.at(1));
        std::vector<float> dwc2_correctedPosition = getDWC2position(dwc2_peakTime, dwc2_offset.at(0), dwc2_offset.at(1));
        dwc1_H=dwc1_correctedPosition.at(0);
        dwc1_V=dwc1_correctedPosition.at(1);
        dwc2_H=dwc2_correctedPosition.at(0);
        dwc2_V=dwc2_correctedPosition.at(1);
        
        // Exercise 2 : Get TBwaveform data and std::vector<short> waveform of module 1 tower 1 C and S channels
        TBwaveform M1T1C_data = anEvt->data(M1T1C_cid); // Your answer here
        TBwaveform M1T1S_data = anEvt->data(M1T1S_cid); // Your answer here
        std::vector<short> M1T1C_waveform = M1T1C_data.waveform(); // Your answer here
        std::vector<short> M1T1S_waveform = M1T1S_data.waveform(); // Your answer here
        TBwaveform M1T2C_data = anEvt->data(M1T2C_cid); // Your answer here
        TBwaveform M1T2S_data = anEvt->data(M1T2S_cid); // Your answer here
        std::vector<short> M1T2C_waveform = M1T2C_data.waveform(); // Your answer here
        std::vector<short> M1T2S_waveform = M1T2S_data.waveform(); // Your answer here
        TBwaveform M1T3C_data = anEvt->data(M1T3C_cid); // Your answer here
        TBwaveform M1T3S_data = anEvt->data(M1T3S_cid); // Your answer here
        std::vector<short> M1T3C_waveform = M1T3C_data.waveform(); // Your answer here
        std::vector<short> M1T3S_waveform = M1T3S_data.waveform(); // Your answer here
        TBwaveform M1T4C_data = anEvt->data(M1T4C_cid); // Your answer here
        TBwaveform M1T4S_data = anEvt->data(M1T4S_cid); // Your answer here
        std::vector<short> M1T4C_waveform = M1T4C_data.waveform(); // Your answer here
        std::vector<short> M1T4S_waveform = M1T4S_data.waveform(); // Your answer here

        TBwaveform mu_data = anEvt->data(mu_cid);
        std::vector<short> mu_waveform = mu_data.waveform();
        
        // PS Int. ADC
        psIntADC = 0.f;
        for (int bin = 220; bin < 390; bin++) {
            int waveformBin = bin + 1;
            psIntADC += psPed - psWaveform[waveformBin];
        }
        //---------------------------------------------------------------------------------------------------------------------
        // Exercise 3 : Get integrated ADC of module 1 tower 1 C and S channels
        //              Look into avgTimeStructure root output and determine appropriate integration range for C and S channels
        //---------------------------------------------------------------------------------------------------------------------
        // M1 C Int. ADC
        M1T1C_IntADC = 0.f;
        M1T2C_IntADC = 0.f;
        M1T3C_IntADC = 0.f;
        M1T4C_IntADC = 0.f;
        for (int bin = 120; bin < 240; bin++) {
            int waveformBin = bin + 1;
            M1T1C_IntADC += M1T1C_ped - M1T1C_waveform[waveformBin];
            M1T2C_IntADC += M1T2C_ped - M1T2C_waveform[waveformBin];
            M1T3C_IntADC += M1T3C_ped - M1T3C_waveform[waveformBin];
            M1T4C_IntADC += M1T4C_ped - M1T4C_waveform[waveformBin];
        }
        // M1 S Int. ADC
        M1T1S_IntADC = 0.f;
        M1T2S_IntADC = 0.f;
        M1T3S_IntADC = 0.f;
        M1T4S_IntADC = 0.f;
        for (int bin = 120; bin < 240; bin++) {
            int waveformBin = bin + 1;
            M1T1S_IntADC += M1T1S_ped - M1T1S_waveform[waveformBin];
            M1T2S_IntADC += M1T2S_ped - M1T2S_waveform[waveformBin];
            M1T3S_IntADC += M1T3S_ped - M1T3S_waveform[waveformBin];
            M1T4S_IntADC += M1T4S_ped - M1T4S_waveform[waveformBin];
        }
        muIntADC = 0.f;
        for (int bin = 840; bin < 990; bin++) {
            int waveformBin = bin+1;
            muIntADC += mu_ped - mu_waveform[waveformBin];
        }
        //---------------------------------------------------------------------------------------------------------------------

        // PS and module Int. ADC plot before PID

        // For DWC corr cut and Pre-shower PID cut
        
        dwcPID=dwcCorrelationPID(dwc1_correctedPosition, dwc2_correctedPosition, 1.5f);
        PID=true;
        if (psIntADC < 4.5*psPIDcut) PID=false;
        if (muIntADC > 0) PID=false;
        event->Fill();

    }



    outputRoot->Write();
    outputRoot->Close();
}