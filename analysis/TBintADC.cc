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
    TFile* pidFile = TFile::Open(("./auxPID/auxPID_Run_" + std::to_string(runNum) + ".root").c_str());
    
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

    // Get channel ID of pre-shower
    TBcid pscid = utility.getcid(TBdetector::detid::preshower);
    // Get pedestal of pre-shower
    float psPed = utility.retrievePed(pscid);
    // Get PS PID cut using (mean value of PS 1-mip fit) * 2.5
    TF1* psFitFunc = (TF1*) pidFile->Get("psFitFunc");
    float psPIDcut = psFitFunc->GetParameter(1) * 2.5f; // This is not the proper way to set cut!! Details will be presented in workshop

    // Exercise 1 : Get cid and pedestal of module 1 tower 1 C ans S channel
    TBcid M1T1C_cid = ; // Your answer here
    TBcid M1T1S_cid = ; // Your answer here
    float M1T1C_ped = ; // Your answer here
    float M1T1S_ped = ; // Your answer here

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
        
        // Exercise 2 : Get TBwaveform data and std::vector<short> waveform of module 1 tower 1 C and S channels
        TBwaveform M1T1C_data = ; // Your answer here
                                ; // Your answer here
        std::vector<short> M1T1C_waveform = ; // Your answer here
                                            ; // Your answer here
        
        // PS Int. ADC
        float psIntADC = 0.f;
        for (int bin = 220; bin < 390; bin++) {
            int waveformBin = bin + 1;
            psIntADC += psPed - psWaveform[waveformBin];
        }
        //---------------------------------------------------------------------------------------------------------------------
        // Exercise 3 : Get integrated ADC of module 1 tower 1 C and S channels
        //              Look into avgTimeStructure root output and determine appropriate integration range for C and S channels
        //---------------------------------------------------------------------------------------------------------------------
        // M1 C Int. ADC
        float M1T1C_IntADC = 0.f;
        for (int bin = ; bin < ; bin++) {
            int waveformBin = bin + 1;
            M1T1C_IntADC += ;
        }
        // M1 S Int. ADC
        float M1T1S_IntADC = 0.f;
        //---------------------------------------------------------------------------------------------------------------------

        // PS and module Int. ADC plot before PID
        psIntHist->Fill(psIntADC);
        M1T1C_Hist->Fill(M1T1C_IntADC);
        M1T1S_Hist->Fill(M1T1S_IntADC);

        // For DWC corr cut and Pre-shower PID cut
        if (! dwcCorrelationPID(dwc1_correctedPosition, dwc2_correctedPosition, 1.5f) ) continue;
        if (psIntADC < psPIDcut) continue;

        // PS and module Int. ADC plot after PID
        psIntHist_PID->Fill(psIntADC);
        M1T1C_Hist_PID->Fill(M1T1C_IntADC);
        M1T1S_Hist_PID->Fill(M1T1S_IntADC);
    }

    TCanvas* c = new TCanvas();

    std::string outFile = "./intADC/intADC_Run_" + std::to_string(runNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    psIntHist->Write();
    psIntHist->Draw();
    c->Update();
    float psHist_min = gPad->GetUymin();
    float psHist_max = gPad->GetUymax();
    psIntHist_PID->GetYaxis()->SetRangeUser(psHist_min, psHist_max);
    psIntHist_PID->Write();

    M1T1C_Hist->Write();
    M1T1S_Hist->Write();

    M1T1C_Hist->Draw();
    c->Update();
    float M1T1C_min = gPad->GetUymin();
    float M1T1C_max = gPad->GetUymax();
    M1T1C_Hist_PID->GetYaxis()->SetRangeUser(M1T1C_min, M1T1C_max);
    M1T1C_Hist_PID->Write();

    M1T1S_Hist->Draw();
    c->Update();
    float M1T1S_min = gPad->GetUymin();
    float M1T1S_max = gPad->GetUymax();
    M1T1S_Hist_PID->GetYaxis()->SetRangeUser(M1T1S_min, M1T1S_max);
    M1T1S_Hist_PID->Write();

    outputRoot->Close();
}
