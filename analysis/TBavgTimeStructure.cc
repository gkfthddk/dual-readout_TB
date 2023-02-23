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

// --------------------------------------------------------------------------------------------------------------------
// Exercise : All info about drawing pre-shower detectors's average time structure is already filled in this script
//            Now, do the same thing for Module 1 Tower 1 S and C channel to draw average time structure of them
//            Store all pre-shower detector, module 1 Tower 1 S and C average time structure in single root file
// --------------------------------------------------------------------------------------------------------------------

// This macro will draw average time structure (== mean of waveforms) of pre-shower, module 1 tower 1 S & C channel from waveform ntuple
// This macro will store 1D histograms in .root format

// How to execute
// On local or batch, run the following command :
// ./TBavgTimeStructure.exe <run number> <max # of events to process>

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

    // Load mapping info
    TButility utility = TButility();
    utility.loading("/gatbawi/dream/mapping/mapping_Aug2022TB.root");
    
    // Preparing histograms
    TH1F* psHist = new TH1F("psAvgTime", "psAvgTime;bin;ADC", 1000, 0, 1000);
    TH1F* M1T1CHist = new TH1F("M1T1CAvgTime", "M1T1CAvgTime;bin;ADC", 1000, 0, 1000);
    TH1F* M1T1SHist = new TH1F("M1T1SAvgTime", "M1T1SAvgTime", 1000, 0, 1000);
    TH2D* dwc1_correctedPos = new TH2D("dwc1_correctedPos", "dwc1_correctedPos;mm;mm;events", 480, -120., 120., 480, -120., 120.);
    TH2D* dwc2_correctedPos = new TH2D("dwc2_correctedPos", "dwc2_correctedPos;mm;mm;events", 480, -120., 120., 480, -120., 120.);
    TH2D* dwc1_correctedPos_PID = new TH2D("dwc1_correctedPos_PID", "dwc1_correctedPos_PID;mm;mm;events", 480, -120., 120., 480, -120., 120.);
    TH2D* dwc2_correctedPos_PID = new TH2D("dwc2_correctedPos_PID", "dwc2_correctedPos_PID;mm;mm;events", 480, -120., 120., 480, -120., 120.);

    /*
    In this exercise, we need info from DWC (For giving correlation cut), pre-shower and M1T1 S, C channels (for drawing average time structure)
    We first get cids for all channels needed
    */

    // Get channel IDs
    // For convenience, we put all DWC cids in DWC cid vector. Order is [R, L, U, D]
    std::vector<TBcid> dwc1_cid;
    dwc1_cid.push_back(TBcid(1,17)); // R
    dwc1_cid.push_back(TBcid(1,19)); // L
    dwc1_cid.push_back(TBcid(1,21)); // U
    dwc1_cid.push_back(TBcid(1,23)); // D
    // Exercise 1 : Get vector of DWC2 with same order as DWC1 using TBcid()
    std::vector<TBcid> dwc2_cid;
    dwc2_cid.push_back(TBcid(1,25)); // R
    dwc2_cid.push_back(TBcid(1,27)); // L
    dwc2_cid.push_back(TBcid(1,29)); // U
    dwc2_cid.push_back(TBcid(1,31)); // D

    // Exercise 2 : Get pre-shower and module 1 tower 1 C, S channel cid using utility.getcid()
    TBcid pscid = utility.getcid(TBdetector::detid::preshower);
    TBcid M1T1C_cid = utility.getcid(); // Your answer here
    TBcid M1T1S_cid = utility.getcid(); // Your answer here

    /*
    When drawing average time structure, we need to remove particles which are likely to miss our module
    This can be done by setting cut using DWC X and Y correlation plot
    For example, we expect particle entering DWC1 center to enter near the DWC2 center
    If not, then particle might be tilted from (DWC1 - DWC2 - module) axis, and will miss the module
    We'll have to exclude such event to reduce events that does not deposit energy in the module
    
    We'll use 1.5mm threshold for this exercise
    That is, abs[ (DWC1 X distance from DWC1 center) - (DWC2 X distance from DWC2 center) ] <= 1.5 mm and same applied for Y position also
    To do so, we need to know the center position of DWC1 and DWC2
    Since they're not "physically" aligned, their center postion is not (0 mm, 0 mm) but shifted
    We'll use the DWC root file saved from previous exercise to get the center position of DWC1 and DWC2
    */
    
    // Get DWC root file created in previous exercise
    TFile* dwcFile = TFile::Open(("./dwc/dwc_Run_" + std::to_string(runNum) + ".root").c_str());
    // Get DWC 1 & 2 2D histograms from DWC root file
    TH2D* dwc1_pos   = (TH2D*) dwcFile->Get("dwc1_pos");
    TH2D* dwc2_pos   = (TH2D*) dwcFile->Get("dwc2_pos");
    // getDWCoffset : Function defined in function.cc for calculating DWC center position.
    // Receive dwc 2D histogram and returns offset vector >> vector[0] = DWC X offset, vector[1] = DWC Y offset
    std::vector<float> dwc1_offset = getDWCoffset(dwc1_pos); // dwc1_offset.at(0) == dwc1 X offset from 0 mm, dwc1_offset.at(1) == dwc1 Y offset from 0 mm
    std::vector<float> dwc2_offset = getDWCoffset(dwc2_pos);

    // Evt Loop start
    // To draw average time structure, we first stack all waveforms into one waveform histogram
    // At the end, we'll divide the waveforms with (# of entries) using histogram->Scale() function
    for (int iEvt = 0; iEvt < totalEntry; iEvt++) {
        if( iEvt % 1000 == 0 ) printProgress(iEvt + 1, totalEntry);

        evtChain->GetEntry(iEvt);
        // Getting TBwaveform data using TBcid, and waveform from TBwaveform
        TBwaveform psData = anEvt->data(pscid);
        std::vector<short> psWaveform = psData.waveform();

        // Exercise 3 : Get waveform of both Module 1 tower 1 S and C channels
        TBwaveform M1T1C_data = ; // Your answer here
        std::vector<short> M1T1C_waveform = ; // Your answer here
                                ; // Your answer here
                                            ; // Your answer here

        // This part is for getting DWC data -> waveform -> peak timing and applying DWC correlation cut
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
        dwc1_correctedPos->Fill(dwc1_correctedPosition.at(0), dwc1_correctedPosition.at(1));
        dwc2_correctedPos->Fill(dwc2_correctedPosition.at(0), dwc2_correctedPosition.at(1));
        // For DWC correlation cut (1.5mm threshold), use function defined in functions.cc : dwcCorrelationPID()
        // This function calculates the X and Y position difference between the two DWCs, and return true if both differences are smaller than 1.5 (mm)
        // If not passed (== X or Y position difference of DWC 1 and 2 is larger than 1.5 mm), skip event loop
        if (! dwcCorrelationPID(dwc1_correctedPosition, dwc2_correctedPosition, 1.5f) ) continue;
        dwc1_correctedPos_PID->Fill(dwc1_correctedPosition.at(0), dwc1_correctedPosition.at(1));
        dwc2_correctedPos_PID->Fill(dwc2_correctedPosition.at(0), dwc2_correctedPosition.at(1));

        for (int bin = 0; bin < 1000; bin++) {
            int waveformBin = bin + 1;
            // Pre-shower
            psHist        ->Fill(bin, psWaveform[waveformBin]);
            // Exercise 4 : Fill M1T1 C, S histogram with corresponding waveforms
            M1T1CHist     ->Fill(, ); // Your answer here
            M1T1SHist     ->Fill(, ); // Your answer here
        }
    }

    // Here we get the total # of events survived after DWC correlation cut
    // We use inverse of (# of events) to scale our histograms to get average time structure ( == average waveform)
    float survivedEntries = (float) psHist->GetEntries() / 1000.;
    float scale = 1. / survivedEntries;

    std::string outFile = "./avgTimeStructure/avg_Run_" + std::to_string(runNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    psHist->Scale( scale );
    psHist->GetYaxis()->SetRangeUser(-100, 4096);
    psHist->Write();

    // M1 C
    M1T1CHist->Scale( scale );
    M1T1CHist->GetYaxis()->SetRangeUser(-100, 4096);
    M1T1CHist->Write();

    // M1 S
    M1T1SHist->Scale( scale );
    M1T1SHist->GetYaxis()->SetRangeUser(-100, 4096);
    M1T1SHist->Write();

    dwc1_correctedPos->Write();
    dwc2_correctedPos->Write();
    dwc1_correctedPos_PID->Write();
    dwc2_correctedPos_PID->Write();

    outputRoot->Close();
}
