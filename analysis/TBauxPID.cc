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

// This macro will draw integrated ADC of pre-shower detector, and fit 1-mip peak using gaussian function
// This macro will store 1D histograms and fit functions in .root format

// How to execute
// On local or batch, run the following command :
// ./TBauxPID.exe <run number> <max # of events to process>

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
    
    // Get DWC center info (offsets) for DWC correlation cut
    TFile* dwcFile = TFile::Open(("./dwc/dwc_Run_" + std::to_string(runNum) + ".root").c_str());
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
    TH2D* dwc1_pos   = (TH2D*) dwcFile->Get("dwc1_pos");
    TH2D* dwc2_pos   = (TH2D*) dwcFile->Get("dwc2_pos");
    std::vector<float> dwc1_offset = getDWCoffset(dwc1_pos); // dwc1_offset.at(0) == dwc1 X offset from (0, 0), dwc1_offset.at(1) == dwc1 Y offset from (0, 0)
    std::vector<float> dwc2_offset = getDWCoffset(dwc2_pos);

    // Preparing histograms
    TH1F* psIntHist = new TH1F("psIntADC", "psIntADC;IntADC;evt", 1000, -5000., 200000.);
    TH1F* psIntHist_PID = new TH1F("psIntADC_PID", "psIntADC_PID;IntADC;evt", 1000, -5000., 200000.);
    TH2D* dwc1_correctedPos = new TH2D("dwc1_correctedPos", "dwc1_correctedPos;mm;mm;events", 480, -120., 120., 480, -120., 120.);
    TH2D* dwc2_correctedPos = new TH2D("dwc2_correctedPos", "dwc2_correctedPos;mm;mm;events", 480, -120., 120., 480, -120., 120.);

    /*
    To get integrated ADC, we must define integration range
    From previous exercise, we decided to use bin 220 ~ 390 for pre-shower detector
    Next we subtract waveform from pedestal, and sum all the ADC values in integral range
    
    Pedestals are already measured & stored in pedestal root file
    One can load pedestal information as follows :
    utility.loadped( ("/gatbawi/dream/ped/mean/Run" + std::to_string(runNum) + "_pedestalHist_mean.root").c_str() );

    To get pedestal value of specific channel from pedestal root file, one can use : utility.retrievePed(TBcid cid)
    For example, to get pedestal value of pre-shower detector, one can do :
    float ps_ped = utility.retrievePed(ps_Cid);
    where ps_Cid is the pre-shower detector channel ID
    */

    // Exercise 1 : Get pre-shower cid
    // Get channel IDs
    TBcid pscid = ; // Your answer here

    // Get pedestals : pedestals are needed to calculate integrated ADC
    utility.loadped( ("/gatbawi/dream/ped/mean/Run" + std::to_string(runNum) + "_pedestalHist_mean.root").c_str() );
    float psPed = utility.retrievePed(pscid);

    // Evt Loop
    for (int iEvt = 0; iEvt < totalEntry; iEvt++) {
        printProgress(iEvt + 1, totalEntry);

        evtChain->GetEntry(iEvt);
        
        // Getting DWC peak timing for DWC correlation cut
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

        // Exercise 2 : Using pre-shower cid, get its waveform
        TBwaveform psData = ; // Your answer here
        std::vector<short> psWaveform = ; // Your answer here

        // To get integrated ADC, we subtract each pre-shower waveform bin value from pre-shower pedestals, and sum them up
        float psIntADC = 0.f;
        for (int bin = 220; bin < 390; bin++) { // Recall that integral range is [220, 390)
            int waveformBin = bin + 1;
            psIntADC += psPed - psWaveform[waveformBin];
        }

        // PS Int. ADC plot before PID
        psIntHist->Fill(psIntADC);
    
        // For DWC correlation cut
        if (! dwcCorrelationPID(dwc1_correctedPosition, dwc2_correctedPosition, 1.5f) ) continue;

        // PS Int. ADC plot after PID
        psIntHist_PID->Fill(psIntADC);
    }

    std::string outFile = "./auxPID/auxPID_Run_" + std::to_string(runNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    // ----------------------------------------------------------------------------------------------------------------------
    // Exercise 3 : First, use the fit function with default range 9000 ~ 15000 to fit the 1-mip peak of pre-shower detector
    //              Then, referring to README.md, find more proper fit function range to get better Chi2/NDF value
    // ----------------------------------------------------------------------------------------------------------------------
    // Using TF1, define gaussian fit function for fitting 1-mip peak of pre-shower detector
    // The range 9000~15000 need to be optimized to give best Chi2/NDF after fitting PS int. ADC
    TF1* psFitFunc = new TF1("psFitFunc", "gaus", 8000, 18000); // Your answer here, change the range after checking the fit result
    
    psIntHist->Write();
    psIntHist_PID->Fit("psFitFunc", "R"); // histogram->Fit("fit function name", "R"); will fit the histogram with corresponding fit function 
    psIntHist_PID->Write();
    psFitFunc->Write();

    dwc1_pos->Write();
    dwc2_pos->Write();
    dwc1_correctedPos->Write();
    dwc2_correctedPos->Write();
    
    outputRoot->Close();
}
