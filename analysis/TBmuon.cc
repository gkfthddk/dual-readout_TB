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

// This macro will draw DWC 1 & 2's position and correlation plots from waveform ntuple
// This macro will store 2D histograms in .root format

// How to execute
// On local or batch, run the following command :
// ./TBdwc.exe <run number> <max # of events to process>

int main(int argc, char** argv) {
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

    // Preparing 2D histograms for DWC 1 & 2 position, correlation plots
    TH1D* muon_pos = new TH1D("muon_pos", "muon", 100, -5000., 50000.);
    /*
    - Structure of DWC
    DWC (Delayed wire chamber) is used to find out the particle's position.
    One DWC uses 4 channels (Right, Left, Up, Down channels) to determine incident particle's position.
    Details about DWC and how to get position of it can be found here : https://twiki.cern.ch/twiki/pub/DREAM/DwcCalibration/DWC_calibration.pdf
    */

    /*
    - Get channel IDs
    To get started, we need to know cids (channel IDs) of each DWC 1 and 2 to draw their information from ntuple
    One can use TBcid(int DAQ_number, int channel_number) to get DWC's cid
    From README.md, open the mapping info link : https://docs.google.com/spreadsheets/d/1ZFLbl-4OOqhWQlKlCitL4Pl__JsoVFMME0hQWy3N-M8/edit#gid=257337700
    and find out which MID(== DAQ number) and channel number the DWC belongs to
    For example, from mapping info one can find that DWC1-R matches to mid 1 channel 17 which means DWC1 right channel id can be declared as : TBcid(1, 17)
    */

    // For example, cid of DWC 1 Right, Left, Up, Down channel can be determined as follows
    //TBcid muon_cid = utility.getcid(TBdetector::detid::preshower);
    TBcid muon_cid = utility.getcid(TBdetector::detid::muon);

    utility.loadped( ("/gatbawi/dream/ped/mean/Run" + std::to_string(runNum) + "_pedestalHist_mean.root").c_str() );
    float muonPed = utility.retrievePed(muon_cid);
    // Exercise 1 : Referring to mapping info, get channel IDs of DWC 2
    float maxADC=0.;
    float minADC=0.;
    // Start Evt Loop to draw DWC plots
    for (int iEvt = 0; iEvt < totalEntry; iEvt++) {
        evtChain->GetEntry(iEvt);

        // Get each DWC1(2) R, L, U, D channels data using evt->data(TBcid cid);
        TBwaveform muon_data = anEvt->data(muon_cid);
        // Exercise 2 : Get DWC2 R, L, U, D channels data
        // Get each DWC1 R, L, U, D channels waveform using data.waveform();
        std::vector<short> muon_waveform = muon_data.waveform();

        /*
        To get DWC position from its calibration value, one needs to know peak timing of the waveform
        Here, peak timing is the time of waveform reaching at its peak position (minimum value)
        We're using 400 ns time window for waveform with 1000 bins
        Therefore, each bin in waveform vector corresponds to 0.4 ns
        So, if the waveform had its peak (minimum value) at N bin, we can get peak timing = (float) (N * 0.4)
        The function for converting bin <-> time is defined in functions.cc : getTime()

        Now all we need to do is to get peak position (== bin number) of each waveform
        This can be done using std::min_element(), which returns the index of vector with minimum value
        We're looking into (waveform.begin() + 1) ~ (waveform.end() - 23), so
        std::min_element( waveform.begin() + 1, waveform.end() - 23) - waveform.begin()
        will return the bin number with the minimum value ( == peak position )
        The function for getting peak position is defined in functions.cc : getPeakBin()

        Combining getTime() and getPeakBin() together, we cat get peak timing of each DWC1(2) R, L, U, D waveforms
        The function for doing this is defined in functions.cc : getPeakTime(std::vector<short> waveform)
        */

        // Get DWC 1 R, L, U, D peak timing from their waveforms
        //float muon_peak = getPeakTime(muon_waveform);
        float muonIntADC = 0.f;
        for (int bin = 0; bin < 1000; bin++) {
            int waveformBin = bin + 1;
            muonIntADC += muonPed - muon_waveform[waveformBin];
        }


        /*
        To get DWC position from timing value, one need to do as follows (refer to DWC info pdf page 3)
        
        X-position = (timing from right - timing from left) * horizontal_Slope + horizontal_Offset
        Y-position = (timing from up - timing from down) * vertical_Slope + vertical_Offset

        Here, each slope and offsets are calibrated constants decided from test beam
        DWC 1 :
        float dwc1_horizontal_Slope = 0.1740806676;
        float dwc1_horizontal_Offset = 0.1680572999;
        float dwc1_Vertical_Slope = -0.17424779576;
        float dwc1_Vertical_Offset = -0.053701300;
        DWC 2 :
        float dwc2_horizontal_Slope = 0.17257273;
        float dwc2_horizontal_Offset = 0.579927452;
        float dwc2_Vertical_Slope = -0.1741203164;
        float dwc2_Vertical_Offset = -0.278179655;
        */


        // Filling 2D histogram to plot DWC 1 and 2 positions
        muon_pos->Fill(muonIntADC);
        if(maxADC<muonIntADC){
           maxADC=muonIntADC;
        }
        if(minADC>muonIntADC){
           minADC=muonIntADC;
        }

        printProgress(iEvt + 1, totalEntry);
    }
    printf("max %g min %g \n ",maxADC,minADC);
    // Saving the DWC position & correlation plots in root file in ./dwc directory.
    std::string outFile = "./muon/muon_Run_" + std::to_string(runNum) + ".root";
    TFile* outputRoot = new TFile(outFile.c_str(), "RECREATE");
    outputRoot->cd();

    muon_pos->Write();

    outputRoot->Close();
}
