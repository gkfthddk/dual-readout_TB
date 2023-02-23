#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "TROOT.h"
#include "TStyle.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

//////////////////////////////////////////////////////////////////////////////
// For loading ntuples
//////////////////////////////////////////////////////////////////////////////

// Change "filePath" to your ntuple file path
TChain* getNtupleChain(int runNumber) {
    TChain* evtChain = new TChain("events");
    for (int fn = 0; fn < 50; fn++) {
        std::string fileName = "ntuple_Run_" + std::to_string(runNumber) + "_Wave_" + std::to_string(fn) + ".root"; // This is your ntuple name, no need to change it
        std::string filePath = "/gatbawi/dream/ntuple/waveform/Run_"  + std::to_string(runNumber) + "/" + fileName; // Change to your own ntuple file path here
        if ( !access(filePath.c_str(), F_OK) ){
            std::cout << "Ntuple # : " << fn << " added to TChain : " << filePath << std::endl; 
            evtChain->Add(filePath.c_str());
        }
    }
    return evtChain;
}

//////////////////////////////////////////////////////////////////////////////
// For DWC position
//////////////////////////////////////////////////////////////////////////////

// Transfer waveform bin index to time (nanosecond)
float getTime(float bin) {
  return 400. * (bin / 1000.);
}

// Get peak bin index from raw waveform (peak == minimum value)
int getPeakBin(std::vector<short> waveform) {
    return ( std::min_element(waveform.begin()+1, waveform.end()-23) - waveform.begin() );
}

// Get peak timing from raw waveform (peak == minimum value)
float getPeakTime(std::vector<short> waveform) {
    int peakBin = getPeakBin(waveform);
    return (getTime(peakBin));
}

// Calculate DWC position using DWC waveform peak timing (timing == nanosecond)
// DWC X-axis position : { (Waveform peak timing from DWC Right ch) - (Waveform peak timing from DWC Left ch) } * (calibrated DWC x-axis slope value) + (calibrated DWC x-axis intercept value)
// DWC Y-axis position : { (Waveform peak timing from DWC Up ch) - (Waveform peak timing from DWC Down ch) } * (calibrated DWC y-axis slope value) + (calibrated DWC y-axis intercept value)
// Subtract x and y offset from calculated DWC x-axis and y-axis position to make mean value of DWC position plot aligns to center. That is, (0 mm, 0 mm) position.
// x and y offset can be obtained from non-corrected DWC position plot's mean value.
std::vector<float> getDWC1position(std::vector<float> dwc1Time, float xOffset = 0.f, float yOffset = 0.f) {
    static float dwc1horizontalSlope = 0.1740806676;
    static float dwc1horizontalOffset = 0.1680572999;
    static float dwc1VerticalSlope = -0.17424779576;
    static float dwc1VerticalOffset = -0.053701300;

    std::vector<float> dwc1Position;

    float horizontalPos =  ( (float)(dwc1Time.at(0) - dwc1Time.at(1)) * dwc1horizontalSlope ) + dwc1horizontalOffset - xOffset;
    float verticalPos = ( (float)(dwc1Time.at(2) - dwc1Time.at(3)) * dwc1VerticalSlope ) + dwc1VerticalOffset - yOffset;

    dwc1Position.push_back(horizontalPos);
    dwc1Position.push_back(verticalPos);

    return dwc1Position;
}

std::vector<float> getDWC2position(std::vector<float> dwc2Time, float xOffset = 0.f, float yOffset = 0.f) {
    static float dwc2horizontalSlope = 0.17257273;
    static float dwc2horizontalOffset = 0.579927452;
    static float dwc2VerticalSlope = -0.1741203164;
    static float dwc2VerticalOffset = -0.278179655;

    std::vector<float> dwc2Position;

    float horizontalPos = ( (float)(dwc2Time.at(0) - dwc2Time.at(1)) * dwc2horizontalSlope ) + dwc2horizontalOffset - xOffset;
    float verticalPos = ( (float)(dwc2Time.at(2) - dwc2Time.at(3)) * dwc2VerticalSlope ) + dwc2VerticalOffset - yOffset;

    dwc2Position.push_back(horizontalPos);
    dwc2Position.push_back(verticalPos);

    return dwc2Position;
}

// Calculate DWC position using DWC waveform peak timing (timing == clock value)
// To translate clock value to nanosecond : (clock * 25. / 1000.) == (nanosecond)
// Clock values can be obtained from fast mode event's evt.timing() function
// Other logics are same as above
std::vector<float> getDWC1positionUsingClock(std::vector<float> dwc1Time) {
    static float dwc1horizontalSlope = -0.1740806676;
    static float dwc1horizontalOffset = -0.1680572999;
    static float dwc1VerticalSlope = -0.17424779576;
    static float dwc1VerticalOffset = -0.053701300;

    std::vector<float> dwc1Position;

    float horizontalPos = -( ( (float)(dwc1Time.at(0) - dwc1Time.at(1)) * 25. / 1000. * dwc1horizontalSlope ) + dwc1horizontalOffset );
    float verticalPos = ( (float)(dwc1Time.at(2) - dwc1Time.at(3)) * 25. / 1000. * dwc1VerticalSlope ) + dwc1VerticalOffset;

    dwc1Position.push_back(horizontalPos);
    dwc1Position.push_back(verticalPos);

    return dwc1Position;
}

std::vector<float> getDWC2positionUsingClock(std::vector<float> dwc2Time) {
    static float dwc2horizontalSlope = -0.17257273;
    static float dwc2horizontalOffset = -0.579927452;
    static float dwc2VerticalSlope = -0.1741203164;
    static float dwc2VerticalOffset = -0.278179655;

    std::vector<float> dwc2Position;

    float horizontalPos = -( ( (float)(dwc2Time.at(0) - dwc2Time.at(1)) * 25. / 1000. * dwc2horizontalSlope ) + dwc2horizontalOffset );
    float verticalPos = ( (float)(dwc2Time.at(2) - dwc2Time.at(3)) * 25. / 1000. * dwc2VerticalSlope ) + dwc2VerticalOffset;

    dwc2Position.push_back(horizontalPos);
    dwc2Position.push_back(verticalPos);

    return dwc2Position;
}

// From non-corrected DWC 2D histogram, get its x and y position mean value and return them
std::vector<float> getDWCoffset(TH2D* dwcHist) {
    float xOffset = dwcHist->GetMean(1);
    float yOffset = dwcHist->GetMean(2);

    std::vector<float> dwcOffset;
    dwcOffset.push_back(xOffset);
    dwcOffset.push_back(yOffset);

    return dwcOffset;
}

//////////////////////////////////////////////////////////////////////////////
// For PID
//////////////////////////////////////////////////////////////////////////////

// Give DWC1, 2 correlation constraint, threshold set to 1.5mm for default.
// This means particle entering DWC1 at (x mm, y mm) should enter DWC2 at (x +- 1.5 mm, y +-1.5 mm)
bool dwcCorrelationPID(std::vector<float> dwc1_correctedPosition, std::vector<float> dwc2_correctedPosition, float threshold = 1.5f) {
    bool passed = false;

    float x_diff = std::abs( dwc1_correctedPosition.at(0) - dwc2_correctedPosition.at(0) );
    float y_diff = std::abs( dwc1_correctedPosition.at(1) - dwc2_correctedPosition.at(1) );

    if ( (x_diff <= threshold) && (y_diff <= threshold) ) passed = true;

    return passed;
}

//////////////////////////////////////////////////////////////////////////////
// print progress
//////////////////////////////////////////////////////////////////////////////

// Simple utility to print progress
void printProgress(const int currentStep, const int totalStep) {
    float progress = (float) currentStep / totalStep;
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; i++) {
        if ( i < pos ) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "]  " << currentStep << "/" << totalStep << "  " << int(progress * 100.0) << "%\r";
    if ( (currentStep + 1) != totalStep) std::cout.flush();
    else {
        std::cout<< std::endl;
        std::cout<< "Finished!" << std::endl;
    }
}
