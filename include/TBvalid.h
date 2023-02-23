#ifndef TBvalid_h
#define TBvalid_h 1

#include "TBmid.h"
#include "TBevt.h"

#include "TROOT.h"
#include "TStyle.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdint>
#include <vector>

#include <boost/python.hpp>

class TBvalid {
    public:
        TBvalid() {};
        ~TBvalid() {};

        void checkRatio(TH1F* den, TH1F* num, const std::string histName, const std::string outDir);
        void drawRatio(TH1F* den, TH1F* num, const std::string histName, const std::string outDir);

        TH1F* drawFastHistFromData(TBcid cid, const std::string histName);
        TH1F* drawFastHistFromData(const std::vector<std::vector<std::string>>& datList, TBcid cid, const std::string histName);
        TH1F* drawFastHistFromNtuple(TBcid cid, const std::string histName);
        TH1F* drawFastHistFromNtuple(const std::vector<std::string>& ntupleList, TBcid cid, const std::string histName);

        TH1F* drawWaveHistFromData(TBcid cid, const std::string histName);
        TH1F* drawWaveHistFromData(const std::vector<std::vector<std::string>>& datList, TBcid cid, const std::string histName);
        TH1F* drawWaveHistFromNtuple(TBcid cid, const std::string histName);
        TH1F* drawWaveHistFromNtuple(const std::vector<std::string>& ntupleList, TBcid cid, const std::string histName);
        
        bool simpleValidFast(const std::vector<std::vector<std::string>>& datList,const std::vector<std::string>& ntupleList, TBcid cid);
        bool simpleValidWave(const std::vector<std::vector<std::string>>& datList,const std::vector<std::string>& ntupleList, TBcid cid);

        bool py_simpleValidFast(TBcid cid);
        bool py_simpleValidWave(TBcid cid);

        void checkFastTrigNum();
        void checkWaveTrigNum();
        
        std::vector<std::vector<std::string>> readDataList(const boost::python::list& py_datList);
        std::vector<std::string> readNtupleList(const boost::python::list& py_ntupleList);

        void printProgress(const int currentStep, const int totalStep);

        void setDataList(const boost::python::list& py_datList) { datList_ = readDataList(py_datList); }
        void setNtupleList(const boost::python::list& py_ntupleList) { ntupleList_ = readNtupleList(py_ntupleList); }

    private:
        std::vector<std::vector<std::string>> datList_;
        std::vector<std::string> ntupleList_;          
};

#endif