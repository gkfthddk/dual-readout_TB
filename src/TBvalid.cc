#include "TBmid.h"
#include "TBevt.h"
#include "TBread.h"
#include "TBvalid.h"
#include "TButility.h"
#include "TBdetector.h"

#include <stdexcept>
#include <boost/python.hpp>

#include <TFile.h>
#include <TTree.h>

void TBvalid::drawRatio(TH1F* num, TH1F* den, const std::string histName, const std::string outDir) {
    TH1F* h_fromNtuple = num;
    TH1F* h_fromData = den;
    TH1F* h_copy = (TH1F*)num->Clone("h_copy");

    TCanvas* c1 = new TCanvas("c", "c", 1000, 800);
    c1->cd();

    TPad* p1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);

    p1->SetGridx();
    p1->SetBottomMargin(0);
    p1->Draw();
    p1->cd();

    h_fromData->SetTitle(histName.c_str());
    h_fromData->SetLineColor(kBlue);
    h_fromData->GetXaxis()->SetLabelSize(0);
    h_fromData->GetYaxis()->SetTitle("Evt");
    h_fromData->GetYaxis()->SetTitleOffset(1.5);
    h_fromData->GetYaxis()->SetTitleSize(0.035);
    h_fromData->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    h_fromData->Draw("Hist");
    p1->Update();
    TPaveStats* stat_Data_C = (TPaveStats*)p1->GetPrimitive("stats");
    stat_Data_C->SetName("Data_C");
    stat_Data_C->SetTextColor(kBlue);
    stat_Data_C->SetY1NDC(0.6); stat_Data_C->SetY2NDC(0.75);
    stat_Data_C->SetX1NDC(0.7); stat_Data_C->SetX2NDC(0.9);

    h_fromNtuple->SetTitle("");
    h_fromNtuple->SetMarkerStyle(21);
    h_fromNtuple->SetMarkerSize(0.6);
    h_fromNtuple->SetMarkerColor(kRed);
    h_fromNtuple->Draw("Hist&sames&p");
    p1->Update();
    TPaveStats* stat_Ntuple_C = (TPaveStats*)p1->GetPrimitive("stats");
    stat_Ntuple_C->SetName("Ntuple_C");
    stat_Ntuple_C->SetTextColor(kRed);
    stat_Ntuple_C->SetY1NDC(0.75); stat_Ntuple_C->SetY2NDC(.9);
    stat_Ntuple_C->SetX1NDC(0.7); stat_Ntuple_C->SetX2NDC(0.9);

    c1->cd();
    TPad* p2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    p2->SetGridx();
    p2->SetGridy();
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0);
    p2->SetBorderSize(0);
    p2->SetBorderMode(0);
    p2->Draw();
    p2->cd();

    h_copy->SetMarkerColor(kBlack);
    h_copy->SetMarkerStyle(21);
    h_copy->SetMarkerSize(0.5);
    h_copy->Sumw2();
    h_copy->SetTitle("");
    h_copy->SetLineWidth(0);
    h_copy->Divide(h_fromData);

    h_copy->GetYaxis()->SetNdivisions(505);
    h_copy->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
    h_copy->GetYaxis()->SetTitle("Ntuple/Data");
    h_copy->GetYaxis()->SetTitleSize(0.08);
    h_copy->GetYaxis()->SetTitleOffset(0.5);
    h_copy->GetYaxis()->SetLabelSize(0.085);

    h_copy->GetXaxis()->SetTitle("Int ADC");
    h_copy->GetXaxis()->SetTitleSize(0.09);
    h_copy->GetXaxis()->SetTitleOffset(0.9);
    h_copy->GetXaxis()->SetLabelSize(0.08);

    h_copy->SetStats(0);
    h_copy->SetMinimum(0.99);
    h_copy->SetMaximum(1.01);

    h_copy->Draw("p");    

    p2->Update();

    gPad->SetBottomMargin(0.2);

    // Part where checking bin contents and drawing ratio if bin contains non-1 values
    int nBins = h_copy->GetNbinsX();
    // bool drawRatio = false;
    bool drawRatio = true;
    for (int bin = 0; bin < nBins; ++bin) {
        float copy_binValue = (float) h_copy->GetBinContent(bin);
        float num_binValue = (float) h_fromNtuple->GetBinContent(bin);
        float den_binValue = (float) h_fromData->GetBinContent(bin);
        if ( copy_binValue != 1.0f && (num_binValue != 0.0f && den_binValue != 0.0f) ) {
            std::cout << "The bin value for bin " << bin << " is not 1, but " << copy_binValue << ", drawing Ratio plot..." << std::endl;
            drawRatio = true;
        }
    }

    if (drawRatio) {
        c1->Draw();
        c1->SaveAs(( outDir + histName + ".png").c_str());
    }
    else {
        std::cout << "All bins are fine, good to proceed..." << std::endl;
    }

    delete p1;
    delete p2;
    delete c1;
}

void TBvalid::checkRatio(TH1F* num, TH1F* den, const std::string histName, const std::string outDir) {
    TH1F* h_fromNtuple = num;
    TH1F* h_fromData = den;
    TH1F* h_copy = (TH1F*)num->Clone("h_copy");

    h_copy->Divide(h_fromData);

    // Part where checking bin contents and drawing ratio if bin contains non-1 values
    int nBins = h_copy->GetNbinsX();
    bool drawRatio = false;
    for (int bin = 0; bin < nBins; ++bin) {
        float copy_binValue = (float) h_copy->GetBinContent(bin);
        float num_binValue = (float) h_fromNtuple->GetBinContent(bin);
        float den_binValue = (float) h_fromData->GetBinContent(bin);
        if ( copy_binValue != 1.0f && (num_binValue != 0.0f && den_binValue != 0.0f) ) {
            std::cout << "The bin value for bin " << bin << " is not 1, but " << copy_binValue << std::endl;
            drawRatio = true;
        }
    }

    if (drawRatio) {
        std::cout << "[ERROR] Contents mismatch between data and ntuple" << std::endl;
    }
    else {
        std::cout << "All bins are fine, good to proceed..." << std::endl;
    }
}

TH1F* TBvalid::drawWaveHistFromData(TBcid cid, const std::string histName) {
    if (this->datList_.size() == 0) {
        std::cout << "No data list found, please use setDataList() first" << std::endl;
        return NULL;
    }
    TH1F* ptr = drawWaveHistFromData(this->datList_, cid, histName);
    return ptr;
}

TH1F* TBvalid::drawWaveHistFromData(const std::vector<std::vector<std::string>>& datList, TBcid cid, const std::string histName) {
    std::cout << "Drawing waveform histogram from data files..." << std::endl;
    TH1F* hist = new TH1F(histName.c_str(), histName.c_str(), 1000, 0., 1000.);
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue);

    const int numOfMID = 15;
    const int numOfFiles = datList.size();

    std::vector<int> entryPerMID(numOfMID);
    std::vector< std::vector<int> > entryPerFile( numOfFiles, std::vector<int>(numOfMID) );
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        for (unsigned MID = 0; MID < numOfMID; MID++) {
        FILE* fp = fopen(datList[fileNum][MID].c_str(), "rb");
        fseek(fp, 0L, SEEK_END);
        unsigned long long file_size = ftell(fp);
        fclose(fp);
        entryPerMID[MID] += (int)(file_size / 65536);
        entryPerFile[fileNum][MID] = (int)(file_size / 65536);
        }
    }
    int totalEntry = *std::min_element(entryPerMID.begin(), entryPerMID.end());
    std::cout << "Total Entry : " << totalEntry << std::endl;

    std::vector<int> currentOpenFileNum(numOfMID);
    std::vector<int> entryCounted(numOfMID);
    std::vector<FILE*> files(numOfMID);
    for (unsigned MID = 0; MID < numOfMID; MID++) {
        files[MID] = fopen(datList[0][MID].c_str(), "rb");
    }
    TBevt<TBwaveform>* anEvtFromData = new TBevt<TBwaveform>();
    // Start evt loop
    for (unsigned iEvt = 0; iEvt < totalEntry; iEvt++){
        // if count[MID] exceeds # of entry of corresponding file, move to next file
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            // currently opened file number for each MID
            if ( entryPerFile[currentOpenFileNum[MID]][MID] == entryCounted[MID] ) {
                fclose(files[MID]);
                currentOpenFileNum[MID]++;
                files[MID] = fopen(datList[currentOpenFileNum[MID]][MID].c_str(), "rb");
                entryCounted[MID] = 0;
            }
        }

        TBread reader = TBread();
        std::vector<TBmid<TBwaveform>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            MIDs[MID] = reader.readWaveform(files[MID]);
        }
        anEvtFromData->setTCB(MIDs[0].evt());
        anEvtFromData->set(MIDs);

        auto dataFromData = anEvtFromData->data(cid);
        auto waveform = dataFromData.waveform();

        for (int bin = 0; bin < 1000; bin++){
            hist->Fill(bin, waveform[bin+1]);
        }

        for (unsigned MID = 0; MID < numOfMID; MID++) entryCounted[MID]++;
        printProgress( (iEvt + 1) , totalEntry);
    }

    TCanvas* c = new TCanvas("","");
    c->cd();
    hist->Draw("Hist");
    // c->SaveAs(("./" + histName + ".png").c_str());

    return hist;
}

TH1F* TBvalid::drawFastHistFromData(TBcid cid, const std::string histName) {
    if (this->datList_.size() == 0) {
        std::cout << "No data list found, please use setDataList() first" << std::endl;
        return NULL;
    }
    TH1F* ptr = drawFastHistFromData(this->datList_, cid, histName);
    return ptr;
}

TH1F* TBvalid::drawFastHistFromData(const std::vector<std::vector<std::string>>& datList, TBcid cid, const std::string histName) {
    std::cout << "Drawing fastmode histogram from data files..." << std::endl;
    TH1F* hist = new TH1F(histName.c_str(), histName.c_str(), 1000, -1000., 50000);
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlue);

    const int numOfMID = 15;
    const int numOfFiles = datList.size();

    std::vector<int> entryPerMID(numOfMID);
    std::vector< std::vector<int> > entryPerFile( numOfFiles, std::vector<int>(numOfMID) );
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        for (unsigned MID = 0; MID < numOfMID; MID++) {
        FILE* fp = fopen(datList[fileNum][MID].c_str(), "rb");
        fseek(fp, 0L, SEEK_END);
        unsigned long long file_size = ftell(fp);
        fclose(fp);
        entryPerMID[MID] += (int)(file_size / 256);
        entryPerFile[fileNum][MID] = (int)(file_size / 256);
        }
    }
    int totalEntry = *std::min_element(entryPerMID.begin(), entryPerMID.end());
    std::cout << "Total Entry : " << totalEntry << std::endl;

    std::vector<int> currentOpenFileNum(numOfMID);
    std::vector<int> entryCounted(numOfMID);
    std::vector<FILE*> files(numOfMID);
    for (unsigned MID = 0; MID < numOfMID; MID++) {
        files[MID] = fopen(datList[0][MID].c_str(), "rb");
    }
    TBevt<TBfastmode>* anEvtFromData = new TBevt<TBfastmode>();
    // Start evt loop
    for (unsigned iEvt = 0; iEvt < totalEntry; iEvt++){
        // if count[MID] exceeds # of entry of corresponding file, move to next file
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            // currently opened file number for each MID
            if ( entryPerFile[currentOpenFileNum[MID]][MID] == entryCounted[MID] ) {
                fclose(files[MID]);
                currentOpenFileNum[MID]++;
                files[MID] = fopen(datList[currentOpenFileNum[MID]][MID].c_str(), "rb");
                entryCounted[MID] = 0;
            }
        }

        TBread reader = TBread();
        std::vector<TBmid<TBfastmode>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            MIDs[MID] = reader.readFastmode(files[MID]);
        }
        anEvtFromData->setTCB(MIDs[0].evt());
        anEvtFromData->set(MIDs);

        auto dataFromData = anEvtFromData->data(cid);

        hist->Fill(dataFromData.adc());

        for (unsigned MID = 0; MID < numOfMID; MID++) entryCounted[MID]++;
        printProgress( (iEvt + 1) , totalEntry);
    }

    TCanvas* c = new TCanvas("","");
    c->cd();
    hist->Draw("Hist");
    // c->SaveAs(("./" + histName + ".png").c_str());

    return hist;
}

TH1F* TBvalid::drawWaveHistFromNtuple(TBcid cid, const std::string histName) {
    if (this->ntupleList_.size() == 0) {
        std::cout << "No ntuple list found, please use setNtupleList() first" << std::endl;
        return NULL;
    }
    TH1F* ptr = drawWaveHistFromNtuple(this->ntupleList_, cid, histName);
    return ptr;
}

TH1F* TBvalid::drawWaveHistFromNtuple(const std::vector<std::string>& ntupleList, TBcid cid, const std::string histName) {
    std::cout << "Drawing waveform histogram from ntuples..." << std::endl;
    TH1F* hist = new TH1F(histName.c_str(), histName.c_str(), 1000, 0., 1000.);
    hist->SetLineWidth(2);
    hist->SetLineColor(kRed);

    const int numOfMID = 15;
    const int numOfFiles = ntupleList.size();

    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntry = fileChain->GetEntries();
    std::cout << "Total entries : " << totalEntry << std::endl;

    TBevt<TBwaveform>* anEvtFromNtuple = new TBevt<TBwaveform>();
    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);

    for (unsigned iEvt = 0; iEvt < totalEntry; ++iEvt) {
        fileChain->GetEntry(iEvt);
        auto dataFromNtuple = anEvtFromNtuple->data(cid);
        auto waveform = dataFromNtuple.waveform();

        for (int bin = 0; bin < 1000; ++bin) {
            hist->Fill(bin, waveform[bin+1]);
        }

        printProgress(iEvt +1, totalEntry);
    }
    TCanvas* c = new TCanvas("","");
    c->cd();
    hist->Draw("Hist");
    // c->SaveAs(("./" + histName + ".png").c_str());

    return hist;
}

TH1F* TBvalid::drawFastHistFromNtuple(TBcid cid, const std::string histName) {
    if (this->ntupleList_.size() == 0) {
        std::cout << "No ntuple list found, please use setNtupleList() first" << std::endl;
        return NULL;
    }
    TH1F* ptr = drawFastHistFromNtuple(this->ntupleList_, cid, histName);
    return ptr;
}

TH1F* TBvalid::drawFastHistFromNtuple(const std::vector<std::string>& ntupleList, TBcid cid, const std::string histName) {
    std::cout << "Drawing fastmode histogram from ntuples..." << std::endl;
    TH1F* hist = new TH1F(histName.c_str(), histName.c_str(), 1000, -1000., 50000);
    hist->SetLineWidth(2);
    hist->SetLineColor(kRed);

    const int numOfMID = 15;
    const int numOfFiles = ntupleList.size();

    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntry = fileChain->GetEntries();
    std::cout << "Total entries : " << totalEntry << std::endl;

    TBevt<TBfastmode>* anEvtFromNtuple = new TBevt<TBfastmode>();
    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);

    for (unsigned iEvt = 0; iEvt < totalEntry; ++iEvt) {
        fileChain->GetEntry(iEvt);
        auto dataFromNtuple = anEvtFromNtuple->data(cid);

        hist->Fill(dataFromNtuple.adc());

        printProgress(iEvt +1, totalEntry);
    }
    TCanvas* c = new TCanvas("","");
    c->cd();
    hist->Draw("Hist");
    // c->SaveAs(("./" + histName + ".png").c_str());

    return hist;
}

void TBvalid::checkFastTrigNum() {
    std::cout << "Checking trigger numbers..." << std::endl;
    const int numOfMID = 15;

    auto ntupleList = this->ntupleList_;

    const int numOfFiles = ntupleList.size();

    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntry = fileChain->GetEntries();
    std::cout << "Total entries : " << totalEntry << std::endl;

    std::vector<int> prev_tcb_trig_number(numOfMID);
    
    TBevt<TBfastmode>* anEvtFromNtuple = new TBevt<TBfastmode>();

    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);
    // Start event loop
    int err_count = 0;
    int war_count = 0;
    for (unsigned iEvt = 0; iEvt < totalEntry; ++iEvt) {
        fileChain->GetEntry(iEvt);

        std::vector<TBmid<TBfastmode>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; ++MID) {
            MIDs[MID] = anEvtFromNtuple->mid(MID);
        }

        // Validation on TCB trig # begin
        for (unsigned MID = 0; MID < numOfMID; ++MID) {

            TBmid<TBfastmode> refMID = MIDs[0];

            // Compare TCB trig # between all MIDs
            if (refMID.tcb_trig_number() != MIDs[MID].tcb_trig_number()) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : TCB trig number mismatch between MID 1 and MID : " << MID+1  << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : trig # of MID 1 : " << refMID.tcb_trig_number() << std::endl;
                std::cout << "[ERROR] : trig # of MID " << MID+1 << " : " << MIDs[MID].tcb_trig_number() << std::endl;
            }
            // Should be previous tcb # < current tcb # (if not, evt order mixed!)
            if ( (iEvt != 0) && (prev_tcb_trig_number[MID] >= MIDs[MID].tcb_trig_number()) ) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : Mixed order in TCB trig number found for MID : " << MID+1 << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : Previous tcb trig number : " << prev_tcb_trig_number[MID] << std::endl;
                std::cout << "[ERROR] : Current tcb trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
            }
            // Should be previous tcb # + 1 == current tcb # for fastmode (if not, missing some events)
            if ( (iEvt != 0) && ((prev_tcb_trig_number[MID] + 1) != MIDs[MID].tcb_trig_number()) ) {
                war_count++;
                std::cout << std::endl;
                std::cout << "[Warning] : Missing TCB trig number found in MID : " << MID+1 << " in event # : " << iEvt << std::endl;
                std::cout << "[Warning] : Previous tcb trig number : " << prev_tcb_trig_number[MID] << std::endl;
                std::cout << "[Warning] : Current tcb trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
            }
            // Sould be tcb # = local # (if not, one of TCB or DAQ malfunction might happened)
            if (MIDs[MID].tcb_trig_number() != MIDs[MID].local_trig_number()) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : Mismatch between local trig number and TCB trig number found in MID : " << MID+1 << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : TCB trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
                std::cout << "[ERROR] : DAQ trig number : " << MIDs[MID].local_trig_number() << std::endl;
            }

            prev_tcb_trig_number[MID] = MIDs[MID].tcb_trig_number();
        }

        printProgress(iEvt +1, totalEntry);
    }
            
    if (! (err_count || war_count) ) {
        std::cout << "No error or warning found, good to proceed" << std::endl;
    }
}

void TBvalid::checkWaveTrigNum() {
    std::cout << "Checking trigger numbers..." << std::endl;
    const int numOfMID = 15;

    auto ntupleList = this->ntupleList_;

    const int numOfFiles = ntupleList.size();

    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntry = fileChain->GetEntries();
    std::cout << "Total entries : " << totalEntry << std::endl;

    std::vector<int> prev_tcb_trig_number(numOfMID);
    
    TBevt<TBwaveform>* anEvtFromNtuple = new TBevt<TBwaveform>();

    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);
    // Start event loop
    int err_count = 0;
    int war_count = 0;
    for (unsigned iEvt = 0; iEvt < totalEntry; ++iEvt) {
        fileChain->GetEntry(iEvt);

        std::vector<TBmid<TBfastmode>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; ++MID) {
            MIDs[MID] = anEvtFromNtuple->mid(MID);
        }

        // Validation on TCB trig # begin
        for (unsigned MID = 0; MID < numOfMID; ++MID) {

            TBmid<TBfastmode> refMID = MIDs[0];

            // Compare TCB trig # between all MIDs
            if (refMID.tcb_trig_number() != MIDs[MID].tcb_trig_number()) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : TCB trig number mismatch between MID 1 and MID : " << MID+1  << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : trig # of MID 1 : " << refMID.tcb_trig_number() << std::endl;
                std::cout << "[ERROR] : trig # of MID " << MID+1 << " : " << MIDs[MID].tcb_trig_number() << std::endl;
            }
            // Should be previous tcb # < current tcb # (if not, evt order mixed!)
            if ( (iEvt != 0) && (prev_tcb_trig_number[MID] >= MIDs[MID].tcb_trig_number()) ) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : Mixed order in TCB trig number found for MID : " << MID+1 << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : Previous tcb trig number : " << prev_tcb_trig_number[MID] << std::endl;
                std::cout << "[ERROR] : Current tcb trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
            }
            // Should be previous tcb # + 1 == current tcb # for fastmode (if not, missing some events)
            // if ( (iEvt != 0) && ((prev_tcb_trig_number[MID] + 1) != MIDs[MID].tcb_trig_number()) ) {
            //     war_count++;
            //     std::cout << std::endl;
            //     std::cout << "[Warning] : Missing TCB trig number found in MID : " << MID+1 << " in event # : " << iEvt << std::endl;
            //     std::cout << "[Warning] : Previous tcb trig number : " << prev_tcb_trig_number[MID] << std::endl;
            //     std::cout << "[Warning] : Current tcb trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
            // }
            // Sould be tcb # = local # (if not, one of TCB or DAQ malfunction might happened)
            if (MIDs[MID].tcb_trig_number() != MIDs[MID].local_trig_number()) {
                err_count++;
                std::cout << std::endl;
                std::cout << "[ERROR] : Mismatch between local trig number and TCB trig number found in MID : " << MID+1 << " in event # : " << iEvt << std::endl;
                std::cout << "[ERROR] : TCB trig number : " << MIDs[MID].tcb_trig_number() << std::endl;
                std::cout << "[ERROR] : DAQ trig number : " << MIDs[MID].local_trig_number() << std::endl;
            }

            prev_tcb_trig_number[MID] = MIDs[MID].tcb_trig_number();
        }

        printProgress(iEvt +1, totalEntry);
    }
            
    if (! (err_count || war_count) ) {
        std::cout << "No error or warning found, good to proceed" << std::endl;
    }
}

bool TBvalid::py_simpleValidFast(TBcid cid) {
    if (this->datList_.size() == 0) {
        std::cout << "No data list found, please use setDataList() first" << std::endl;
        return false;
    }
    if (this->ntupleList_.size() == 0) {
        std::cout << "No ntuple list found, please use setNtupleList() first" << std::endl;
        return false;
    }
    bool valid = simpleValidFast(this->datList_, this->ntupleList_, cid);
    return valid;
}

bool TBvalid::simpleValidFast(const std::vector<std::vector<std::string>>& datList, const std::vector<std::string>& ntupleList, TBcid cid) {
    std::cout << "Quick fastmode validation..." << std::endl;
    const int numOfMID = 15;
    const int numOfFiles = datList.size();

    // Counting total entries from data file
    std::vector<int> entryPerMID(numOfMID);
    std::vector< std::vector<int> > entryPerFile( numOfFiles, std::vector<int>(numOfMID) );
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        for (unsigned MID = 0; MID < numOfMID; MID++) {
        FILE* fp = fopen(datList[fileNum][MID].c_str(), "rb");
        fseek(fp, 0L, SEEK_END);
        unsigned long long file_size = ftell(fp);
        fclose(fp);
        entryPerMID[MID] += (int)(file_size / 256);
        entryPerFile[fileNum][MID] = (int)(file_size / 256);
        }
    }
    int totalEntryfromData = *std::min_element(entryPerMID.begin(), entryPerMID.end());

    // Set TChain, count total entries from ntuple
    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntryfromNtuple = fileChain->GetEntries();

    // Error if total entry is different
    if (totalEntryfromData != totalEntryfromNtuple) {
        std::cout << "[ERROR] Ntuple total entry : " << totalEntryfromNtuple << " != Data total entry : " << totalEntryfromData << std::endl;
    }

    std::cout << "Total Entry : " << totalEntryfromData << std::endl;

    // Open data files and prepare reader
    TBevt<TBfastmode>* anEvtFromData = new TBevt<TBfastmode>();
    std::vector<int> currentOpenFileNum(numOfMID);
    std::vector<int> entryCounted(numOfMID);
    std::vector<FILE*> files(numOfMID);
    for (unsigned MID = 0; MID < numOfMID; MID++) {
        files[MID] = fopen(datList[0][MID].c_str(), "rb");
    }
    TBread reader = TBread();

    // Opening & reading ntuple with TChain
    TBevt<TBfastmode>* anEvtFromNtuple = new TBevt<TBfastmode>();
    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);
    
    // For error count
    int errorCount = 0;

    // Start evt loop
    for (unsigned iEvt = 0; iEvt < totalEntryfromData; iEvt++){
        // Reading data
        // if count[MID] exceeds # of entry of corresponding file, move to next file
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            // currently opened file number for each MID
            if ( entryPerFile[currentOpenFileNum[MID]][MID] == entryCounted[MID] ) {
                fclose(files[MID]);
                currentOpenFileNum[MID]++;
                files[MID] = fopen(datList[currentOpenFileNum[MID]][MID].c_str(), "rb");
                entryCounted[MID] = 0;
            }
        }
        std::vector<TBmid<TBfastmode>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            MIDs[MID] = reader.readFastmode(files[MID]);
        }
        anEvtFromData->setTCB(MIDs[0].evt());
        anEvtFromData->set(MIDs);
        auto dataFromData = anEvtFromData->data(cid);
        // Reading ntuple
        fileChain->GetEntry(iEvt);
        auto dataFromNtuple = anEvtFromNtuple->data(cid);

        // Compare contents
        if (dataFromData.adc() != dataFromNtuple.adc()) {
            errorCount++;
            std::cout << "[ERROR] ADC from data does not match with ntuple" << std::endl;
            std::cout << "Evt : " << iEvt << " Data ADC : " << dataFromData.adc() << " Ntuple ADC : " << dataFromNtuple.adc() << std::endl;
        }

        printProgress( (iEvt + 1) , totalEntryfromData);
    }

    if (errorCount) {
        std::cout << "[ERROR] Event content mismatch error in MID : " << cid.mid() << " Ch : " << cid.channel() << std::endl;
        return false;
    }
    else {
        std::cout << "[OK] Event content well matched in MID : " << cid.mid() << " Ch : " << cid.channel() << std::endl;
        return true;
    }
}

bool TBvalid::py_simpleValidWave(TBcid cid) {
    if (this->datList_.size() == 0) {
        std::cout << "No data list found, please use setDataList() first" << std::endl;
        return false;
    }
    if (this->ntupleList_.size() == 0) {
        std::cout << "No ntuple list found, please use setNtupleList() first" << std::endl;
        return false;
    }
    bool valid = simpleValidWave(this->datList_, this->ntupleList_, cid);
    return valid;
}

bool TBvalid::simpleValidWave(const std::vector<std::vector<std::string>>& datList,const std::vector<std::string>& ntupleList, TBcid cid) {
    std::cout << "Quick waveform validation..." << std::endl;
    const int numOfMID = 15;
    const int numOfFiles = datList.size();

    // Counting total entries from data file
    std::vector<int> entryPerMID(numOfMID);
    std::vector< std::vector<int> > entryPerFile( numOfFiles, std::vector<int>(numOfMID) );
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        for (unsigned MID = 0; MID < numOfMID; MID++) {
        FILE* fp = fopen(datList[fileNum][MID].c_str(), "rb");
        fseek(fp, 0L, SEEK_END);
        unsigned long long file_size = ftell(fp);
        fclose(fp);
        entryPerMID[MID] += (int)(file_size / 65536);
        entryPerFile[fileNum][MID] = (int)(file_size / 65536);
        }
    }
    int totalEntryfromData = *std::min_element(entryPerMID.begin(), entryPerMID.end());

    // Set TChain, count total entries from ntuple
    TChain* fileChain = new TChain("events");
    for (std::string fName : ntupleList) {
        fileChain->Add(fName.c_str());
    }
    int totalEntryfromNtuple = fileChain->GetEntries();

    // Error if total entry is different
    if (totalEntryfromData != totalEntryfromNtuple) {
        std::cout << "[ERROR] Ntuple total entry : " << totalEntryfromNtuple << " != Data total entry : " << totalEntryfromData << std::endl;
    }

    std::cout << "Total Entry : " << totalEntryfromData << std::endl;

    // Open data files and prepare reader
    TBevt<TBwaveform>* anEvtFromData = new TBevt<TBwaveform>();
    std::vector<int> currentOpenFileNum(numOfMID);
    std::vector<int> entryCounted(numOfMID);
    std::vector<FILE*> files(numOfMID);
    for (unsigned MID = 0; MID < numOfMID; MID++) {
        files[MID] = fopen(datList[0][MID].c_str(), "rb");
    }
    TBread reader = TBread();

    // Opening & reading ntuple with TChain
    TBevt<TBwaveform>* anEvtFromNtuple = new TBevt<TBwaveform>();
    fileChain->SetBranchAddress("TBevt", &anEvtFromNtuple);
    
    // For error count
    int errorCount = 0;

    // Start evt loop
    for (unsigned iEvt = 0; iEvt < totalEntryfromData; iEvt++){
        // Reading data
        // if count[MID] exceeds # of entry of corresponding file, move to next file
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            // currently opened file number for each MID
            if ( entryPerFile[currentOpenFileNum[MID]][MID] == entryCounted[MID] ) {
                fclose(files[MID]);
                currentOpenFileNum[MID]++;
                files[MID] = fopen(datList[currentOpenFileNum[MID]][MID].c_str(), "rb");
                entryCounted[MID] = 0;
            }
        }
        std::vector<TBmid<TBwaveform>> MIDs(numOfMID);
        for (unsigned MID = 0; MID < numOfMID; MID++) {
            MIDs[MID] = reader.readWaveform(files[MID]);
        }
        anEvtFromData->setTCB(MIDs[0].evt());
        anEvtFromData->set(MIDs);
        auto dataFromData = anEvtFromData->data(cid);
        // Reading ntuple
        fileChain->GetEntry(iEvt);
        auto dataFromNtuple = anEvtFromNtuple->data(cid);

        std::vector<short> waveFromData = dataFromData.waveform();
        std::vector<short> waveFromNtuple = dataFromNtuple.waveform();
        
        for (int i = 1; i < 1001; i++) {
            if (dataFromData.waveform()[i] != dataFromNtuple.waveform()[i]) {
                errorCount++;
                std::cout << "[ERROR] Waveform from data does not match with ntuple Evt : " << iEvt << " bin : " << i << std::endl;
                std::cout << "Data wave : " << dataFromData.waveform()[i] << std::endl;
                std::cout << "Ntuple wave : " << dataFromNtuple.waveform()[i] << std::endl;
            }
        }

        printProgress( (iEvt + 1) , totalEntryfromData);
    }

    if (errorCount) {
        std::cout << "[ERROR] Event content mismatch error in MID : " << cid.mid() << " Ch : " << cid.channel() << std::endl;
        return false;
    }
    else {
        std::cout << "[OK] Event content well matched in MID : " << cid.mid() << " Ch : " << cid.channel() << std::endl;
        return true;
    }
}

std::vector<std::vector<std::string>> TBvalid::readDataList(const boost::python::list& py_datList) {
    const int numOfFiles = boost::python::len(py_datList);

    std::vector< std::vector<std::string> > totalFiles(numOfFiles, std::vector<std::string>(boost::python::len(py_datList[0])));
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        for (unsigned MID = 0; MID < boost::python::len(py_datList[fileNum]); MID++) {
            totalFiles[fileNum][MID] = boost::python::extract<std::string>(py_datList[fileNum][MID]);
        }
    }

    return totalFiles;
}

std::vector<std::string> TBvalid::readNtupleList(const boost::python::list& py_ntupleList) {
    const int numOfFiles = boost::python::len(py_ntupleList);

    std::vector<std::string> totalFiles(numOfFiles);
    for (unsigned fileNum = 0; fileNum < numOfFiles; fileNum++) {
        totalFiles[fileNum] = boost::python::extract<std::string>(py_ntupleList[fileNum]);
    }

    return totalFiles;
}

void TBvalid::printProgress(const int currentStep, const int totalStep) {
    // print progress
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