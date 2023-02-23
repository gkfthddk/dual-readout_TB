#include "TButility.h"

#include <iostream>
#include <stdexcept>
#include <fstream>

TBdetector::detid TButility::detid(int tid) const {
  if ( static_cast<int>(TBdetector::detid::DWC1analogue) <= tid &&
       static_cast<int>(TBdetector::detid::DWC1digital) > tid )
    return TBdetector::detid::DWC1analogue;
  if ( static_cast<int>(TBdetector::detid::DWC1digital) <= tid &&
       static_cast<int>(TBdetector::detid::DWC2analogue) > tid )
    return TBdetector::detid::DWC1digital;
  if ( static_cast<int>(TBdetector::detid::DWC2analogue) <= tid &&
       static_cast<int>(TBdetector::detid::DWC2digital) > tid )
    return TBdetector::detid::DWC2analogue;
  if ( static_cast<int>(TBdetector::detid::DWC2digital) <= tid &&
       static_cast<int>(TBdetector::detid::extTrig1) > tid )
    return TBdetector::detid::DWC2digital;
  if ( static_cast<int>(TBdetector::detid::extTrig1) == tid )
    return TBdetector::detid::extTrig1;
  if ( static_cast<int>(TBdetector::detid::extTrig2) == tid )
    return TBdetector::detid::extTrig2;
  if ( static_cast<int>(TBdetector::detid::veto) == tid )
    return TBdetector::detid::veto;
  if ( static_cast<int>(TBdetector::detid::preshower) == tid )
    return TBdetector::detid::preshower;
  if ( static_cast<int>(TBdetector::detid::tail) == tid )
    return TBdetector::detid::tail;
  if ( static_cast<int>(TBdetector::detid::muon) == tid )
    return TBdetector::detid::muon;
  if ( static_cast<int>(TBdetector::detid::SiPM) == tid )
    return TBdetector::detid::SiPM;
  if ( static_cast<int>(TBdetector::detid::PMT) == tid )
    return TBdetector::detid::PMT;
  if ( static_cast<int>(TBdetector::detid::Module3D_C) == tid )
    return TBdetector::detid::Module3D_C;
  if ( static_cast<int>(TBdetector::detid::Module3D_S) == tid )
    return TBdetector::detid::Module3D_S;

  return TBdetector::detid::nulldet;
}

void TButility::loading(const std::string& path) {
  int mid, ch;
  int module, tower;
  int isCeren, isSiPM;
  int plate, column; 

  std::cout << "Loading mapping file : " << path << std::endl;

  TChain* mapChain = new TChain("mapping");
  mapChain->Add((TString)path);

  mapChain->SetBranchAddress("mid",&mid);
  mapChain->SetBranchAddress("ch",&ch);
  mapChain->SetBranchAddress("module",&module);
  mapChain->SetBranchAddress("tower",&tower);
  mapChain->SetBranchAddress("isCeren",&isCeren);
  mapChain->SetBranchAddress("isSiPM",&isSiPM);
  mapChain->SetBranchAddress("plate",&plate);
  mapChain->SetBranchAddress("column",&column);


  for ( int i = 0; i < mapChain->GetEntries(); i++ ) {
    mapChain->GetEntry(i);

    auto cid = TBcid(mid,ch);

    if ( module==3 ) {
      TBdetector::detid adetid = this->detid(tower);
      auto det = TBdetector(adetid);

      // column starts from 0 for all auxiliary detectors
      int col = tower - static_cast<int>(adetid);
      det.encodeColumn(col);

      mapping_.insert(std::make_pair(cid,det));

      continue;
    }

    if ( isSiPM ) {
      auto det = TBdetector(TBdetector::detid::SiPM);
      det.encodeModule(module,tower,static_cast<bool>(isCeren));
      det.encodeSiPM(plate,column);

      mapping_.insert(std::make_pair(cid,det));

      continue;
    }

    auto det = TBdetector(TBdetector::detid::PMT);
    det.encodeModule(module,tower,static_cast<bool>(isCeren));

    mapping_.insert(std::make_pair(cid,det));
  }

  delete mapChain;
}

TBdetector TButility::find(const TBcid& cid) const {
  if (mapping_.find(cid)==mapping_.end())
    return TBdetector(TBdetector::detid::nulldet);

  return mapping_.at(cid);
}

void TButility::loadped(const std::string& path) {

  TH2F* pedMap = (TH2F*)(TFile::Open((TString)path)->Get("pedestal"));

  float ped;

  for ( int mid = 1; mid <= 15; mid++ ) {
    for ( int ch = 1; ch <= 32; ch++ ) {

      ped = pedMap->GetBinContent(mid, ch);

      auto cid = TBcid(mid,ch);
      pedmap_.insert(std::make_pair(cid,ped));

    }
  }

  delete pedMap;
};

float TButility::retrievePed(const TBcid& cid) const {
  if (pedmap_.find(cid)==pedmap_.end()) {
    cid.print();
    throw std::runtime_error("TButility - cannot find pedestal of the TBcid!");
  }

  return pedmap_.at(cid);
}

int TButility::pid(float psadc, float muadc=0.) const {
  int pid = 0; // ped

  if ( psadc < PSpedcut_ )
    return pid;

  if ( psadc > PSpedcut_ && psadc < PS1mipcut_ ) {
    if( muadc < muoncut_ )
      pid = 211; // hadron
    else
      pid = 13; // muon
  }

  if ( psadc > PS3mipcut_ )
    pid = 11; // electron

  return pid;
}

TBcid TButility::getcid(TBdetector::detid did) const {

  for ( auto detInfo : mapping_ ) 
    if ( detInfo.second.det() == did )
      return detInfo.first;

  return TBcid(0, 0);
}

TBcid TButility::getcid(TBdetector::detid did, int module, int tower, bool isCeren) const {

  for ( auto detInfo : mapping_ ) 
    if ( detInfo.second.det() == did 
      && detInfo.second.module() == module
      && detInfo.second.tower() == tower
      && detInfo.second.isCeren() == isCeren )
      return detInfo.first;

  return TBcid(0, 0);
}

TBcid TButility::getcid(int module, int tower, bool isCeren) const {

  for ( auto detInfo : mapping_ ) 
    if ( detInfo.second.module() == module
      && detInfo.second.tower() == tower
      && detInfo.second.isCeren() == isCeren )
      return detInfo.first;
  
  return TBcid(0, 0);
}

TBcid TButility::getcid(int did, int module, int tower, bool isCeren) const {

    TBdetector::detid detid = TButility::detid(did);

    for ( auto detInfo : mapping_ ) 
    if ( detInfo.second.det() == detid 
      && detInfo.second.module() == module
      && detInfo.second.tower() == tower
      && detInfo.second.isCeren() == isCeren )
      return detInfo.first;

  return TBcid(0, 0);
}

TBcid TButility::getcid(TBdetector::detid did, int module, int tower, int column, int plate, bool isCeren) const {

  for ( auto detInfo : mapping_ ) 
    if ( detInfo.second.det() == did 
      && detInfo.second.module() == module
      && detInfo.second.tower() == tower
      && detInfo.second.column() == column
      && detInfo.second.plate() == plate
      && detInfo.second.isCeren() == isCeren )
      return detInfo.first;

  return TBcid(0, 0);
}