#ifndef TButility_h
#define TButility_h 1

#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include "TBdetector.h"

#include <TChain.h>
#include <TString.h>
#include <TH2.h>
#include <TFile.h>

class TButility {
public:
  TButility() {}
  ~TButility() {}

  void loading(const std::string& path);
  void loadped(const std::string& path);

  TBdetector::detid detid(int tid) const;
  TBdetector find(const TBcid& cid) const;
  float retrievePed(const TBcid& cid) const;
  int pid(float psadc, float muadc) const;

  TBcid getcid(TBdetector::detid did) const;
  TBcid getcid(TBdetector::detid did, int module, int tower, bool isCeren) const;

  TBcid getcid(int did, int module, int tower, bool isCeren) const;
  TBcid getcid(int module, int tower, bool isCeren) const;

  TBcid getcid(TBdetector::detid did, int module, int tower, int column, int plate, bool isCeren) const;

  void setPSpedcut(float in) { PSpedcut_ = in; }
  void setPS1mipcut(float in) { PS1mipcut_ = in; }
  void setPS3mipcut(float in) { PS3mipcut_ = in; }
  void setMuoncut(float in) { muoncut_ = in; }

private:
  std::map<TBcid, TBdetector> mapping_;
  std::map<TBcid, float> pedmap_;

  float PSpedcut_;
  float PS1mipcut_;
  float PS3mipcut_;
  float muoncut_;
};

#endif
