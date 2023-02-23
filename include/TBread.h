#ifndef TBread_h
#define TBread_h 1

#include "TBmid.h"
#include <string>

#include <boost/python.hpp>

class TBread {
public:
  TBread() {}
  ~TBread() {}

  void ntuplizeWaveform(const boost::python::list& alist, const std::string& output, const int maxEntry, const int entryPerSplit);
  void ntuplizeFastmode(const boost::python::list& alist, const std::string& output, const int maxEntry, const int entryPerSplit);

  std::string mappingpath() const { return mappingpath_; }
  void setMappingPath(const std::string& path) { mappingpath_ = path; }

  std::string pedestalpath() const { return pedestalpath_; }
  void setPedestalPath(const std::string& path) { pedestalpath_ = path; }

  FILE* py_readFile(char* filename);
  TBmid<TBwaveform> readWaveform(FILE* fp);
  TBmid<TBfastmode> readFastmode(FILE* fp);
private:
  TBmidbase readMetadata(FILE* fp);
  void printProgress(const int currentStep, const int progressStep);

  std::string mappingpath_;
  std::string pedestalpath_;
};

#endif
