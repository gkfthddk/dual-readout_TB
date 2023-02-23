#include "TBread.h"
#include "TBdetector.h"
#include "TButility.h"
#include "TBmonit.h"
#include "TBplot.h"
#include "TBvalid.h"
#include "TBmid.h"

#include <TH1.h>

#include <boost/python.hpp>

BOOST_PYTHON_MODULE(pydrcTB) {
  boost::python::class_<TBread>("TBread")
    .def("ntuplizeWaveform", &TBread::ntuplizeWaveform)
    .def("ntuplizeFastmode", &TBread::ntuplizeFastmode)
    .def("readWaveform", &TBread::readWaveform)
    .def("readFastmode", &TBread::readFastmode)
    .def("setMappingPath", &TBread::setMappingPath)
    .def("setPedestalPath", &TBread::setPedestalPath);

  boost::python::class_<TH1F>("TH1F"); // For returning TH1F
  boost::python::class_<TH1D>("TH1D"); // For returning TH1D

  boost::python::class_<TBvalid>("TBvalid")
    .def("drawRatio", &TBvalid::drawRatio)
    .def("checkRatio", &TBvalid::checkRatio)
    .def("drawFastHistFromData", static_cast<TH1F* (TBvalid::*)(TBcid cid, const std::string histName)>(&TBvalid::drawFastHistFromData), boost::python::return_internal_reference<>())
    .def("drawFastHistFromNtuple", static_cast<TH1F* (TBvalid::*)(TBcid cid, const std::string histName)>(&TBvalid::drawFastHistFromNtuple), boost::python::return_internal_reference<>())
    .def("drawWaveHistFromData", static_cast<TH1F* (TBvalid::*)(TBcid cid, const std::string histName)>(&TBvalid::drawWaveHistFromData), boost::python::return_internal_reference<>())
    .def("drawWaveHistFromNtuple", static_cast<TH1F* (TBvalid::*)(TBcid cid, const std::string histName)>(&TBvalid::drawWaveHistFromNtuple), boost::python::return_internal_reference<>())
    .def("simpleValidFast", &TBvalid::py_simpleValidFast)
    .def("simpleValidWave", &TBvalid::py_simpleValidWave)
    .def("checkFastTrigNum", &TBvalid::checkFastTrigNum)
    .def("checkWaveTrigNum", &TBvalid::checkWaveTrigNum)
    .def("setDataList", &TBvalid::setDataList)
    .def("setNtupleList", &TBvalid::setNtupleList);

  boost::python::class_<TBcid>("TBcid", boost::python::init<int, int>())
    .def("mid", &TBcid::mid)
    .def("channel", &TBcid::channel)
    .def("print", &TBcid::print);

  boost::python::class_<TBdetector>("TBdetector")
    .def("detType", &TBdetector::detType)
    .def("isSiPM", &TBdetector::isSiPM)
    .def("isModule", &TBdetector::isModule)
    .def("isNull", &TBdetector::isNull)
    .def("module", &TBdetector::module)
    .def("tower", &TBdetector::tower)
    .def("isCeren", &TBdetector::isCeren)
    .def("plate", &TBdetector::plate)
    .def("column", &TBdetector::column);

  boost::python::class_<TButility>("TButility")
    .def("loading", &TButility::loading)
    .def("loadped", &TButility::loadped)
    .def("find", &TButility::find)
    .def("retrievePed", &TButility::retrievePed)
    .def("setPSpedcut", &TButility::setPSpedcut)
    .def("setPS1mipcut", &TButility::setPS1mipcut)
    .def("setPS3mipcut", &TButility::setPS3mipcut)
    .def("setMuoncut", &TButility::setMuoncut)
    .def("getcid", static_cast<TBcid (TButility::*)(TBdetector::detid did) const>(&TButility::getcid))
    .def("getcid", static_cast<TBcid (TButility::*)(TBdetector::detid did, int module, int tower, bool isCeren) const>(&TButility::getcid))
    .def("getcid", static_cast<TBcid (TButility::*)(int module, int tower, bool isCeren) const>(&TButility::getcid))
    .def("getcid", static_cast<TBcid (TButility::*)(int did, int module, int tower, bool isCeren) const>(&TButility::getcid))
    .def("getcid", static_cast<TBcid (TButility::*)(TBdetector::detid did, int module, int tower, int column, int plate, bool isCeren) const>(&TButility::getcid))
    .def("pid", &TButility::pid);

  boost::python::class_<TBmonit>("TBmonit", boost::python::init<const std::string&>())
    .def("setFastmodeFiles", &TBmonit::setFastmodeFiles)
    .def("setWaveformFiles", &TBmonit::setWaveformFiles)
    .def("setPedestalPath", &TBmonit::setPedestalPath)
    .def("setMappingPath", &TBmonit::setMappingPath)
    .def("SetADCmax", &TBmonit::SetADCmax)
    .def("SetADCbin", &TBmonit::SetADCbin)
    .def("SetOutputName", &TBmonit::SetOutputName)
    .def("MonitPlots", &TBmonit::MonitPlots);

  boost::python::class_<TBplot>("TBplot", boost::python::init<int, int, const std::string&, const std::string&>())
    .def("openFile", &TBplot::openFile)
    .def("closeFile", &TBplot::closeFile)
    .def("loadTH", &TBplot::loadTH)
    .def("loadTH1D", &TBplot::loadTH1D)
    .def("loadTH2D", &TBplot::loadTH2D)
    .def("Draw", &TBplot::Draw)
    .def("setSaveName", &TBplot::setSaveName)
    .def("SavePng", &TBplot::SavePng);
}
