#include <boost/test/unit_test.hpp>

#include "modelsettings.hpp"
#include "xmlmodelfile.hpp"
#include "seismic_parameters.hpp"
#include "seismic_regridding.hpp"
#include "seismic_forward.hpp"

#include "nrlib/segy/segytrace.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/segy/traceheader.hpp"

#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/iotools/logkit.hpp"

#include <iostream>

BOOST_AUTO_TEST_CASE( test_nmo_pp_seis )
{

  NRLib::LogKit::SetScreenLog(NRLib::LogKit::L_Low);
  NRLib::LogKit::StartBuffering();

  bool failedModelFile = false;
  //std::string data_dir;
  //if (boost::unit_test::framework::master_test_suite().argc > 1)
  //  data_dir = boost::unit_test::framework::master_test_suite().argv[1];
  //std::string inputfile = NRLib::PrependDir(data_dir, "modelfile_boost_test.xml");
  std::string inputfile = "modelfile_nmo_pp.xml";
  XmlModelFile modelFile(inputfile);
  ModelSettings *model_settings = modelFile.getModelSettings();
  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
  if (!failedModelFile) {
    SeismicParameters seismic_parameters = SeismicParameters(model_settings);
    SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
    SeismicForward::DoSeismicForward(seismic_parameters);
  }


  //compare traces
  std::string filename = "answ_nmo_pp_seismic_time.segy";
  NRLib::TraceHeaderFormat thf(2);
  NRLib::SegY segy_answ(filename, 0, thf);
  filename = "nmo_pp_seismic_time.segy";
  NRLib::SegY segy_test(filename, 0, thf);
  for (size_t tr = 0; tr < 4; ++tr){
    NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
    NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
    float value_answ, value_test;
    for (size_t i = 1; i < 126; ++i){
      value_answ = trace_answ->GetValue(i);
      value_test = trace_test->GetValue(i);
      BOOST_CHECK_CLOSE(value_answ, value_test, 1);
    }
  }

  if (true) {
    filename = "answ_nmo_pp_seismic_timeshift.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_pp_seismic_timeshift.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_nmo_pp_seismic_depth.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_pp_seismic_depth.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_nmo_pp_seismic_prenmo_time.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_pp_seismic_prenmo_time.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( test_nmo_ps_seis )
{
  bool failedModelFile = false;
  std::string inputfile = "modelfile_nmo_ps.xml";
  XmlModelFile modelFile(inputfile);
  ModelSettings *model_settings = modelFile.getModelSettings();
  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
  if (!failedModelFile) {
    SeismicParameters seismic_parameters = SeismicParameters(model_settings);
    SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
    SeismicForward::DoSeismicForward(seismic_parameters);
  }


  //compare traces
  std::string filename = "answ_nmo_ps_seismic_time.segy";
  NRLib::TraceHeaderFormat thf(2);
  NRLib::SegY segy_answ(filename, 0, thf);
  filename = "nmo_ps_seismic_time.segy";
  NRLib::SegY segy_test(filename, 0, thf);
  for (size_t tr = 0; tr < 4; ++tr){
    NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
    NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
    float value_answ, value_test;
    for (size_t i = 1; i < 126; ++i){
      value_answ = trace_answ->GetValue(i);
      value_test = trace_test->GetValue(i);
      BOOST_CHECK_CLOSE(value_answ, value_test, 1);
    }
  }

  if (true) {
    filename = "answ_nmo_ps_seismic_timeshift.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_ps_seismic_timeshift.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_nmo_ps_seismic_depth.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_ps_seismic_depth.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_nmo_ps_seismic_prenmo_time.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "nmo_ps_seismic_prenmo_time.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
}


BOOST_AUTO_TEST_CASE( test_pp_seis ) 
{

  bool failedModelFile = false;
  std::string inputfile = "modelfile_pp.xml";
  XmlModelFile modelFile(inputfile);
  ModelSettings *model_settings = modelFile.getModelSettings();
  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
  if (!failedModelFile) {
    SeismicParameters seismic_parameters = SeismicParameters(model_settings);
    SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
    SeismicForward::DoSeismicForward(seismic_parameters);
  }


  //compare traces
  std::string filename = "answ_pp_seismic_time.segy";
  NRLib::TraceHeaderFormat thf(2);
  NRLib::SegY segy_answ(filename, 0, thf);
  filename = "pp_seismic_time.segy";
  NRLib::SegY segy_test(filename, 0, thf);
  for (size_t tr = 0; tr < 4; ++tr){
    NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
    NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
    float value_answ, value_test;
    for (size_t i = 1; i < 126; ++i){
      value_answ = trace_answ->GetValue(i);
      value_test = trace_test->GetValue(i);
      BOOST_CHECK_CLOSE(value_answ, value_test, 1);
    }
  }

  if (true) {
    filename = "answ_pp_seismic_timeshift.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "pp_seismic_timeshift.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_pp_seismic_depth.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "pp_seismic_depth.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( test_ps_seis ) 
{
  bool failedModelFile = false;
  std::string inputfile = "modelfile_ps.xml";
  XmlModelFile modelFile(inputfile);
  ModelSettings *model_settings = modelFile.getModelSettings();
  if (modelFile.getParsingFailed()) {
    failedModelFile = true;
  }
  if (!failedModelFile) {
    SeismicParameters seismic_parameters = SeismicParameters(model_settings);
    SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
    SeismicForward::DoSeismicForward(seismic_parameters);
  }


  //compare traces
  std::string filename = "answ_ps_seismic_time.segy";
  NRLib::TraceHeaderFormat thf(2);
  NRLib::SegY segy_answ(filename, 0, thf);
  filename = "ps_seismic_time.segy";
  NRLib::SegY segy_test(filename, 0, thf);
  for (size_t tr = 0; tr < 4; ++tr){
    NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
    NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
    float value_answ, value_test;
    for (size_t i = 1; i < 126; ++i){
      value_answ = trace_answ->GetValue(i);
      value_test = trace_test->GetValue(i);
      BOOST_CHECK_CLOSE(value_answ, value_test, 1);
    }
  }

  if (true) {
    filename = "answ_ps_seismic_timeshift.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "ps_seismic_timeshift.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  if (true) {
    filename = "answ_ps_seismic_depth.segy";
    NRLib::SegY segy_answ(filename, 0, thf);
    filename = "ps_seismic_depth.segy";
    NRLib::SegY segy_test(filename, 0, thf);
    for (size_t tr = 0; tr < 4; ++tr){
      NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
      NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
      float value_answ, value_test;
      for (size_t i = 1; i < 126; ++i){
        value_answ = trace_answ->GetValue(i);
        value_test = trace_test->GetValue(i);
        BOOST_CHECK_CLOSE(value_answ, value_test, 1);
      }
    }
  }
  NRLib::LogKit::EndLog();
}