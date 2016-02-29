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
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;
    NRLib::LogKit::SetScreenLog(NRLib::LogKit::L_Low);
    NRLib::LogKit::StartBuffering();

    std::string inputfile = "modelfile_nmo_pp.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 4;
    std::string prefix = "nmo_pp";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";
    filenames[3] = prefix + "_seismic_time_prenmo.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_"+filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_nmo_pp_seis_noise)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_nmo_pp_noise.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 4;
    std::string prefix = "nmo_pp_noise";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";
    filenames[3] = prefix + "_seismic_time_prenmo.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( test_nmo_ps_seis )
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_nmo_ps.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 4;
    std::string prefix = "nmo_ps";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";
    filenames[3] = prefix + "_seismic_time_prenmo.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }

  }
}


BOOST_AUTO_TEST_CASE( test_pp_seis ) 
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_pp.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 3;
    std::string prefix = "pp";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( test_ps_seis )
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_ps.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 3;
    std::string prefix = "ps";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
  }
}


BOOST_AUTO_TEST_CASE(test_ps_seis_noise)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_ps_noise.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 3;
    std::string prefix = "ps_noise";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          if (check_small == false && value_answ < check_value && value_test < check_value) {
            value_test = value_answ;
          }
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
  }
}


BOOST_AUTO_TEST_CASE(test_snorre)
{
  if (false) {
    std::string inputfile = "modelfile_snorre_test.xml";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters);
      SeismicForward::DoSeismicForward(seismic_parameters);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 4;
    std::string prefix = "ps_noise";
    std::vector<std::string> filenames(n_output);
    filenames[0] = prefix + "_seismic_time.segy";
    filenames[1] = prefix + "_seismic_timeshift.segy";
    filenames[2] = prefix + "_seismic_depth.segy";
    filenames[3] = prefix + "_seismic_time_prenmo.segy";

    for (size_t fi = 0; fi < n_output; ++fi) {
      NRLib::TraceHeaderFormat thf(2);
      NRLib::SegY segy_test(filenames[fi], 0, thf);
      NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
      for (size_t tr = 0; tr < n_traces; ++tr) {
        NRLib::SegYTrace *trace_answ = segy_answ.GetNextTrace();
        NRLib::SegYTrace *trace_test = segy_test.GetNextTrace();
        float value_answ, value_test;
        for (size_t i = 1; i < n_samples; ++i) {
          value_answ = trace_answ->GetValue(i);
          value_test = trace_test->GetValue(i);
          BOOST_CHECK_CLOSE(value_answ, value_test, 1);
        }
      }
    }
    NRLib::LogKit::EndLog();
  }
}

