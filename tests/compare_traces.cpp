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
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_nmo_pp_seis_noise)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_nmo_pp_noise.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE( test_nmo_ps_seis )
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_nmo_ps.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_off_pp_seis)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_off_pp.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 3;
    std::string prefix = "off_pp";
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_off_ps_seis)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_off_ps.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
    }

    size_t n_traces = 16, n_samples = 126, n_output = 3;
    std::string prefix = "off_ps";
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
    std::cout << "\n";
  }
}
BOOST_AUTO_TEST_CASE( test_pp_seis )
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_pp.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE( test_ps_seis )
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_ps.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_ps_seis_noise)
{
  if (true) {
    bool check_small = true;
    double check_value = 1e-4;

    std::string inputfile = "modelfile_ps_noise.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);
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
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_rem_negz)
{
  if (true) {
    std::string inputfile = "modelfile_rem_negz.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);

      size_t n_traces = 3, n_samples = 154, n_output = 1;
      std::string prefix = "rem_negz";
      std::vector<std::string> filenames(n_output);
      filenames[0] = prefix + "_seismic_time.segy";
      for (size_t fi = 0; fi < n_output; ++fi) {
        NRLib::TraceHeaderFormat thf(2);
        NRLib::SegY segy_test(filenames[fi], 0, thf);
        NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
        //for (size_t tr = 0; tr < n_traces; ++tr) {
        size_t tr = 0;
        NRLib::SegYTrace *trace_answ;
        NRLib::SegYTrace *trace_test;
        while (tr < 1243) {
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
          tr++;
        }
        float value_answ, value_test;
        std::vector<double> answ_vec(n_samples);
        for (size_t j = 0; j < n_traces; ++j) {
          for (size_t i = 1; i < n_samples; ++i) {
            value_answ = trace_answ->GetValue(i);
            answ_vec[i] = value_answ;
            value_test = trace_test->GetValue(i);
            BOOST_CHECK_CLOSE(value_answ, value_test, 1);
          }
          //seismic_parameters.GetSeismicOutput()->PrintVector(answ_vec, NRLib::ToString(j) + "_answ_vec.txt");
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
        }
        n_samples = 259;
      }
    }
    std::cout << "\n";
  }
}

BOOST_AUTO_TEST_CASE(test_keep_negz)
{
  if (true) {
    std::string inputfile = "modelfile_keep_negz.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);

      size_t n_traces = 3, n_samples = 154, n_output = 1;
      std::string prefix = "keep_negz";
      std::vector<std::string> filenames(n_output);
      filenames[0] = prefix + "_seismic_time.segy";
      for (size_t fi = 0; fi < n_output; ++fi) {
        NRLib::TraceHeaderFormat thf(2);
        NRLib::SegY segy_test(filenames[fi], 0, thf);
        NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
        //for (size_t tr = 0; tr < n_traces; ++tr) {
        size_t tr = 0;
        NRLib::SegYTrace *trace_answ;
        NRLib::SegYTrace *trace_test;
        while (tr < 1243) {
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
          tr++;
        }
        float value_answ, value_test;
        std::vector<double> answ_vec(n_samples);
        for (size_t j = 0; j < n_traces; ++j) {
          for (size_t i = 1; i < n_samples; ++i) {
            value_answ = trace_answ->GetValue(i);
            answ_vec[i] = value_answ;
            value_test = trace_test->GetValue(i);
            BOOST_CHECK_CLOSE(value_answ, value_test, 1);
          }
          //seismic_parameters.GetSeismicOutput()->PrintVector(answ_vec, NRLib::ToString(j) + "_answ_vec.txt");
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
        }
        n_samples = 259;
      }
    }
    std::cout << "\n";
  }
}


BOOST_AUTO_TEST_CASE(test_cornerpt)
{
  if (true) {
    std::string inputfile = "modelfile_cornerpt.xml";
    std::cout << "Modelfile: " << inputfile << "\n";
    XmlModelFile modelFile(inputfile);
    ModelSettings *model_settings = modelFile.getModelSettings();
    if (modelFile.getParsingFailed() == false) {
      SeismicParameters seismic_parameters = SeismicParameters(model_settings);
      SeismicRegridding::MakeSeismicRegridding(seismic_parameters, model_settings, 1);
      SeismicForward::DoSeismicForward(seismic_parameters, model_settings);

      size_t n_traces = 3, n_samples = 154, n_output = 1;
      std::string prefix = "cornerpt";
      std::vector<std::string> filenames(n_output);
      filenames[0] = prefix + "_seismic_time.segy";
      for (size_t fi = 0; fi < n_output; ++fi) {
        NRLib::TraceHeaderFormat thf(2);
        NRLib::SegY segy_test(filenames[fi], 0, thf);
        NRLib::SegY segy_answ("answ_" + filenames[fi], 0, thf);
        //for (size_t tr = 0; tr < n_traces; ++tr) {
        size_t tr = 0;
        NRLib::SegYTrace *trace_answ;
        NRLib::SegYTrace *trace_test;
        while (tr < 1243) {
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
          tr++;
        }
        float value_answ, value_test;
        std::vector<double> answ_vec(n_samples);
        for (size_t j = 0; j < n_traces; ++j) {
          for (size_t i = 1; i < n_samples; ++i) {
            value_answ = trace_answ->GetValue(i);
            answ_vec[i] = value_answ;
            value_test = trace_test->GetValue(i);
            BOOST_CHECK_CLOSE(value_answ, value_test, 1);
          }
          //seismic_parameters.GetSeismicOutput()->PrintVector(answ_vec, NRLib::ToString(j) + "_answ_vec.txt");
          trace_answ = segy_answ.GetNextTrace();
          trace_test = segy_test.GetNextTrace();
        }
        n_samples = 259;
      }
      NRLib::LogKit::EndLog();
    }
    std::cout << "\n";
  }
}
