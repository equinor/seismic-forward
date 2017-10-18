/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "nr/nrlib/exception/exception.hpp"

#include "nr/nrlib/iotools/stringtools.hpp"
#include "nr/nrlib/iotools/logkit.hpp"

#include "nr/nrlib/segy/traceheader.hpp"
#include "nr/nrlib/segy/segygeometry.hpp"
#include "nr/nrlib/segy/segytrace.hpp"
#include "nr/nrlib/segy/segy.hpp"

#include <iostream>

//-------------------------------------------------------------------
void WriteDifferencesToFile(const std::string & filename,
                            int                 grid_defs_are_equal,
                            float               avg_diff,
                            float               bias,
                            float               max_amp,
                            float               max_diff,
                            int                 max_trace,
                            int                 max_sample)
//-------------------------------------------------------------------
{
  std::cout
    << "\n"
    << "equal_defs = " << grid_defs_are_equal << "\n"
    << "max_amp    = " << max_amp    << "\n"
    << "max_diff   = " << max_diff   << "\n"
    << "avg_diff   = " << avg_diff   << "\n"
    << "bias       = " << bias       << "\n"
    << "max_trace  = " << max_trace  << "\n"
    << "max_sample = " << max_sample << "\n"
    << std::endl;

  std::ofstream outfile;
  NRLib::OpenWrite(outfile, filename);
  outfile
    << grid_defs_are_equal << " "
    << max_amp             << " "
    << max_diff            << " "
    << avg_diff            << " "
    << bias                << " "
    << max_trace           << " "
    << max_sample
    << std::endl;

  outfile.close();
}

//--------------------------
void WriteComparisonFailed()
//--------------------------
{
  std::string diff_file           = "segy_amplitude_difference.txt";
  int         grid_defs_are_equal = -1;
  float       avg_diff            = -1;
  float       bias                = -1;
  float       max_amp             = -1;
  float       max_diff            = -1;
  int         max_trace           = -1;
  int         max_sample          = -1;

  WriteDifferencesToFile(diff_file,
                         grid_defs_are_equal,
                         avg_diff,
                         bias,
                         max_amp,
                         max_diff,
                         max_trace,
                         max_sample);
}

//--------------------------------------------------
NRLib::SegY * ReadSegY(const std::string & filename)
//--------------------------------------------------
{
  // Check that file is present
  if (!NRLib::FileExists(filename)) {
    std::cout << "Cannot find the SegY file : \'" << filename << "\'" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    WriteComparisonFailed();
    exit(1);
  }

  std::cout << "IN READSEGY" << std::endl;

  // Read and create grid

  NRLib::SegY * grid = NULL;
  try {
    // Because of a bug, a thf has to be specified when Segy is read. However,
    // all formats recognized by NRLib will be tried if SEISWORKS fails.
    NRLib::TraceHeaderFormat thf(NRLib::TraceHeaderFormat::SEISWORKS);
    double z0 = 0.0; // Dummy value. We are comparing two equal volumes
    grid = new NRLib::SegY(filename, z0, thf);
  }
  catch (NRLib::Exception & e) {
    std::cout << e.what() << std::endl;
    std::cout << "The file \'" << filename << "\' is probably not a SegY grid file" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    WriteComparisonFailed();
    exit(1);
  }
  return grid;
}



//------------------------------------------------------------------
NRLib::SegyGeometry * ReadSegyGeometry(const std::string & filename)
//------------------------------------------------------------------
{
  NRLib::SegyGeometry * geometry = NULL;

  try {
    // Because of a bug, a thf has to be specified when Segy is read. However,
    // all formats recognized by NRLib will be tried if SEISWORKS fails.
    const NRLib::TraceHeaderFormat thf(NRLib::TraceHeaderFormat::SEISWORKS);
    geometry = NRLib::SegY::FindGridGeometry(filename, &thf);
  }
  catch (NRLib::Exception & e) {
    std::cout << e.what() << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    WriteComparisonFailed();
    exit(1);
  }
  return geometry;
}

//-------------------------------------------------------
bool CompareGridDefinitions(NRLib::SegyGeometry * output,
                            NRLib::SegyGeometry * answer)
//-------------------------------------------------------
{
  std::cout << "ANSWER-CUBE\n"
            << "MinIL()  : " << answer->GetMinIL()  << "\n"
            << "MaxIL()  : " << answer->GetMaxIL()  << "\n"
            << "MinXL()  : " << answer->GetMinXL()  << "\n"
            << "MaxXL()  : " << answer->GetMaxXL()  << "\n"
            << "ILStep() : " << answer->GetILStep() << "\n"
            << "XLStep() : " << answer->GetXLStep() << "\n"
            << "X0()     : " << answer->GetX0()     << "\n"
            << "Y0()     : " << answer->GetY0()     << "\n"
            << "Nx()     : " << answer->GetNx()     << "\n"
            << "Ny()     : " << answer->GetNy()     << "\n"
            << "Dx()     : " << answer->GetDx()     << "\n"
            << "Dy()     : " << answer->GetDy()     << "\n"
            << "lx()     : " << answer->Getlx()     << "\n"
            << "ly()     : " << answer->Getly()     << "\n"
            << "Angle()  : " << answer->GetAngle()  << "\n"
            << std::endl;

  std::cout << "OUTPUT-CUBE\n"
            << "MinIL()  : " << output->GetMinIL()  << "\n"
            << "MaxIL()  : " << output->GetMaxIL()  << "\n"
            << "MinXL()  : " << output->GetMinXL()  << "\n"
            << "MaxXL()  : " << output->GetMaxXL()  << "\n"
            << "ILStep() : " << output->GetILStep() << "\n"
            << "XLStep() : " << output->GetXLStep() << "\n"
            << "X0()     : " << output->GetX0()     << "\n"
            << "Y0()     : " << output->GetY0()     << "\n"
            << "Nx()     : " << output->GetNx()     << "\n"
            << "Ny()     : " << output->GetNy()     << "\n"
            << "Dx()     : " << output->GetDx()     << "\n"
            << "Dy()     : " << output->GetDy()     << "\n"
            << "lx()     : " << output->Getlx()     << "\n"
            << "ly()     : " << output->Getly()     << "\n"
            << "Angle()  : " << output->GetAngle()  << "\n"
            << std::endl;

  if (answer->GetMinIL()  == output->GetMinIL()  &&
      answer->GetMaxIL()  == output->GetMaxIL()  &&
      answer->GetMinXL()  == output->GetMinXL()  &&
      answer->GetMaxXL()  == output->GetMaxXL()  &&
      answer->GetILStep() == output->GetILStep() &&
      answer->GetXLStep() == output->GetXLStep() &&
      answer->GetX0()     == output->GetX0()     &&
      answer->GetY0()     == output->GetY0()     &&
      answer->GetNx()     == output->GetNx()     &&
      answer->GetNy()     == output->GetNy()     &&
      answer->GetDx()     == output->GetDx()     &&
      answer->GetDy()     == output->GetDy()     &&
      answer->Getlx()     == output->Getlx()     &&
      answer->Getly()     == output->Getly()     &&
      answer->GetAngle()  == output->GetAngle())
    return true;
  else
    return false;
}

//--------------------------------------------
void CompareTraces(NRLib::SegY * segy_output,
                   NRLib::SegY * segy_answer,
                   float       & avg_abs_diff,
                   float       & bias,
                   float       & max_amp,
                   float       & max_diff,
                   int         & max_trace,
                   int         & max_sample)
//--------------------------------------------
{
  float  undef    = 99999;
  size_t n_tot    = 0;
  size_t n_traces = 0;

  try {
    n_traces = segy_answer->GetNTraces();
    std::cout << "Number of traces = " << n_traces << std::endl;
  }
  catch (NRLib::Exception & e) {
    std::cout << "Could not find number of traces in answer volume!\nAborting ...\n" << std::endl;
    WriteComparisonFailed();
    exit(1);
  }

  for (size_t t = 0 ; t < n_traces ; ++t) {

    NRLib::SegYTrace * trace_output = segy_output->GetNextTrace();
    NRLib::SegYTrace * trace_answer = segy_answer->GetNextTrace();

    if (trace_answer == NULL) {
      std::cout << "\nAnswer trace number " << t << " is null!\nAborting ...\n" << std::endl;
      WriteComparisonFailed();
      exit(1);
    }
    if (trace_output == NULL) {
      std::cout << "\nOutput trace number " << t << " is null!\nAborting ...\n" << std::endl;
      WriteComparisonFailed();
      exit(1);
    }

    if (t==0) {
      std::cout << "Number of samples = " << trace_answer->GetTrace().size() << std::endl;
    }

    for (size_t i = 0 ; i < trace_answer->GetTrace().size() ; ++i) {
      float amp_output = std::abs(trace_output->GetValue(i));   // Abs value to avoid seismic amplitudes to cancel
      float amp_answer = std::abs(trace_answer->GetValue(i));   // Abs value to avoid seismic amplitudes to cancel
      float abs_diff   = std::abs(amp_output - amp_answer);
      float diff       = amp_output - amp_answer;

      if (amp_answer != undef && amp_answer > max_amp) {
        max_amp = amp_answer;
      }
      if (abs_diff > max_diff) {
        max_diff   = abs_diff;
        max_trace  = t;
        max_sample = i;
      }
      avg_abs_diff += abs_diff;
      bias         += diff;
      n_tot        += 1;
    }
  }
  avg_abs_diff /= n_tot;
  bias         /= n_tot;
}

//-----------------------------
int main(int argc, char** argv)
//-----------------------------
{
  //
  // Read options
  // ------------
  //
  if (argc != 3) {
    std::cout << "Usage: ./compare answer-file output-file\n" << std::endl;
    WriteComparisonFailed();
    exit(1);
  }
  std::string file_answer = std::string(argv[1]);
  std::string file_output = std::string(argv[2]);

  //
  // Read grids (both volumes are assumed to have the same z0)
  // ---------------------------------------------------------
  //
  NRLib::SegY * segy_output = ReadSegY(file_output);
  NRLib::SegY * segy_answer = ReadSegY(file_answer);

  //
  // Check that grids are equal
  // --------------------------
  //
  NRLib::SegyGeometry * geometry_output = ReadSegyGeometry(file_output);
  NRLib::SegyGeometry * geometry_answer = ReadSegyGeometry(file_answer);

  int grid_defs_are_equal = CompareGridDefinitions(geometry_output,
                                                   geometry_answer);

  //
  // Compare traces
  // --------------
  //
  float  avg_diff   = 0.0;
  float  bias       = 0.0;
  float  max_amp    = 0.0;
  float  max_diff   = 0.0;
  int    max_trace  =  -1;
  int    max_sample =  -1;

  if (grid_defs_are_equal == 1) {
    CompareTraces(segy_output,
                  segy_answer,
                  avg_diff,
                  bias,
                  max_amp,
                  max_diff,
                  max_trace,
                  max_sample);
  }

  //
  // Write info to file for Perl import
  // ----------------------------------
  //
  std::string diff_file = "segy_amplitude_difference.txt";
  WriteDifferencesToFile(diff_file,
                         grid_defs_are_equal,
                         avg_diff,
                         bias,
                         max_amp,
                         max_diff,
                         max_trace,
                         max_sample);

  if (segy_answer != NULL)
    delete segy_answer;
  if (segy_output != NULL)
    delete segy_output;
}
