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

//--------------------------------------------------
NRLib::SegY * ReadSegY(const std::string & filename)
//--------------------------------------------------
{
  // Check that file is present
  if (!NRLib::FileExists(filename)) {
    std::cout << "Cannot find the SegY file : \'" << filename << "\'" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    exit(1);
  }

  // Read and create grid

  NRLib::SegY * grid = NULL;
  try {
    // Because of a bug, a thf has to be specified when Segy is read. However,
    // all formats recognized by NRLib will be tried if SEISWORKS fails.

    //NRLib::TraceHeaderFormat thf(NRLib::TraceHeaderFormat::SEISWORKS);

    double z0 = 0.0; // Dummy value. We are comparing two equal volumes
    grid = new NRLib::SegY(filename, z0);
  }

  catch (NRLib::Exception & e) {
    std::cout << e.what() << std::endl;
    std::cout << "The file \'" << filename << "\' is probably not a SegY grid file" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    exit(1);
  }
  return grid;
}

//--------------------------------------------------------
bool CompareGridDefinitions(const NRLib::SegyGeometry * a,
                            const NRLib::SegyGeometry * b)
//--------------------------------------------------------
{
  std::cout << "MinIL()  : " << a->GetMinIL()  << "\n"
            << "MaxIL()  : " << a->GetMaxIL()  << "\n"
            << "MinXL()  : " << a->GetMinXL()  << "\n"
            << "MaxXL()  : " << a->GetMaxXL()  << "\n"
            << "ILStep() : " << a->GetILStep() << "\n"
            << "XLStep() : " << a->GetXLStep() << "\n"
            << "X0()     : " << a->GetX0()     << "\n"
            << "Y0()     : " << a->GetY0()     << "\n"
            << "Nx()     : " << a->GetNx()     << "\n"
            << "Ny()     : " << a->GetNy()     << "\n"
            << "Dx()     : " << a->GetDx()     << "\n"
            << "Dy()     : " << a->GetDy()     << "\n"
            << "lx()     : " << a->Getlx()     << "\n"
            << "ly()     : " << a->Getly()     << "\n"
            << "Angle()  : " << a->GetAngle()  << "\n"
            << std::endl;

  if (a->GetMinIL()  == b->GetMinIL()  &&
      a->GetMaxIL()  == b->GetMaxIL()  &&
      a->GetMinXL()  == b->GetMinXL()  &&
      a->GetMaxXL()  == b->GetMaxXL()  &&
      a->GetILStep() == b->GetILStep() &&
      a->GetXLStep() == b->GetXLStep() &&
      a->GetX0()     == b->GetX0()     &&
      a->GetY0()     == b->GetY0()     &&
      a->GetNx()     == b->GetNx()     &&
      a->GetNy()     == b->GetNy()     &&
      a->GetDx()     == b->GetDx()     &&
      a->GetDy()     == b->GetDy()     &&
      a->Getlx()     == b->Getlx()     &&
      a->Getly()     == b->Getly()     &&
      a->GetAngle()  == b->GetAngle())
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
    std::cout << "n_traces = " << n_traces << std::endl;
  }
  catch (NRLib::Exception & e) {
    std::cout << "Could not find number of traces in answer volume!\nAborting ...\n" << std::endl;
    exit(1);
  }

  for (size_t t = 0 ; t < n_traces ; ++t) {

    NRLib::SegYTrace * trace_output = segy_output->GetNextTrace();
    NRLib::SegYTrace * trace_answer = segy_answer->GetNextTrace();

    if (trace_answer == NULL) {
      std::cout << "Answer trace number " << t << " is null!\nAborting ...\n" << std::endl;
      exit(1);
    }
    if (trace_output == NULL) {
      std::cout << "Output trace number " << t << " is null!\nAborting ...\n" << std::endl;
      exit(1);
    }

    std::cout << trace_answer->GetTrace().size() << std::endl;

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

//---------------------------------------------------------
void WriteDifferencesToFile(const std::string & filename,
                            float               avg_diff,
                            float               bias,
                            float               max_amp,
                            float               max_diff,
                            int                 max_trace,
                            int                 max_sample)
//---------------------------------------------------------
{
  std::cout
    << "max_amp    = " << max_amp    << "\n"
    << "max_diff   = " << max_diff   << "\n"
    << "avg_diff   = " << avg_diff   << "\n"
    << "bias       = " << bias       << "\n"
    << "max_trace  = " << max_trace  << "\n"
    << "max_sample = " << max_sample << "\n"
    << std::endl;

  int grid_defs_are_equal = 1; // Not implemented check for this yet ...

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

  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

  //
  // Check that grids are equal
  // --------------------------
  //

  /*
  NRLib::SegyGeometry * answer_geo = NULL;
  NRLib::SegyGeometry * output_geo = NULL;
  try {
    answer_geo = segy_answer->FindGridGeometry();
    output_geo = segy_output->FindGridGeometry();
  }
  catch (NRLib::Exception & e) {
    std::cout << e.what() << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    exit(1);
  }

  std::cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY" << std::endl;


  bool equal = CompareGridDefinitions(answer_geo, output_geo);
*/

  bool equal = true;

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

  CompareTraces(segy_output,
                segy_answer,
                avg_diff,
                bias,
                max_amp,
                max_diff,
                max_trace,
                max_sample);

  //
  // Write info to file for Perl import
  // ----------------------------------
  //
  std::string diff_file = "segy_amplitude_difference.txt";
  WriteDifferencesToFile(diff_file,
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
