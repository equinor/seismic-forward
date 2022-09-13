/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "nr/nrlib/stormgrid/stormcontgrid.hpp"
#include "nr/nrlib/iotools/fileio.hpp"

#include <stdlib.h>

#include <iostream>
#include <string>
#include <cmath>

//----------------------------------------------------------------
NRLib::StormContGrid * readStormGrid(const std::string & filename)
//----------------------------------------------------------------
{
  // Check that file is present

  if (!NRLib::FileExists(filename)) {
    std::cout << "Cannot find the STORM file : \'" << filename << "\'" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    exit(1);
  }
  // Read and create grid

  NRLib::StormContGrid * grid = NULL;
  try {
    grid = new NRLib::StormContGrid(filename);
  }
  catch (NRLib::Exception & e) {
    std::cout << e.what() << std::endl;
    std::cout << "The file \'" << filename << "\' is probably not a STORM grid file" << std::endl;
    std::cout << "Aborting ...\n" << std::endl;
    exit(1);
  }

  return grid;
}

//---------------------------------------------------
bool compareGridDefinitions(NRLib::StormContGrid * a,
                            NRLib::StormContGrid * b)
//---------------------------------------------------
{
  if (a->GetNI()    == b->GetNI()   &&
      a->GetNJ()    == b->GetNJ()   &&
      a->GetNK()    == b->GetNK()   &&
      a->GetXMin()  == b->GetXMin() &&
      a->GetYMin()  == b->GetYMin() &&
      a->GetLX()    == b->GetLX()   &&
      a->GetLY()    == b->GetLY()   &&
      a->GetAngle() == b->GetAngle())
    return true;
  else
    return false;
}

//-----------------------------------------------------------
void compareGrids(NRLib::StormContGrid * a,
                  NRLib::StormContGrid * b,
                  float                & largest_difference,
                  float                & mean_abs_difference,
                  float                & mean_val_difference,
                  float                & mean_val_both_cubes)
//-----------------------------------------------------------
{
  size_t n = a->GetN();

  for (size_t i = 0; i < n ; i++) {
    float ai       = std::abs((*a)(i)); // Abs value to avoid seismic amplitudes to cancel in sum
    float bi       = std::abs((*b)(i)); // Abs value to avoid seismic amplitudes to cancel in sum
    float diff     = ai - bi;
    float abs_diff = std::abs(diff);

    mean_val_difference += diff;
    mean_abs_difference += abs_diff;
    mean_val_both_cubes += ai + bi;

    if (abs_diff > largest_difference) {
      largest_difference = abs_diff;
    }
  }
  mean_val_difference /= n;
  mean_abs_difference /= n;
  mean_val_both_cubes /= 2*n;
}

//------------------------------------------------------------------
void writeDifferencesToFile(const std::string & filename,
                            bool                equal,
                            float               largest_difference,
                            float               mean_abs_difference,
                            float               mean_val_difference,
                            float               mean_val_both_cubes)
//------------------------------------------------------------------
{
  std::ofstream outfile;
  NRLib::OpenWrite(outfile, filename);
  outfile
    << equal               << " "
    << largest_difference  << " "
    << mean_abs_difference << " "
    << mean_val_difference << " "
    << mean_val_both_cubes
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

  if (argc != 3) {
    std::cout << "Usage: ./compare answer-file output-file\n" << std::endl;
    exit(1);
  }
  std::string answer_file = std::string(argv[1]);
  std::string output_file = std::string(argv[2]);

  //
  // Read grids
  // ----------

  NRLib::StormContGrid * answer_grid = readStormGrid(answer_file);
  NRLib::StormContGrid * output_grid = readStormGrid(output_file);

  //
  // Check that grids are equal
  // --------------------------
  bool equal = compareGridDefinitions(answer_grid, output_grid);

  //
  // Compare grids
  // -------------

  float largest_difference  = 0.0f;
  float mean_abs_difference = 0.0f;
  float mean_val_difference = 0.0f;
  float mean_val_both_cubes = 0.0f;

  if (equal == true) {
    compareGrids(answer_grid,
                 output_grid,
                 largest_difference,
                 mean_abs_difference,
                 mean_val_difference,
                 mean_val_both_cubes);

    //
    // Write info to file for Perl import
    // ----------------------------------
  }
    std::string diff_file = "storm_volume_difference.txt";
    writeDifferencesToFile(diff_file,
                           equal,
                           largest_difference,
                           mean_abs_difference,
                           mean_val_difference,
                           mean_val_both_cubes);

  if (answer_grid != NULL)
    delete answer_grid;
  if (output_grid != NULL)
    delete output_grid;
}
