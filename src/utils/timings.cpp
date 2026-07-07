
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "utils/timings.hpp"

//---------------------------------------
void Timings::reportAll(double threshold)
//---------------------------------------
{
  NRLib::LogKit::WriteHeader("Timings Summary");

  calculateRest();

  if (c_total_ < 0.00001)
    c_total_ = 0.00001;
  if (w_total_ < 0.00001)
    w_total_ = 0.00001;

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSection                                          CPU Time/s            Real Time/s");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n----------------------------------------------------------------------------------");
  report("Loading corner-point grid              ", c_load_cornerpoint_grid_  , w_load_cornerpoint_grid_  , threshold);
  report("Make regular grid                      ", c_find_zvalues_           , w_find_zvalues_           , threshold);
  report("Resample elastic parameters            ", c_find_elastic_parameters_, w_find_elastic_parameters_, threshold);
  report("Forward modelling                      ", c_forward_modelling_      , w_forward_modelling_      , threshold);
  report("Write SegY files                       ", c_write_segy_             , w_write_segy_             , threshold);
  report("Dummy                                  ", c_dummy_                  , w_dummy_                  , threshold);
  report("Miscellaneous                          ", c_rest_                   , w_rest_                   , threshold);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n----------------------------------------------------------------------------------");
  report("Total                                  ", c_total_                  , w_total_                 , threshold);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
}

//-----------------------------
void Timings::reportTotal(void)
//-----------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nTotal CPU  time used by Seismic Forward: %8.2f s"  , c_total_);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nTotal real time used in Seismic Forward: %8.2f s\n", w_total_);
}

//-------------------------------------------------
void Timings::report(const std::string & text,
                     double              cpuThis,
                     double              wallThis,
                     double              threshold)
//-------------------------------------------------
{
  if (wallThis < 0.00001) // To omit stupit zero-treatment in ToString()
    wallThis = 0.00001;

  double percentCPU  = 100.0*cpuThis/c_total_;
  double percentWall = 100.0*wallThis/w_total_;

  if (percentCPU > 100.0)
    percentCPU = 100.0;
  if (percentWall > 100.0)
    percentWall = 100.0;

  if (cpuThis > threshold && percentCPU > threshold) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%s %9.2f   %6.2f%%    %9.2f   %6.2f%%", text.c_str(), cpuThis, percentCPU, wallThis, percentWall);
  }
}

void
Timings::calculateRest(void)
{
  c_rest_ = c_total_ - (c_load_cornerpoint_grid_
                        + c_find_zvalues_
                        + c_find_elastic_parameters_
                        + c_forward_modelling_
                        + c_write_segy_
                        + c_dummy_);
  w_rest_ = w_total_ - (w_load_cornerpoint_grid_
                        + w_find_zvalues_
                        + w_find_elastic_parameters_
                        + w_forward_modelling_
                        + w_write_segy_
                        + w_dummy_);
}

void
Timings::setTimeTotal(Timer & timer)
{
  c_total_ = timer.CPU();
  w_total_ = timer.Clock();
}

void
Timings::setTimeLoadCornerPointGrid(Timer & timer)
{
  c_load_cornerpoint_grid_ = timer.CPU();
  w_load_cornerpoint_grid_ = timer.Clock();
}

void
Timings::setTimeFindZValues(Timer & timer)
{
  c_find_zvalues_ = timer.CPU();
  w_find_zvalues_ = timer.Clock();
}

void
Timings::setTimeFindElasticParameters(Timer & timer)
{
  c_find_elastic_parameters_ = timer.CPU();
  w_find_elastic_parameters_ = timer.Clock();
}

void
Timings::setTimeForwardModelling(Timer & timer)
{
  c_forward_modelling_ = timer.CPU();
  w_forward_modelling_ = timer.Clock();
}

void
Timings::setTimeWriteSegy(Timer & timer)
{
  c_write_segy_ = timer.CPU();
  w_write_segy_ = timer.Clock();
}

void
Timings::addTimeDummy(Timer & timer)
{
  c_dummy_ += timer.CPU();
  w_dummy_ += timer.Clock();
}

double Timings::c_total_                   = 0.0;
double Timings::w_total_                   = 0.0;

double Timings::c_rest_                    = 0.0;
double Timings::w_rest_                    = 0.0;

double Timings::c_load_cornerpoint_grid_   = 0.0;
double Timings::w_load_cornerpoint_grid_   = 0.0;

double Timings::c_find_zvalues_            = 0.0;
double Timings::w_find_zvalues_            = 0.0;

double Timings::c_find_elastic_parameters_ = 0.0;
double Timings::w_find_elastic_parameters_ = 0.0;

double Timings::c_forward_modelling_       = 0.0;
double Timings::w_forward_modelling_       = 0.0;

double Timings::c_write_segy_              = 0.0;
double Timings::w_write_segy_              = 0.0;

double Timings::c_dummy_                   = 0.0;
double Timings::w_dummy_                   = 0.0;
