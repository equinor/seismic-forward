// $Id: modelsettings.cpp 39 2014-02-06 08:53:13Z vigsnes $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// o  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// o  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "nrlib/math/mathutility.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "modelsettings.hpp"
#include "tasklist.hpp"

#include <iostream>
#include <cstdlib>

ModelSettings::ModelSettings(void)
  : constvp_(3),
    constvs_(3),
    constrho_(3),
    parameter_names_(3),
    output_segy_file_format_(NRLib::TraceHeaderFormat(NRLib::TraceHeaderFormat::SIP))
{
  log_file_name_                       = "Logfile";
  log_level_                           = NRLib::LogKit::L_Low;

  prefix_                              = "";
  suffix_                              = "";

  zero_thickness_limit_                = 0.1;
  twt_file_name_                       = "";
  traces_in_memory_                    = 100000;
  max_threads_                         = 100;

  add_noise_to_refl_coef_              = false;
  add_white_noise_                     = false;
  equal_noise_for_offsets_             = true;
  standard_deviation_1_                = 1.0;
  standard_deviation_2_                = 1.0;
  seed_1_                              = static_cast<unsigned long>(std::time(0)    );
  seed_2_                              = static_cast<unsigned long>(std::time(0) + 1);

  default_underburden_                 = false;
  ps_seismic_                          = false;
  nmo_corr_                            = false;

  theta_0_                             = 0.0;  // default values
  dtheta_                              = 0.0;
  theta_max_                           = 0.0;
  offset_0_                            = 0.0;
  doffset_                             = 0.0;
  offset_max_                          = 0.0;

  top_time_surface_                    = "";
  top_time_constant_                   = 1000.0;
  dx_                                  = 25.0;
  dy_                                  = 25.0;
  dz_                                  = 4.0;
  dt_                                  = 4.0;

  area_given_                          = false;
  area_from_surface_                   = "";
  area_from_segy_                      = "";

  il0_loc_                             = 189;
  xl0_loc_                             = 193;
  utmx_loc_                            = 181;
  utmy_loc_                            = 185;
  scalco_loc_                          =  71;
  start_time_loc_                      = 109;
  utm_precision_                       = -10;

  inline_start_                        = 0;
  xline_start_                         = 0;
  inline_direction_                    = "y";
  inline_step_                         = 1;
  xline_step_                          = 1;

  wavelet_scale_                       = 1.0;
  wavelet_length_                      = -999.0;
  wavelet_length_factor_               = 1.0;
  use_zero_time_from_header_           = false;
  z_wavelet_top_                       = 0.0;
  z_wavelet_bot_                       = 0.0;
  z_extrapol_factor_                   = 50.0;
  offset_without_stretch_              = false;

  time_window_specified_               = false;
  depth_window_specified_              = false;

  top_time_window_                     = -9999.0;
  bot_time_window_                     = -9999.0;
  top_depth_window_                    = -9999.0;
  bot_depth_window_                    = -9999.0;

  padding_factor_                      = 1.0;

  use_cornerpoint_interpol_            = false;
  cornerpoint_interpol_at_faults_      = false;
  use_fixed_triangularization_         = true;
  use_horizontal_interpolation_        = false;
  use_vertical_interpolation_          = true;
  use_bilinear_interpolation_          = false;
  use_active_pillars_                  = false;
  remove_negative_delta_z_             = false;
  elastic_parameters_time_segy_        = false;
  elastic_parameters_depth_segy_       = false;
  extra_parameters_time_segy_          = false;
  extra_parameters_depth_segy_         = false;
  seismic_stack_time_storm_            = false;
  seismic_stack_time_shift_storm_      = false;
  seismic_stack_depth_storm_           = false;
  seismic_stack_time_segy_             = false;
  seismic_stack_time_shift_segy_       = false;
  seismic_stack_depth_segy_            = false;

  v_w_                                 = 0.0;
  z_w_                                 = 0.0;

  resampl_param_to_segy_with_interpol_ = false;

  output_vp_                           = false;
  output_reflections_                  = false;
  output_zvalues_                      = false;
  output_seismic_time_                 = false;
  output_seismic_depth_                = false;
  output_seismic_timeshift_            = false;
  output_time_surfaces_                = false;
  output_depth_surfaces_               = false;
  output_wavelet_                      = false;
  output_twt_                          = false;
  output_vrms_                         = false;
  output_twt_offset_                   = false;
  output_time_segy_                    = false;
  output_depth_segy_                   = false;
  output_timeshift_segy_               = false;
  output_prenmo_time_segy_             = false;
}

ModelSettings::~ModelSettings(void)
{
}

//--------------------------------------------------------
void ModelSettings::CheckConsistency(std::string & errTxt)
//--------------------------------------------------------
{
  if (GetOutputPrenmoTimeSegy() && !GetNMOCorr()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: You cannot ask for pre NMO time output without specifying NMO correction. Output has been turned off.\n\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
    SetOutputPrenmoTimeSegy(false);
  }
  if (GetOutputTwtOffset() && !GetNMOCorr()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: You cannot ask for TWT offset output without specifying NMO correction. Output has been turned off.\n\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
    SetOutputTwtOffset(false);
  }
  if (GetUseHorizontalInterpolation() && GetUseVerticalInterpolation()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Warning, "WARNING: You cannot ask for both horizontal and vertical interpolation of layers.\n         Horizontal interpolation has been turned off.\n\n");
    TaskList::AddTask("Inconsistent XML model file specified. See beginning of log file.");
    SetUseHorizontalInterpolation(false);
  }

  //Eli Zachariassen
  //if (!(GetTimeOutput() || GetDepthOutput() || GetTimeshiftOutput())) {
  //  errTxt += "No seisimc output has been asked for. ";
  //}
}

void ModelSettings::SetDerivedVariables(void)
{

  //
  // Log file
  //
  if (GetPrefix() != "")
    log_file_name_ = GetPrefix() + "_" + log_file_name_;
  if (GetSuffix() != "")
    log_file_name_ += "_" + GetSuffix();

  log_file_name_ += ".txt";

  //
  // Setup angle or offset vectors
  //
  if (nmo_corr_) {
    size_t noffset = 1u;
    if (doffset_ > 0.0) {
      noffset  = size_t((offset_max_ - offset_0_) / doffset_) + 1u;
    }
    offset_vec_.resize(noffset);
    for (size_t i = 0; i < noffset; ++i) {
      offset_vec_[i] = offset_0_ + i*doffset_;
    }
  }
  else {
    size_t ntheta = 1u;

    if (dtheta_ > 0.0) {
      ntheta = static_cast<size_t>((theta_max_ - theta_0_) / dtheta_ + 1.01);
    }
    theta_vec_.resize(ntheta);
    for (size_t i = 0; i < ntheta; ++i) {
      theta_vec_[i] = theta_0_ + i*dtheta_;
    }
  }
}

bool ModelSettings::GetTimeOutput() const
{
  return (GetOutputTimeSegy()
          || GetOutputSeismicStackTimeSegy()
          || GetOutputSeismicTime()
          || GetOutputSeismicStackTimeStorm()
          || GetOutputPrenmoTimeSegy());
}

bool ModelSettings::GetDepthOutput() const
{
  return (GetOutputDepthSegy()
          || GetOutputSeismicStackDepthSegy()
          || GetOutputSeismicDepth()
          || GetOutputSeismicStackDepthStorm());
}

bool ModelSettings::GetTimeshiftOutput() const
{
  return (GetOutputTimeshiftSegy()
          || GetOutputSeismicStackTimeShiftSegy()
          || GetOutputSeismicTimeshift()
          || GetOutputSeismicStackTimeShiftStorm());
}

bool ModelSettings::GetStackOutput() const
{
  return (GetOutputSeismicStackTimeSegy()
          || GetOutputSeismicStackDepthSegy()
          || GetOutputSeismicStackTimeShiftSegy()
          || GetOutputSeismicStackTimeStorm()
          || GetOutputSeismicStackDepthStorm()
          || GetOutputSeismicStackTimeShiftStorm());
}

bool ModelSettings::GetSegyOutput() const
{
  return (GetOutputTimeSegy()
          || GetOutputSeismicStackTimeSegy()
          || GetOutputDepthSegy()
          || GetOutputSeismicStackDepthSegy()
          || GetOutputTimeshiftSegy()
          || GetOutputSeismicStackTimeShiftSegy()
          || GetOutputPrenmoTimeSegy());
}

bool ModelSettings::GetTimeStormOutput() const
{
  return (GetOutputSeismicStackTimeStorm()
          || GetOutputSeismicTime());
}

bool ModelSettings::GetDepthStormOutput() const
{
  return (GetOutputSeismicStackDepthStorm()
          || GetOutputSeismicDepth());
}

bool ModelSettings::GetTimeshiftStormOutput() const
{
  return (GetOutputSeismicStackTimeShiftStorm()
          || GetOutputSeismicTimeshift());
}

bool ModelSettings::GetStormOutput() const
{
  return (GetTimeStormOutput()
          || GetDepthStormOutput()
          || GetTimeshiftStormOutput());
}

void ModelSettings::PrintSettings(void)
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Max threads                                         : %10d\n", GetMaxThreads());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Traces in memory                                    : %10d\n", GetTracesInMemory());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Log level                                           : %10d\n", GetLogLevel());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Prefix                                              : %10s\n", GetPrefix().c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Suffix                                              : %10s\n", GetSuffix().c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic data type                                   : %10s\n"  , GetPSSeismic()                        ? "PS"  : "PP");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "NMO correction                                      : %10s\n"  , GetNMOCorr()                          ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Offset without stretch                              : %10s\n"  , GetOffsetWithoutStretch()             ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Use default underburden                             : %10s\n"  , GetDefaultUnderburden()               ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Resample parameters to Segy with interpolation      : %10s\n"  , GetResamplParamToSegyInterpol()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Regular grid\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Remove negative thicknesses                       : %10s\n"  , GetRemoveNegativeDeltaZ()             ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Minimum thickness for Eclipse grid cells          : %10.1f\n", GetZeroThicknessLimit());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use corner-point interpolation                    : %10s\n"  , GetUseCornerpointInterpol()           ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use fixed triangularization of Eclipe grid        : %10s\n"  , GetUseFixedTriangularization()        ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Interpolate at faults when corner-point           : %10s\n"  , GetCornerpointInterpolationAtFaults() ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use horizontal interpolation of layers            : %10s\n"  , GetUseHorizontalInterpolation()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use vertical interpolation of layers              : %10s\n"  , GetUseVerticalInterpolation()         ? "yes" : "no");
  if (GetUseBilinearInterpolation())
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use bilinear interpolation in regridding          :        yes\n");
  else
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use triangular interpolation                      :        yes\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use active pillars for center-point interpolation : %10s\n"  , GetUseActivePillars()                 ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Add white noise to seismic data                     : %10s\n"  , GetAddWhiteNoise()                    ? "yes" : "no");
  if (GetAddWhiteNoise()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use equal noise for offsets                       : %10s\n"  , GetUseEqualNoiseForOffsets()        ? "yes" : "no");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Standard deviation                                : %10.1f\n", GetStandardDeviation1());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Seed                                              : %10lu\n" , GetSeed1());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Add noise to reflection coefficients                : %10s\n"  , GetAddNoiseToReflCoef()               ? "yes" : "no");
  if (GetAddNoiseToReflCoef()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Standard deviation                                : %10.1f\n", GetStandardDeviation2());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Seed                                              : %10lu\n" , GetSeed2());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  ANGLES and OFFSETS
  //
  if (GetNMOCorr()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "NMO correction settings\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Sea floor depth                                   : %10.1f\n", GetZw());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp velocity in water                              : %10.1f\n", GetVw());
/*
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Offset span\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "    Minimum                                         : %10.1f\n", GetOffset0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "    Delta                                           : %10.1f\n", GetDOffset());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "    Maximum                                         : %10.1f\n", GetOffsetMax());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n");
*/
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Offsets                                           :  ");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," %.1f",offset_vec_[0]);
    for (size_t i = 1 ; i < GetOffsetVec().size() ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," -> %.1f",offset_vec_[i]);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"\n");
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "AVA angle span\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Minimum                                           : %10.1f\n", GetTheta0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Delta                                             : %10.1f\n", GetDTheta());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Maximum                                           : %10.1f\n", GetThetaMax());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Angles                                            :  ");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  %.1f", theta_vec_[0]);
    for (size_t i = 1 ; i < GetThetaVec().size() ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, " -> %.1f", theta_vec_[i]);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
  }

  //
  //  WAVELET
  //
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nWavelet\n");
  if (GetRicker()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Ricker with peak frequency                        : %10.1f\n", GetPeakFrequency());
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File                                              : %10s\n", GetWaveletFileName().c_str());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Format                                            : %10s\n", GetWaveletFileFormat().c_str());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Use zero time from header                         : %10s\n", GetUseZeroTimeFromHeader() ? "yes" : "no");
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Scale                                             : %10.1f\n", GetWaveletScale());
  if (GetWaveletLength() != -999.0)
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Length                                            : %10.1f\n", GetWaveletLength());
  if (GetWaveletLengthFactor() != -999.0)
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Length Factor                                     : %10.1f\n", GetWaveletLengthFactor());
  if (GetZWaveletTop() != 0.0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  z-extension above grid                            : %10.1f\n", GetZWaveletTop());
  }
  if (GetZWaveletBot() != 0.0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  z-extension belov grid                            : %10.1f\n", GetZWaveletBot());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  const std::vector<double> & vp  = GetConstVp();
  const std::vector<double> & vs  = GetConstVs();
  const std::vector<double> & rho = GetConstRho();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Parameters names\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp                                                : %10s\n", GetParameterNames()[0].c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vs                                                : %10s\n", GetParameterNames()[1].c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Rho                                               : %10s\n", GetParameterNames()[2].c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Default values for overburden, reservoir and underburden\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp                                                :%7.1f ->%7.1f ->%7.1f\n", vp[0] , vp[1] , vp[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vs                                                :%7.1f ->%7.1f ->%7.1f\n", vs[0] , vs[1] , vs[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Rho                                               :%7.1f ->%7.1f ->%7.1f\n", rho[0], rho[1], rho[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  size_t n = GetExtraParameterNames().size();
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Extra parameters                                    : %10s\n", n > 0 ? "yes" : "no");
  if (n > 0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Names and default values\n");
    for (size_t i = 0 ; i < n ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "    %20s                            : %10.1f\n",extra_parameter_names_[i].c_str(), extra_parameter_default_values_[i]);
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  AREA
  //
  if (GetAreaFromSegy() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Output area is taken from SegY file\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                                         : %10s\n", GetAreaFromSegy().c_str());
  }
  else if (GetAreaGiven()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Output area is specified in model file\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  x-start                                           : %10.1f\n", GetX0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  y-start                                           : %10.1f\n", GetY0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  x-length                                          : %10.1f\n", GetLx());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  y-length                                          : %10.1f\n", GetLy());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  angle                                             : %10.3f\n", NRLib::RadToDeg(GetAngle()));
  }
  else if (GetAreaFromSurface() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Output area is taken from surface\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                                         : %10s\n", GetAreaFromSurface().c_str());
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Output area is taken from Eclipse grid\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                                         : %10s\n", GetEclipseFileName().c_str());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Cell size\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  dx                                                : %10.1f\n", GetDx());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  dy                                                : %10.1f\n", GetDy());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  dt                                                : %10.1f\n", GetDt());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  dz                                                : %10.1f\n", GetDz());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  if (GetTopTimeSurfaceFile() != "")
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Top time surface                                    : %10s\n"  , GetTopTimeSurfaceFile().c_str());
  else
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Top time value                                      : %10.1f\n", GetTopTimeConstant());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  if (GetTimeWindowSpecified()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Specified time window\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Top:                                              : %10.1f\n", GetTopTimeWindow());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Base:                                             : %10.1f\n", GetBotTimeWindow());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  if (GetDepthWindowSpecified()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Specified depth window\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Top:                                              : %10.1f\n", GetTopDepthWindow());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Base:                                             : %10.1f\n", GetBotDepthWindow());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  OUTPUT
  //
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic data output\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputSeismicTime()                ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                                        : %10s\n", GetOutputSeismicTimeshift()           ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputSeismicDepth()               ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic stack time/depth in SEGY format\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Format name                                       : %10s\n", GetOutputSegyFileFormat().GetFormatName().c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputSeismicStackTimeSegy()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                                        : %10s\n", GetOutputSeismicStackTimeShiftSegy()  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputSeismicStackDepthSegy()      ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic stack time/depth in STORM format\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputSeismicStackTimeStorm()      ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                                        : %10s\n", GetOutputSeismicStackTimeShiftStorm() ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputSeismicStackDepthStorm()     ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Elastic parameters output\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputElasticParametersTimeSegy()  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputElasticParametersDepthSegy() ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Extra parameters output\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputExtraParametersTimeSegy()    ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputExtraParametersDepthSegy()   ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Other output\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  TWT                                               : %10s\n", GetOutputTwt()                        ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  TWT offset                                        : %10s\n", GetOutputTwtOffset()                  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vrms                                              : %10s\n", GetOutputVrms()                       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp                                                : %10s\n", GetOutputVp()                         ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Reflections                                       : %10s\n", GetOutputReflections()                ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time surfaces                                     : %10s\n", GetOutputTimeSurfaces()               ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth surfaces                                    : %10s\n", GetOutputDepthSurfaces()              ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                              : %10s\n", GetOutputTimeSegy()                   ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time pre NMO                                      : %10s\n", GetOutputPrenmoTimeSegy()             ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                                        : %10s\n", GetOutputTimeshiftSegy()              ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                             : %10s\n", GetOutputDepthSegy()                  ? "yes" : "no");
}
