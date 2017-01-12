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

#include "nrlib/iotools/logkit.hpp"

#include "modelsettings.hpp"
#include <iostream>
#include <cstdlib>

ModelSettings::ModelSettings(void)
{
  constvp_.resize(3);
  constvs_.resize(3);
  constrho_.resize(3);
  parameter_names_.resize(3);

  log_file_name_                       = "Logfile";
  log_level_                           = NRLib::LogKit::L_Low;

  prefix_                              = "";
  suffix_                              = "";

  seed_                                = static_cast<unsigned long>(std::time(0));
  standard_deviation_                  = 1.0;
  zero_thickness_limit_                = 0.1;
  twt_file_name_                       = "";
  traces_in_memory_                    = 100000;
  max_threads_                         = 100;
  white_noise_                         = false;
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

  il0_in_                              = 189;
  xl0_in_                              = 193;
  utmx_in_                             = 181;
  utmy_in_                             = 185;
  utm_precision_                       = -10;

  inline_start_                        = 0;
  xline_start_                         = 0;
  inline_direction_                    = "y";
  inline_step_                         = 1;
  xline_step_                          = 1;

  wavelet_scale_                       = 1.0;
  z_wavelet_top_                       = 0.0;
  z_wavelet_bot_                       = 0.0;
  z_extrapol_factor_                   = 50.0;
  offset_without_stretch_              = false;

  top_time_window_                     = -9999;
  bot_time_window_                     = -9999;
  top_depth_window_                    = -9999;
  bot_depth_window_                    = -9999;

  time_window_specified_               = false;
  depth_window_specified_              = false;

  use_cornerpoint_interpol_            = false;
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

void ModelSettings::PrintSettings(void)
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Max threads                               : %10d\n", GetMaxThreads());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Traces in memory                          : %10d\n", GetTracesInMemory());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Log level                                 : %10d\n", GetLogLevel());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Prefix                                    : %10s\n", GetPrefix().c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Suffix                                    : %10s\n", GetSuffix().c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seed                                      : %10d\n", GetSeed());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic data type                         : %10s\n", GetPSSeismic()                  ? "PS"  : "PP");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "NMO correction                            : %10s\n", GetNMOCorr()                    ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Offset without stretch                    : %10s\n", GetOffsetWithoutStretch()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Remove negative thicknesses               : %10s\n", GetRemoveNegativeDeltaZ()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Use corner-point interpolation            : %10s\n", GetUseCornerpointInterpol()     ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Use default underburden                   : %10s\n", GetDefaultUnderburden()         ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Resample parameters to Segy with interpol.: %10s\n", GetResamplParamToSegyInterpol() ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Add white noise                           : %10s\n", GetWhiteNoise()                 ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Standard deviation                      : %10.1f\n", GetStandardDeviation());

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  WAVELET
  //
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Wavelet\n");
  if (GetRicker()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Ricker with peak frequency              : %10.1f\n", GetPeakFrequency());
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File                                    : %10.1f\n", GetWaveletFileName().c_str());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Format                                  : %10.1f\n", GetWaveletFileFormat().c_str());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Scale                                   : %10.1f\n", GetWaveletScale());
  if (GetZWaveletTop() != 0.0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  z-extension above grid                  : %10.1f\n", GetZWaveletTop());
  }
  if (GetZWaveletBot() != 0.0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  z-extension belov grid                  : %10.1f\n", GetZWaveletBot());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  ANGLES and OFFSETS
  //
  if (GetNMOCorr()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Offset span\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Minimum                                 : %10.1f\n", GetOffset0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Delta                                   : %10.1f\n", GetDOffset());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Maximum                                 : %10.1f\n", GetOffsetMax());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Offsets                                 :  ");
    size_t n = GetOffsetVec().size();
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"{%.1f,",offset_vec_[0]);
    for (size_t i = 1 ; i < n - 1 ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," %.1f,",offset_vec_[i]);
    if (n > 1)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," %.1f",offset_vec_[n - 1]);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"}\n\n");
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "AVA angle span\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Minimum                                 : %10.1f\n", GetTheta0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Delta                                   : %10.1f\n", GetDTheta());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Maximum                                 : %10.1f\n", GetThetaMax());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Angles                                    :  ");
    size_t n = GetThetaVec().size();
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, " {%.1f,", theta_vec_[0]);
    for (size_t i = 1 ; i < n - 1 ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, " %.1f,", theta_vec_[i]);
    if (n > 1)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, " %.1f", theta_vec_[n - 1]);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "}\n\n");
  }

  //
  //  AREA
  //
  if (GetAreaFromSegy() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Area is taken from SegY file\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                               : %10s\n\n", GetAreaFromSegy().c_str());
  }
  else if (GetAreaGiven()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Area is specified in model file:\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  x-start                                 : %10.1f\n", GetX0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  y-start                                 : %10.1f\n", GetY0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  x-length                                : %10.1f\n", GetLx());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  y-length                                : %10.1f\n", GetLy());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  angle                                   : %10.1f\n", GetAngle());
  }
  else if (GetAreaFromSurface() != "") {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Area is taken from surface\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                               : %10s\n", GetAreaFromSurface().c_str());
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Area is taken from Eclipse grid\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  File name                               : %10s\n", GetEclipseFileName().c_str());
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  const std::vector<double> & vp  = GetConstVp();
  const std::vector<double> & vs  = GetConstVs();
  const std::vector<double> & rho = GetConstRho();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Default values for overburden, reservoir and underburden: \n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp                                      : %7.1f ->%7.1f ->%7.1f\n", vp[0] , vp[1] , vp[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vs                                      : %7.1f ->%7.1f ->%7.1f\n", vs[0] , vs[1] , vs[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Rho                                     : %7.1f ->%7.1f ->%7.1f\n", rho[0], rho[1], rho[2]);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  size_t n = GetExtraParameterNames().size();
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Extra parameters                          : %10s\n", n > 0 ? "yes" : "no");
  if (n > 0) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Names and default values:\n");
    for (size_t i = 0 ; i < n ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "    %20s                  : %10.1f\n",extra_parameter_names_[i].c_str(), extra_parameter_default_values_[i]);
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");

  //
  //  OUTPUT
  //
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic data:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputSeismicTime()                ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                              : %10s\n", GetOutputSeismicTimeshift()           ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputSeismicDepth()               ? "yes" : "no");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic stack time/depth in SEGY format:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputSeismicStackTimeSegy()       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                              : %10s\n", GetOutputSeismicStackTimeShiftSegy()  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputSeismicStackDepthSegy()      ? "yes" : "no");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seismic stack time/depth in STORM format:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputSeismicStackTimeStorm()      ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                              : %10s\n", GetOutputSeismicStackTimeShiftStorm() ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputSeismicStackDepthStorm()     ? "yes" : "no");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Elastic parameters output:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputElasticParametersTimeSegy()  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputElasticParametersDepthSegy() ? "yes" : "no");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Extra parameters output:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputExtraParametersTimeSegy()    ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputExtraParametersDepthSegy()   ? "yes" : "no");

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Other output:\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  TWT                                     : %10s\n", GetOutputTwt()                        ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  TWT offset                              : %10s\n", GetOutputTwtOffset()                  ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vrms                                    : %10s\n", GetOutputVrms()                       ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Vp                                      : %10s\n", GetOutputVp()                         ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Reflections                             : %10s\n", GetOutputReflections()                ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time surfaces                           : %10s\n", GetOutputTimeSurfaces()               ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth surfaces                          : %10s\n", GetOutputDepthSurfaces()              ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time                                    : %10s\n", GetOutputTimeSegy()                   ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time pre NMO                            : %10s\n", GetOutputPrenmoTimeSegy()             ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Time shift                              : %10s\n", GetOutputTimeshiftSegy()              ? "yes" : "no");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Depth                                   : %10s\n", GetOutputDepthSegy()                  ? "yes" : "no");


/*

  int                       GetSegyInlineStart()                      { return inline_start_                   ;}
  int                       GetSegyXlineStart()                       { return xline_start_                    ;}
  std::string               GetSegyInlineDirection()                  { return inline_direction_               ;}
  int                       GetSegyInlineStep()                       { return inline_step_                    ;}
  int                       GetSegyXlineStep()                        { return xline_step_                     ;}

  double                    GetTopTimeConstant()                      { return top_time_constant_              ;}
  std::string               GetTopTimeSurfaceFile()                   { return top_time_surface_               ;}
  double                    GetZExtrapolFactor()                      { return z_extrapol_factor_              ;}
  double                    GetZeroThicknessLimit()                   { return zero_thickness_limit_           ;}

  double                    GetDx()                                   { return dx_                             ;}
  double                    GetDy()                                   { return dy_                             ;}
  double                    GetDz()                                   { return dz_                             ;}
  double                    GetDt()                                   { return dt_                             ;}

  int                       GetIL0In()                                { return il0_in_                         ;}
  int                       GetXL0In()                                { return xl0_in_                         ;}
  int                       GetUtmxIn()                               { return utmx_in_                        ;}
  int                       GetUtmyIn()                               { return utmy_in_                        ;}
  short                     GetUtmPrecision()                         { return utm_precision_                  ;}

  double                    GetTopTimeWindow()                        { return top_time_window_                ;}
  double                    GetBotTimeWindow()                        { return bot_time_window_                ;}
  double                    GetTopDepthWindow()                       { return top_depth_window_               ;}
  double                    GetBotDepthWindow()                       { return bot_depth_window_               ;}

  bool                      GetTimeWindowSpecified()                  { return time_window_specified_          ;}
  bool                      GetDepthWindowSpecified()                 { return depth_window_specified_         ;}

  std::string               GetTwtFileName()                    const { return twt_file_name_                  ;}
  double                    GetVw()                                   { return v_w_                            ;}
  double                    GetZw()                                   { return z_w_                            ;}

*/

}
