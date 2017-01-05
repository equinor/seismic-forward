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

ModelSettings::ModelSettings(void)
{
  constvp_.resize(3);
  constvs_.resize(3);
  constrho_.resize(3);
  parameter_names_.resize(3);

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
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "Seed                                      : %10d\n", GetSeed());

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nNMO correction                            : %10s\n", GetNMOCorr() ? "yes" : "no");

  if (GetNMOCorr()) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nOffset span\n");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Minimum                                 : %10.1f\n", GetOffset0());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Delta                                   : %10.1f\n", GetDOffset());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "  Maximum                                 : %10.1f\n", GetOffsetMax());
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nOffsets                                   :  ");
    size_t n = GetOffsetVec().size();
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"{%.1f,",offset_vec_[0]);
    for (size_t i = 1 ; i < n - 1 ; i++)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," %.1f,",offset_vec_[i]);
    if (n > 1)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low," %.1f",offset_vec_[n - 1]);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low,"}\n");
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nAVA angle span\n");
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
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "}\n");
  }
}
