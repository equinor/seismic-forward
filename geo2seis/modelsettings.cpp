// $Id: modelsettings.cpp 39 2014-02-06 08:53:13Z vigsnes $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •    Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •    Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include <stdio.h>
#include <string.h>
#include "nrlib/iotools/logkit.hpp"
#include "modelsettings.hpp"

ModelSettings::ModelSettings(void) {
    constvp_.resize(3);
    constvs_.resize(3);
    constrho_.resize(3);
    parameter_names_.resize(3);
    theta_0_ = 0.0;  // default values
    dtheta_ = 0.0;
    theta_max_ = 0.0;
    offset_0_ = 0.0;
    doffset_ = 0.0;
    offset_max_ = 0.0;
    top_time_surface_ = "";
    top_time_constant_ = 0.0;
    dx_ = 25.0;
    dy_ = 25.0;
    dz_ = 4.0;
    dt_ = 4.0;
    wavelet_scale_ = 1.0;
    output_vp_ = false;
    output_reflections_ = false;
    output_zvalues_ = false;
    output_seismic_time_ = false; //true when not nmo seismic
    output_seismic_depth_ = false;//true when not nmo seismic
    output_seismic_timeshift_ = false;
    output_time_surfaces_ = false;
    output_depth_surfaces_ = false;
    output_twt_ = false;
    output_vrms_ = false;
    output_twt_offset_ = false;
    nlayers_file_name_ = "";
    prefix_ = "";
    suffix_ = "";
    zero_thickness_limit_ = 0.1;
    output_time_segy_ = false;
    output_depth_segy_ = false;
    output_timeshift_segy_ = false;
    output_prenmo_time_segy_ = false;
    use_cornerpoint_interpol_ = false;
    area_from_surface_ = "";
    elastic_parameters_time_segy_ = false;
    elastic_parameters_depth_segy_ = false;
    extra_parameters_time_segy_ = false;
    extra_parameters_depth_segy_ = false;
    inline_start_ = 0;
    xline_start_ = 0;
    inline_direction_ = "y";
    inline_step_ = 1;
    xline_step_ = 1;
    seismic_stack_time_storm_ = false;
    seismic_stack_time_shift_storm_ = false;
    seismic_stack_depth_storm_ = false;
    seismic_stack_time_segy_ = false;
    seismic_stack_time_shift_segy_ = false;
    seismic_stack_depth_segy_ = false;
    white_noise_ = false;
    standard_deviation_ = 1.0;
    seed_ = static_cast<unsigned long>(std::time(0));
    il0_in_ = 189;
    xl0_in_ = 193;
    utmx_in_ = 181;
    utmy_in_ = 185;
    area_from_segy_ = "";
    utm_precision_ = -10;
    memory_limit_ = 8000000000.0;
    twt_file_name_ = "";
    ps_seismic_ = false;
    old_mod_ = false;
    nmo_corr_ = false;
    v_w_ = 0.0;
    z_w_ = 0.0;
    top_time_window_ = -9999;
    bot_time_window_ = -9999;
    top_depth_window_ = -9999;
    bot_depth_window_ = -9999;
    time_window_specified_ = false;
    depth_window_specified_ = false;
}

ModelSettings::~ModelSettings(void) {
}


