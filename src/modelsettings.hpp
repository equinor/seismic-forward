// $Id: modelsettings.hpp 41 2014-03-28 09:42:21Z vigsnes $

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

#ifndef MODELSETTINGS_HPP
#define MODELSETTINGS_HPP

#include "nrlib/segy/traceheader.hpp"

#include "nrlib/iotools/logkit.hpp"

#include "nrlib/math/constants.hpp"

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <ctime>

class ModelSettings
{
public:
  ModelSettings(void);

  ~ModelSettings(void);

  void                      CheckConsistency(std::string & errTxt);
  void                      SetDerivedVariables(void);

  void                      PrintSettings(void);

  size_t                    GetTracesInMemory(void)                   const { return traces_in_memory_               ;}
  size_t                    GetMaxThreads(void)                       const { return max_threads_                    ;}

  std::string               GetTwtFileName()                    const { return twt_file_name_                  ;}
  bool                      GetPSSeismic()                            const { return ps_seismic_                     ;}
  bool                      GetResamplParamToSegyInterpol()           const { return resampl_param_to_segy_with_interpol_ ;}
  double                    GetZExtrapolFactor()                      const { return z_extrapol_factor_              ;}
  double                    GetZeroThicknessLimit()                   const { return zero_thickness_limit_           ;}
  bool                      GetRemoveNegativeDeltaZ()                 const { return remove_negative_delta_z_        ;}
  bool                      GetUseCornerpointInterpol()               const { return use_cornerpoint_interpol_       ;}
  bool                      GetCornerpointInterpolationAtFaults()     const { return cornerpoint_interpol_at_faults_ ;}
  bool                      GetUseFixedTriangularization()            const { return use_fixed_triangularization_    ;}
  bool                      GetUseHorizontalInterpolation()           const { return use_horizontal_interpolation_   ;}
  bool                      GetUseVerticalInterpolation()             const { return use_vertical_interpolation_     ;}
  bool                      GetUseBilinearInterpolation()             const { return use_bilinear_interpolation_     ;}
  bool                      GetUseActivePillars()                     const { return use_active_pillars_             ;}

  int                       GetIL0Loc()                               const { return il0_loc_                        ;}
  int                       GetXL0Loc()                               const { return xl0_loc_                        ;}
  int                       GetUtmxLoc()                              const { return utmx_loc_                       ;}
  int                       GetUtmyLoc()                              const { return utmy_loc_                       ;}
  int                       GetScalcoLoc()                            const { return scalco_loc_                     ;}
  int                       GetStartTimeLoc()                         const { return start_time_loc_                 ;}
  short                     GetUtmPrecision()                         const { return utm_precision_                  ;}

  int                       GetSegyInlineStart()                      const { return inline_start_                   ;}
  int                       GetSegyXlineStart()                       const { return xline_start_                    ;}
  std::string               GetSegyInlineDirection()                  const { return inline_direction_               ;}
  int                       GetSegyInlineStep()                       const { return inline_step_                    ;}
  int                       GetSegyXlineStep()                        const { return xline_step_                     ;}

  bool                      GetAddNoiseToReflCoef()                   const { return add_noise_to_refl_coef_         ;}
  bool                      GetAddWhiteNoise()                        const { return add_white_noise_                ;}
  bool                      GetUseEqualNoiseForOffsets()              const { return equal_noise_for_offsets_        ;}
  double                    GetStandardDeviation1()                   const { return standard_deviation_1_           ;}
  double                    GetStandardDeviation2()                   const { return standard_deviation_2_           ;}
  unsigned long             GetSeed1()                                const { return seed_1_                         ;}
  unsigned long             GetSeed2()                                const { return seed_2_                         ;}

  std::string               GetEclipseFileName()                      const { return eclipse_file_name_              ;}

  std::vector<std::string>  GetParameterNames()                       const { return parameter_names_                ;}
  std::vector<double>       GetConstVp()                              const { return constvp_                        ;}
  std::vector<double>       GetConstVs()                              const { return constvs_                        ;}
  std::vector<double>       GetConstRho()                             const { return constrho_                       ;}
  bool                      GetDefaultUnderburden(void)               const { return default_underburden_            ;}

  std::vector<std::string>  GetExtraParameterNames()                  const { return extra_parameter_names_          ;}
  std::vector<double>       GetExtraParameterDefaultValues()          const { return extra_parameter_default_values_ ;}

  bool                      GetNMOCorr()                              const { return nmo_corr_                       ;}
  bool                      GetOffsetWithoutStretch()                 const { return offset_without_stretch_         ;}
  double                    GetVw()                                   const { return v_w_                            ;}
  double                    GetZw()                                   const { return z_w_                            ;}

  const std::vector<double> & GetThetaVec()                             const { return theta_vec_                      ;}
  const std::vector<double> & GetOffsetVec()                            const { return offset_vec_                     ;}

  double                    GetTheta0()                               const { return theta_0_                        ;}
  double                    GetDTheta()                               const { return dtheta_                         ;}
  double                    GetThetaMax()                             const { return theta_max_                      ;}

  bool                      GetRicker()                               const { return ricker_                         ;}
  double                    GetPeakFrequency()                        const { return peak_f_                         ;}
  std::string               GetWaveletFileFormat()                    const { return wavelet_file_format_            ;}
  std::string               GetWaveletFileName()                      const { return wavelet_file_name_              ;}
  double                    GetWaveletScale()                         const { return wavelet_scale_                  ;}
  double                    GetWaveletLength()                        const { return wavelet_length_                 ;}
  double                    GetWaveletLengthFactor()                  const { return wavelet_length_factor_          ;}
  bool                      GetUseZeroTimeFromHeader()                const { return use_zero_time_from_header_      ;}
  double                    GetZWaveletTop()                          const { return z_wavelet_top_                  ;}
  double                    GetZWaveletBot()                          const { return z_wavelet_bot_                  ;}

  //
  // Output parameters
  //
  const std::string       & GetLogFileName()                          const { return log_file_name_                  ;}
  int                       GetLogLevel()                             const { return log_level_                      ;}

  std::string               GetPrefix()                               const { return prefix_                         ;}
  std::string               GetSuffix()                               const { return suffix_                         ;}

  double                    GetTopTimeConstant()                      const { return top_time_constant_              ;}
  std::string               GetTopTimeSurfaceFile()                   const { return top_time_surface_               ;}

  bool                      GetTimeWindowSpecified()                  const { return time_window_specified_          ;}
  bool                      GetDepthWindowSpecified()                 const { return depth_window_specified_         ;}
  double                    GetTopTimeWindow()                        const { return top_time_window_                ;}
  double                    GetBotTimeWindow()                        const { return bot_time_window_                ;}
  double                    GetTopDepthWindow()                       const { return top_depth_window_               ;}
  double                    GetBotDepthWindow()                       const { return bot_depth_window_               ;}

  double                    GetPaddingFactorSeismicModelling()        const { return padding_factor_                 ;}

  bool                      GetAreaGiven()                            const { return area_given_                     ;}
  std::string               GetAreaFromSegy()                         const { return area_from_segy_                 ;}
  std::string               GetAreaFromSurface()                      const { return area_from_surface_              ;}
  double                    GetX0()                                   const { return x0_                             ;}
  double                    GetY0()                                   const { return y0_                             ;}
  double                    GetLx()                                   const { return lx_                             ;}
  double                    GetLy()                                   const { return ly_                             ;}
  double                    GetAngle()                                const { return angle_                          ;}

  double                    GetDx()                                   const { return dx_                             ;}
  double                    GetDy()                                   const { return dy_                             ;}
  double                    GetDz()                                   const { return dz_                             ;}
  double                    GetDt()                                   const { return dt_                             ;}

  bool                      GetOutputVp()                             const { return output_vp_                      ;}
  bool                      GetOutputReflections()                    const { return output_reflections_             ;}
  bool                      GetOutputZvalues()                        const { return output_zvalues_                 ;}
  bool                      GetOutputSeismicDepth()                   const { return output_seismic_depth_           ;}
  bool                      GetOutputSeismicTime()                    const { return output_seismic_time_            ;}
  bool                      GetOutputSeismicTimeshift()               const { return output_seismic_timeshift_       ;}
  bool                      GetOutputDepthSurfaces()                  const { return output_depth_surfaces_          ;}
  bool                      GetOutputTimeSurfaces()                   const { return output_time_surfaces_           ;}
  bool                      GetOutputWavelet()                        const { return output_wavelet_                 ;}
  bool                      GetOutputTwt()                            const { return output_twt_                     ;}
  bool                      GetOutputVrms()                           const { return output_vrms_                    ;}
  bool                      GetOutputTwtOffset()                      const { return output_twt_offset_              ;}
  bool                      GetOutputTimeSegy()                       const { return output_time_segy_               ;}
  bool                      GetOutputTimeshiftSegy()                  const { return output_timeshift_segy_          ;}
  bool                      GetOutputDepthSegy()                      const { return output_depth_segy_              ;}
  bool                      GetOutputPrenmoTimeSegy()                 const { return output_prenmo_time_segy_        ;}
  bool                      GetOutputElasticParametersTimeSegy()      const { return elastic_parameters_time_segy_   ;}
  bool                      GetOutputElasticParametersDepthSegy()     const { return elastic_parameters_depth_segy_  ;}
  bool                      GetOutputExtraParametersTimeSegy()        const { return extra_parameters_time_segy_     ;}
  bool                      GetOutputExtraParametersDepthSegy()       const { return extra_parameters_depth_segy_    ;}
  bool                      GetOutputSeismicStackTimeStorm()          const { return seismic_stack_time_storm_       ;}
  bool                      GetOutputSeismicStackTimeShiftStorm()     const { return seismic_stack_time_shift_storm_ ;}
  bool                      GetOutputSeismicStackDepthStorm()         const { return seismic_stack_depth_storm_      ;}
  bool                      GetOutputSeismicStackTimeSegy()           const { return seismic_stack_time_segy_        ;}
  bool                      GetOutputSeismicStackTimeShiftSegy()      const { return seismic_stack_time_shift_segy_  ;}
  bool                      GetOutputSeismicStackDepthSegy()          const { return seismic_stack_depth_segy_       ;}
  NRLib::TraceHeaderFormat  GetOutputSegyFileFormat()                 const { return output_segy_file_format_        ;}

  void SetLogLevel(int level)                          { log_level_                           = level    ;}

  void SetEclipseGrid(std::string filename)            { eclipse_file_name_                   = filename ;}
  void SetTwtFileName(std::string name)                { twt_file_name_                       = name     ;}

  void SetPrefix(std::string val)                      { prefix_                              = val      ;}
  void SetSuffix(std::string val)                      { suffix_                              = val      ;}
  void SetTracesInMemory(size_t value)                 { traces_in_memory_                    = value    ;}
  void SetMaxThreads(size_t value)                     { max_threads_                         = value    ;}

  void SetZeroThicknessLimit(double val)               { zero_thickness_limit_                = val     ;}

  void SetVw(double value)                             { v_w_                                 = value   ;}
  void SetZw(double value)                             { z_w_                                 = value   ;}
  void SetZExtrapolFactor(double fact)                 { z_extrapol_factor_                   = fact    ;}

  void SetOffset0(double value)                        { offset_0_                            = value   ;}
  void SetDOffset(double value)                        { doffset_                             = value   ;}
  void SetOffsetMax(double value)                      { offset_max_                          = value   ;}
  void SetOffsetWithoutStretch(bool value)             { offset_without_stretch_              = value   ;}

  void SetAddNoiseToReflCoef(void)                     { add_noise_to_refl_coef_              = true    ;}
  void SetAddWhiteNoise(void)                          { add_white_noise_                     = true    ;}
  void SetUseEqualNoiseForOffsets(bool equal)          { equal_noise_for_offsets_             = equal   ;}
  void SetStandardDeviation1(double value)             { standard_deviation_1_                = value   ;}
  void SetStandardDeviation2(double value)             { standard_deviation_2_                = value   ;}
  void SetSeed1(unsigned long value)                   { seed_1_                              = value   ;}
  void SetSeed2(unsigned long value)                   { seed_2_                              = value   ;}

  void SetUseCornerpointInterpol(bool value)           { use_cornerpoint_interpol_            = value   ;}
  void SetCornerpointInterpolationAtFaults(bool value) { cornerpoint_interpol_at_faults_      = value   ;}
  void SetUseFixedTriangularization(bool value)        { use_fixed_triangularization_         = value   ;}
  void SetUseHorizontalInterpolation(bool value)       { use_horizontal_interpolation_        = value   ;}
  void SetUseVerticalInterpolation(bool value)         { use_vertical_interpolation_          = value   ;}
  void SetUseBilinearInterpolation(bool value)         { use_bilinear_interpolation_          = value   ;}
  void SetUseActivePillars(bool value)                 { use_active_pillars_                  = value   ;}
  void SetRemoveNegativeDeltaZ(bool value)             { remove_negative_delta_z_             = value   ;}
  void SetPSSeismic(bool ps)                           { ps_seismic_                          = ps      ;}
  void SetDefaultUnderburden(bool value)               { default_underburden_                 = value   ;}
  void SetResamplParamToSegyInterpol(bool value)       { resampl_param_to_segy_with_interpol_ = value   ;}
  void SetNMOCorr(bool nmo)                            { nmo_corr_                            = nmo     ;}

  void SetVpTop(double vptop)                          { constvp_[0]               = vptop  ;}
  void SetVpMid(double vpmid)                          { constvp_[1]               = vpmid  ;}
  void SetVpBot(double vpbot)                          { constvp_[2]               = vpbot  ;}

  void SetVsTop(double vstop)                          { constvs_[0]               = vstop  ;}
  void SetVsMid(double vsmid)                          { constvs_[1]               = vsmid  ;}
  void SetVsBot(double vsbot)                          { constvs_[2]               = vsbot  ;}

  void SetRhoTop(double rhotop)                        { constrho_[0]              = rhotop ;}
  void SetRhoMid(double rhomid)                        { constrho_[1]              = rhomid ;}
  void SetRhoBot(double rhobot)                        { constrho_[2]              = rhobot ;}

  void SetVpName(std::string & name)                   { parameter_names_[0]       = name ;}
  void SetVsName(std::string & name)                   { parameter_names_[1]       = name ;}
  void SetRhoName(std::string & name)                  { parameter_names_[2]       = name ;}

  void SetTheta0(double theta)                         { theta_0_                  = NRLib::Degree * theta ;}    // in radians
  void SetDTheta(double theta)                         { dtheta_                   = NRLib::Degree * theta ;}
  void SetThetaMax(double theta)                       { theta_max_                = NRLib::Degree * theta ;}
  void SetAreaAngle(double angle)                      { angle_                    = NRLib::Degree * angle ;}    // in radians
  void SetAreaGiven(bool areagiven)                    { area_given_               = areagiven             ;}
  void SetAreaFromSurface(std::string val)             { area_from_surface_        = val                   ;}
  void SetAreaFromSegy(std::string val)                { area_from_segy_           = val                   ;}
  void SetUtmPrecision(short value)                    { utm_precision_            = value                 ;}

  void SetTopTimeSurface(std::string toptime)          { top_time_surface_         = toptime ;}
  void SetTopTimeConstant(double value)                { top_time_constant_        = value   ;}
  void SetTopTimeWindow(double value)                  { top_time_window_          = value   ;}
  void SetBotTimeWindow(double value)                  { bot_time_window_          = value   ;}
  void SetTopDepthWindow(double value)                 { top_depth_window_         = value   ;}
  void SetBotDepthWindow(double value)                 { bot_depth_window_         = value   ;}
  void SetTimeWindowSpecified(bool value)              { time_window_specified_    = value   ;}
  void SetDepthWindowSpecified(bool value)             { depth_window_specified_   = value   ;}
  void SetPaddingFactorSeismicModelling(double value)  { padding_factor_           = value   ;}

  void SetX0(double x0)                                { x0_                       = x0 ;}
  void SetY0(double y0)                                { y0_                       = y0 ;}
  void SetLx(double lx)                                { lx_                       = lx ;}
  void SetLy(double ly)                                { ly_                       = ly ;}

  void SetDx(double dx)                                { dx_                       = dx ;}
  void SetDy(double dy)                                { dy_                       = dy ;}
  void SetDz(double dz)                                { dz_                       = dz ;}
  void SetDt(double dt)                                { dt_                       = dt ;}

  void SetIL0Loc(int value)                            { il0_loc_                  = value  ;}
  void SetXL0Loc(int value)                            { xl0_loc_                  = value  ;}
  void SetUtmxLoc(int value)                           { utmx_loc_                 = value  ;}
  void SetUtmyLoc(int value)                           { utmy_loc_                 = value  ;}
  void SetScalcoLoc(int value)                         { scalco_loc_               = value  ;}
  void SetStartTimeLoc(int value)                      { start_time_loc_           = value  ;}

  void SetSegyInlineStart(int value)                   { inline_start_             = value  ;}
  void SetSegyXlineStart(int value)                    { xline_start_              = value  ;}
  void SetSegyInlineDirection(std::string value)       { inline_direction_         = value  ;}
  void SetSegyInlineStep(int value)                    { inline_step_              = value  ;}
  void SetSegyXlineStep(int value)                     { xline_step_               = value  ;}

  void SetRicker(bool ricker)                          { ricker_                   = ricker ;}
  void SetPeakF(double peakf)                          { peak_f_                   = peakf  ;}

  void SetWaveletFileFormat(std::string value)         { wavelet_file_format_      = value  ;}
  void SetWaveletFileName(std::string value)           { wavelet_file_name_        = value  ;}
  void SetWaveletScale(double scale)                   { wavelet_scale_            = scale  ;}
  void SetWaveletLength(double length)                 { wavelet_length_           = length ;}
  void SetWaveletLengthFactor(double value)            { wavelet_length_factor_    = value  ;}
  void SetUseZeroTimeFromHeader(bool value)            { use_zero_time_from_header_= value  ;}
  void SetZWaveletTop(double wave)                     { z_wavelet_top_            = wave   ;}
  void SetZWaveletBot(double wave)                     { z_wavelet_bot_            = wave   ;}

  void SetOutputVp(bool value)                         { output_vp_                = value  ;}
  void SetOutputReflections(bool value)                { output_reflections_       = value  ;}
  void SetOutputZvalues(bool value)                    { output_zvalues_           = value  ;}
  void SetOutputSeismicDepth(bool value)               { output_seismic_depth_     = value  ;}
  void SetOutputSeismicTime(bool value)                { output_seismic_time_      = value  ;}
  void SetOutputSeismicTimeshift(bool value)           { output_seismic_timeshift_ = value  ;}
  void SetOutputDepthSurfaces(bool value)              { output_depth_surfaces_    = value  ;}
  void SetOutputTimeSurfaces(bool value)               { output_time_surfaces_     = value  ;}
  void SetOutputTwt(bool value)                        { output_twt_               = value  ;}
  void SetOutputVrms(bool value)                       { output_vrms_              = value  ;}
  void SetOutputTwtOffset(bool value)                  { output_twt_offset_        = value  ;}
  void SetOutputTimeSegy(bool value)                   { output_time_segy_         = value  ;}
  void SetOutputTimeshiftSegy(bool value)              { output_timeshift_segy_    = value  ;}
  void SetOutputDepthSegy(bool value)                  { output_depth_segy_        = value  ;}
  void SetOutputPrenmoTimeSegy(bool value)             { output_prenmo_time_segy_  = value  ;}
  void SetOutputWavelet(bool value)                    { output_wavelet_           = value  ;}

  void SetOutputSeismicStackTimeStorm(bool value)      { seismic_stack_time_storm_       = value ;}
  void SetOutputSeismicStackTimeShiftStorm(bool value) { seismic_stack_time_shift_storm_ = value ;}
  void SetOutputSeismicStackDepthStorm(bool value)     { seismic_stack_depth_storm_      = value ;}
  void SetOutputSeismicStackTimeSegy(bool value)       { seismic_stack_time_segy_        = value ;}
  void SetOutputSeismicStackTimeShiftSegy(bool value)  { seismic_stack_time_shift_segy_  = value ;}
  void SetOutputSeismicStackDepthSegy(bool value)      { seismic_stack_depth_segy_       = value ;}

  void SetOutputSegyFileFormat(NRLib::TraceHeaderFormat thf) { output_segy_file_format_  = thf   ;}

  void SetOutputElasticParametersTimeSegy(bool value)  { elastic_parameters_time_segy_   = value ;}
  void SetOutputElasticParametersDepthSegy(bool value) { elastic_parameters_depth_segy_  = value ;}
  void SetOutputExtraParametersTimeSegy(bool value)    { extra_parameters_time_segy_     = value ;}
  void SetOutputExtraParametersDepthSegy(bool value)   { extra_parameters_depth_segy_    = value ;}

  void AddExtraParameterName(std::string name)         { extra_parameter_names_.push_back(name);}
  void AddExtraParameterDefaultValue(double value)     { extra_parameter_default_values_.push_back(value);}

  bool GetTimeOutput() const;
  bool GetDepthOutput() const;
  bool GetTimeshiftOutput() const;
  bool GetStackOutput();
  bool GetSegyOutput();
  bool GetTimeStormOutput();
  bool GetDepthStormOutput();
  bool GetTimeshiftStormOutput();
  bool GetStormOutput();

private:

  std::string               prefix_;
  std::string               suffix_;
  std::string               log_file_name_;
  int                       log_level_;

  short                     utm_precision_;
  size_t                    traces_in_memory_;
  size_t                    max_threads_;
  double                    zero_thickness_limit_;

  bool                      add_noise_to_refl_coef_;
  bool                      add_white_noise_;
  bool                      equal_noise_for_offsets_;
  double                    standard_deviation_1_;
  double                    standard_deviation_2_;
  unsigned long             seed_1_;
  unsigned long             seed_2_;

  bool                      default_underburden_;
  bool                      ps_seismic_;
  bool                      nmo_corr_;
  bool                      resampl_param_to_segy_with_interpol_;

  std::vector<double>       constvp_;
  std::vector<double>       constvs_;
  std::vector<double>       constrho_;

  std::vector<std::string>  parameter_names_;
  double                    theta_0_;        // seismic angle
  double                    dtheta_;         // seismic angle
  double                    theta_max_;      // seismic angle
  std::vector<double>       theta_vec_;      // seismic angle
  double                    offset_0_;       // seismic offset
  double                    doffset_;        // seismic offset
  double                    offset_max_;     // seismic offset

  std::vector<double>       offset_vec_;

  std::string               eclipse_file_name_;
  double                    x0_;             // area parameters
  double                    y0_;             // area parameters
  double                    lx_;             // area parameters
  double                    ly_;             // area parameters
  double                    angle_;          // area parameters

  int                       inline_start_;
  int                       xline_start_;
  std::string               inline_direction_;
  int                       inline_step_;
  int                       xline_step_;

  double                    dx_;
  double                    dy_;
  double                    dz_;
  double                    dt_; // cell size

  int                       il0_loc_;
  int                       xl0_loc_;
  int                       utmx_loc_;
  int                       utmy_loc_;
  int                       scalco_loc_;
  int                       start_time_loc_;

  std::string               wavelet_file_format_;
  std::string               wavelet_file_name_;
  double                    peak_f_;
  bool                      ricker_;
  double                    wavelet_scale_;
  double                    wavelet_length_;
  double                    wavelet_length_factor_;   // Multiply (1/2 wavelet) length by this number.
  bool                      use_zero_time_from_header_;
  double                    z_wavelet_top_;
  double                    z_wavelet_bot_;
  double                    z_extrapol_factor_;
  bool                      offset_without_stretch_;

  std::string               twt_file_name_;
  std::string               area_from_segy_;
  std::string               area_from_surface_;
  bool                      area_given_;
  std::string               top_time_surface_;
  double                    top_time_constant_;

  NRLib::TraceHeaderFormat  output_segy_file_format_;
  bool                      output_time_segy_;
  bool                      output_depth_segy_;
  bool                      output_timeshift_segy_;
  bool                      output_prenmo_time_segy_;
  bool                      output_vp_;
  bool                      output_reflections_;
  bool                      output_zvalues_;
  bool                      output_seismic_time_;
  bool                      output_seismic_depth_;
  bool                      output_seismic_timeshift_;
  bool                      output_time_surfaces_;
  bool                      output_depth_surfaces_;
  bool                      output_wavelet_;
  bool                      output_twt_;
  bool                      output_vrms_;
  bool                      output_twt_offset_;

  bool                      use_cornerpoint_interpol_;
  bool                      cornerpoint_interpol_at_faults_;
  bool                      use_fixed_triangularization_;
  bool                      use_horizontal_interpolation_;
  bool                      use_vertical_interpolation_;
  bool                      use_bilinear_interpolation_;
  bool                      use_active_pillars_;
  bool                      remove_negative_delta_z_;

  bool                      elastic_parameters_time_segy_;
  bool                      elastic_parameters_depth_segy_;
  bool                      extra_parameters_time_segy_;
  bool                      extra_parameters_depth_segy_;

  std::vector<std::string>  extra_parameter_names_;
  std::vector<double>       extra_parameter_default_values_;

  bool                      seismic_stack_time_storm_;
  bool                      seismic_stack_time_shift_storm_;
  bool                      seismic_stack_depth_storm_;
  bool                      seismic_stack_time_segy_;
  bool                      seismic_stack_time_shift_segy_;
  bool                      seismic_stack_depth_segy_;

  double                    v_w_;
  double                    z_w_;

  double                    top_time_window_;
  double                    bot_time_window_;
  double                    top_depth_window_;
  double                    bot_depth_window_;

  double                    padding_factor_;    // Padding factor for seismic modelling

  bool                      time_window_specified_;
  bool                      depth_window_specified_;
};

#endif
