// $Id: modelsettings.hpp 41 2014-03-28 09:42:21Z vigsnes $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// �    Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// �    Redistributions in binary form must reproduce the above copyright notice, this list of
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

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/math/constants.hpp"

class ModelSettings {
  public:
    ModelSettings(void);

    ~ModelSettings();

    void SetEclipseGrid(std::string filename) {
        eclipse_file_name_ = filename;
    }

    void SetVpTop(double vptop) {
        constvp_[0] = vptop;
    }

    void SetVpMid(double vpmid) {
        constvp_[1] = vpmid;
    }

    void SetVpBot(double vpbot) {
        constvp_[2] = vpbot;
    }

    void SetVsTop(double vstop) {
        constvs_[0] = vstop;
    }

    void SetVsMid(double vsmid) {
        constvs_[1] = vsmid;
    }

    void SetVsBot(double vsbot) {
        constvs_[2] = vsbot;
    }

    void SetRhoTop(double rhotop) {
        constrho_[0] = rhotop;
    }

    void SetRhoMid(double rhomid) {
        constrho_[1] = rhomid;
    }

    void SetRhoBot(double rhobot) {
        constrho_[2] = rhobot;
    }

    void SetVpName(std::string &name) {
        parameter_names_[0] = name;
    }

    void SetVsName(std::string &name) {
        parameter_names_[1] = name;
    }

    void SetRhoName(std::string &name) {
        parameter_names_[2] = name;
    }

    void SetTheta0(double theta) {
        theta_0_ = NRLib::Degree * theta;    // in radians
    }

    void SetDTheta(double theta) {
        dtheta_ = NRLib::Degree * theta;
    }

    void SetThetaMax(double theta) {
        theta_max_ = NRLib::Degree * theta;
    }

    void SetRicker(bool ricker) {
        ricker_ = ricker;
    }

    void SetPeakF(double peakf) {
        peak_f_ = peakf;
    }

    void SetX0(double x0) {
        x0_ = x0;
    }

    void SetY0(double y0) {
        y0_ = y0;
    }

    void SetLx(double lx) {
        lx_ = lx;
    }

    void SetLy(double ly) {
        ly_ = ly;
    }

    void SetAreaAngle(double angle) {
        angle_ = NRLib::Degree * angle;    // in radians
    }

    //void SetTopSurfaceFile(std::string top) { top_surface_file_ = top; }
    //void SetBotSurfaceFile(std::string bot) { bot_surface_file_ = bot; }
    void SetTopTimeSurface(std::string toptime) {
        top_time_surface_ = toptime;
    }

    void SetTopTimeConstant(double value) {
        top_time_constant_ = value;
    }

    //void SetConstantTop(double top) { constant_top_ = top;}
    //void SetConstantBot(double bot) { constant_bot_ = bot;}
    void SetDx(double dx) {
        dx_ = dx;
    }

    void SetDy(double dy) {
        dy_ = dy;
    }

    void SetDz(double dz) {
        dz_ = dz;
    }

    void SetDt(double dt) {
        dt_ = dt;
    }

    void SetAreaGiven(bool areagiven) {
        area_given_ = areagiven;
    }

    void SetWaveletScale(double scale) {
        wavelet_scale_ = scale;
    }

    void SetOutputVp(bool value) {
        output_vp_ = value;
    }

    void SetOutputReflections(bool value) {
        output_reflections_ = value;
    }

    void SetOutputZvalues(bool value) {
        output_zvalues_ = value;
    }

    void SetOutputSeismicDepth(bool value) {
        output_seismic_depth_ = value;
    }

    void SetOutputSeismicTime(bool value) {
        output_seismic_time_ = value;
    }

    void SetOutputSeismicTimeshift(bool value) {
        output_seismic_timeshift_ = value;
    }

    void SetOutputDepthSurfaces(bool value) {
        output_depth_surfaces_ = value;
    }

    void SetOutputTimeSurfaces(bool value) {
        output_time_surfaces_ = value;
    }

    void SetOutputTwt(bool value) {
        output_twt_ = value;
    }

    void SetNLayersFile(std::string val) {
        nlayers_file_name_ = val;
    }

    void SetPrefix(std::string val) {
        prefix_ = val;
    }

    void SetSuffix(std::string val) {
        suffix_ = val;
    }

    void SetZeroThicknessLimit(double val) {
        zero_thickness_limit_ = val;
    }

    void SetOutputTimeSegy(bool value) {
        output_time_segy_ = value;
    }

    void SetOutputTimeshiftSegy(bool value) {
        output_timeshift_segy_ = value;
    }

    void SetOutputDepthSegy(bool value) {
        output_depth_segy_ = value;
    }

    void SetUseCornerpointInterpol(bool value) {
        use_cornerpoint_interpol_ = value;
    }

    void SetAreaFromSurface(std::string val) {
        area_from_surface_ = val;
    }

    ///--------------------------------------------------------------------------------------------------
    void SetOutputElasticParametersTimeSegy(bool value) {
        elastic_parameters_time_segy_ = value;
    }

    void SetOutputElasticParametersDepthSegy(bool value) {
        elastic_parameters_depth_segy_ = value;
    }

    void SetOutputExtraParametersTimeSegy(bool value) {
        extra_parameters_time_segy_ = value;
    }

    void SetOutputExtraParametersDepthSegy(bool value) {
        extra_parameters_depth_segy_ = value;
    }

    void SetSegyInlineStart(int value) {
        inline_start_ = value;
    }

    void SetSegyXlineStart(int value) {
        xline_start_ = value;
    }

    void SetSegyInlineDirection(std::string value) {
        inline_direction_ = value;
    }

    void SetSegyInlineStep(int value) {
        inline_step_ = value;
    }

    void SetSegyXlineStep(int value) {
        xline_step_ = value;
    }

    void SetWaveletFileFormat(std::string value) {
        wavelet_file_format_ = value;
    }

    void SetWaveletFileName(std::string value) {
        wavelet_file_name_ = value;
    }

    void AddExtraParameterName(std::string name) {
        extra_parameter_names_.push_back(name);
    }

    void AddExtraParameterDefaultValue(double value) {
        extra_parameter_default_values_.push_back(value);
    }

    void SetOutputSeismicStackTimeStorm(bool value) {
        seismic_stack_time_storm_ = value;
    }

    void SetOutputSeismicStackTimeShiftStorm(bool value) {
        seismic_stack_time_shift_storm_ = value;
    }

    void SetOutputSeismicStackDepthStorm(bool value) {
        seismic_stack_depth_storm_ = value;
    }

    void SetOutputSeismicStackTimeSegy(bool value) {
        seismic_stack_time_segy_ = value;
    }

    void SetOutputSeismicStackTimeShiftSegy(bool value) {
        seismic_stack_time_shift_segy_ = value;
    }

    void SetOutputSeismicStackDepthSegy(bool value) {
        seismic_stack_depth_segy_ = value;
    }

    void SetWhiteNoise() {
        white_noise_ = true;
    }

    void SetStandardDeviation(double value) {
        standard_deviation_ = value;
    }

    void SetSeed(double value) {
        seed_ = static_cast<unsigned long>(value);
    }

    void SetIL0In(int value) {
        il0_in_ = value;
    }

    void SetXL0In(int value) {
        xl0_in_ = value;
    }

    void SetUtmxIn(int value) {
        utmx_in_ = value;
    }

    void SetUtmyIn(int value) {
        utmy_in_ = value;
    }

    void SetAreaFromSegy(std::string val) {
        area_from_segy_ = val;
    }

    void SetUtmPrecision(short value) {
        utm_precision_ = value;
    }

    void SetTopTimeWindow(double value) {
        top_time_window_ = value;
    }

    void SetBotTimeWindow(double value) {
        bot_time_window_ = value;
    }

    void SetTopDepthWindow(double value) {
        top_depth_window_ = value;
    }

    void SetBotDepthWindow(double value) {
        bot_depth_window_ = value;
    }

    void SetTimeWindowSpecified(bool value) {
        time_window_specified_ = value;
    }

    void SetDepthWindowSpecified(bool value) {
        depth_window_specified_ = value;
    }

    void SetMemoryLimit(double value) {
        memory_limit_ = value;
    }

    void SetTwtFileName(std::string name) {
        twt_file_name_ = name;
    }

    void SetPSSeismic(bool ps) {
        ps_seismic_ = ps;
    }

    void SetNMOCorr(bool nmo) {
        nmo_corr_ = nmo;
    }

    void SetVw(double value) {
        v_w_ = value;
    }

    void SetZw(double value) {
        z_w_ = value;
    }

    void SetOffset0(double value) {
      offset_0_ = value;
    }

    void SetDOffset(double value) {
      doffset_ = value;
    }

    void SetOffsetMax(double value) {
      offset_max_ = value;
    }



    bool GetOutputElasticParametersTimeSegy() {
        return elastic_parameters_time_segy_;
    }

    bool GetOutputElasticParametersDepthSegy() {
        return elastic_parameters_depth_segy_;
    }

    bool GetOutputExtraParametersTimeSegy() {
        return extra_parameters_time_segy_;
    }

    bool GetOutputExtraParametersDepthSegy() {
        return extra_parameters_depth_segy_;
    }

    int GetSegyInlineStart() {
        return inline_start_;
    }

    int GetSegyXlineStart() {
        return xline_start_;
    }

    std::string GetSegyInlineDirection() {
        return inline_direction_;
    }

    int GetSegyInlineStep() {
        return inline_step_;
    }

    int GetSegyXlineStep() {
        return xline_step_;
    }

    std::vector<std::string> GetExtraParameterNames() {
        return extra_parameter_names_;
    }

    std::vector<double> GetExtraParameterDefaultValues() {
        return extra_parameter_default_values_;
    }

    bool GetRicker() {
        return ricker_;
    }

    std::string GetWaveletFileFormat() {
        return wavelet_file_format_;
    }

    std::string GetWaveletFileName() {
        return wavelet_file_name_;
    }

    bool GetOutputSeismicStackTimeStorm() {
        return seismic_stack_time_storm_;
    }

    bool GetOutputSeismicStackTimeShiftStorm() {
        return seismic_stack_time_shift_storm_;
    }

    bool GetOutputSeismicStackDepthStorm() {
        return seismic_stack_depth_storm_;
    }

    bool GetOutputSeismicStackTimeSegy() {
        return seismic_stack_time_segy_;
    }

    bool GetOutputSeismicStackTimeShiftSegy() {
        return seismic_stack_time_shift_segy_;
    }

    bool GetOutputSeismicStackDepthSegy() {
        return seismic_stack_depth_segy_;
    }

    bool GetWhiteNoise() {
        return white_noise_;
    }

    double GetStandardDeviation() {
        return standard_deviation_;
    }

    unsigned long GetSeed() {
        return seed_;
    }

    ///--------------------------------------------------------------------------------------------------


    std::string GetEclipseFileName() {
        return eclipse_file_name_;
    }

    std::vector<std::string> GetParameterNames() {
        return parameter_names_;
    }

    std::vector<double> GetConstVp() {
        return constvp_;
    }

    std::vector<double> GetConstVs() {
        return constvs_;
    }

    std::vector<double> GetConstRho() {
        return constrho_;
    }

    double GetTheta0() {
        return theta_0_;
    }

    double GetDTheta() {
        return dtheta_;
    }

    double GetThetaMax() {
        return theta_max_;
    }

    double GetOffset0() {
        return offset_0_;
    }

    double GetDOffset() {
        return doffset_;
    }

    double GetOffsetMax() {
        return offset_max_;
    }

    double GetPeakFrequency() {
        return peak_f_;
    }

    double GetX0() {
        return x0_;
    }

    double GetY0() {
        return y0_;
    }

    double GetLx() {
        return lx_;
    }

    double GetLy() {
        return ly_;
    }

    double GetAngle() {
        return angle_;
    }

    //double GetConstantTop() { return constant_top_; }
    //double GetConstantBot() { return constant_bot_; }
    double GetDx() {
        return dx_;
    }

    double GetDy() {
        return dy_;
    }

    double GetDz() {
        return dz_;
    }

    double GetDt() {
        return dt_;
    }

    bool GetAreaGiven() {
        return area_given_;
    }

    double GetTopTimeConstant() {
        return top_time_constant_;
    }

    std::string GetTopTimeSurfaceFile() {
        return top_time_surface_;
    }

    double GetWaveletScale() {
        return wavelet_scale_;
    }

    bool GetOutputVp() {
        return output_vp_;
    }

    bool GetOutputReflections() {
        return output_reflections_;
    }

    bool GetOutputZvalues() {
        return output_zvalues_;
    }

    bool GetOutputSeismicDepth() {
        return output_seismic_depth_;
    }

    bool GetOutputSeismicTime() {
        return output_seismic_time_;
    }

    bool GetOutputSeismicTimeshift() {
        return output_seismic_timeshift_;
    }

    bool GetOutputDepthSurfaces() {
        return output_depth_surfaces_;
    }

    bool GetOutputTimeSurfaces() {
        return output_time_surfaces_;
    }

    bool GetOutputTwt() {
        return output_twt_;
    }

    std::string GetNLayersFileName() {
        return nlayers_file_name_;
    }

    std::string GetPrefix() {
        return prefix_;
    }

    std::string GetSuffix() {
        return suffix_;
    }

    double GetZeroThicknessLimit() {
        return zero_thickness_limit_;
    }

    bool GetOutputTimeSegy() {
        return output_time_segy_;
    }

    bool GetOutputTimeshiftSegy() {
        return output_timeshift_segy_;
    }

    bool GetOutputDepthSegy() {
        return output_depth_segy_;
    }

    bool GetUseCornerpointInterpol() {
        return use_cornerpoint_interpol_;
    }

    std::string GetAreaFromSurface() {
        return area_from_surface_;
    }

    int GetIL0In() {
        return il0_in_;
    }

    int GetXL0In() {
        return xl0_in_;
    }

    int GetUtmxIn() {
        return utmx_in_;
    }

    int GetUtmyIn() {
        return utmy_in_;
    }

    std::string GetAreaFromSegy() {
        return area_from_segy_;
    }

    short GetUtmPrecision() {
        return utm_precision_;
    }

    double GetTopTimeWindow() {
        return top_time_window_;
    }

    double GetBotTimeWindow() {
        return bot_time_window_;
    }

    double GetTopDepthWindow() {
        return top_depth_window_;
    }

    double GetBotDepthWindow() {
        return bot_depth_window_;
    }

    bool GetTimeWindowSpecified() {
        return time_window_specified_;
    }

    bool GetDepthWindowSpecified() {
        return depth_window_specified_;
    }

    double GetMemoryLimit() {
        return memory_limit_;
    }

    std::string GetTwtFileName() const {
        return twt_file_name_;
    }

    bool GetPSSeismic() {
        return ps_seismic_;
    }

    bool GetNMOCorr() {
        return nmo_corr_;
    }

    double GetVw() {
        return v_w_;
    }

    double GetZw() {
        return z_w_;
    }
  private:

    std::vector<double> constvp_;
    std::vector<double> constvs_;
    std::vector<double> constrho_;
    std::vector<std::string> parameter_names_;
    double theta_0_, dtheta_, theta_max_;  // seismic angle
    double offset_0_, doffset_, offset_max_;  // seismic offset
    //no default start
    std::string eclipse_file_name_;
    double x0_, y0_, lx_, ly_, angle_; // area parameters
    bool ricker_;
    bool area_given_;
    double peak_f_;
    bool topsurface_, botsurface_;
    //no default end
    std::string top_time_surface_;
    double top_time_constant_;
    double dx_, dy_, dz_, dt_; // cell size
    double wavelet_scale_;
    bool output_vp_, output_reflections_, output_zvalues_, output_seismic_time_, output_seismic_depth_, output_seismic_timeshift_;
   
    bool output_time_surfaces_, output_depth_surfaces_, output_twt_;
    std::string nlayers_file_name_;
    std::string prefix_, suffix_;
    double zero_thickness_limit_;
    bool output_time_segy_, output_depth_segy_, output_timeshift_segy_;
    bool use_cornerpoint_interpol_;
    std::string area_from_surface_;

    bool elastic_parameters_time_segy_, elastic_parameters_depth_segy_, extra_parameters_time_segy_, extra_parameters_depth_segy_;

    int inline_start_, xline_start_;
    std::string inline_direction_;
    int inline_step_, xline_step_;

    //no default start
    std::string wavelet_file_format_, wavelet_file_name_;
    std::vector<std::string> extra_parameter_names_;
    std::vector<double> extra_parameter_default_values_;
    //no default end
    bool seismic_stack_time_storm_, seismic_stack_time_shift_storm_, seismic_stack_depth_storm_, seismic_stack_time_segy_, seismic_stack_time_shift_segy_, seismic_stack_depth_segy_;

    bool white_noise_;
    double standard_deviation_;
    unsigned long seed_;

    int il0_in_, xl0_in_, utmx_in_, utmy_in_;
    std::string area_from_segy_;
    short utm_precision_;
    double memory_limit_;

    std::string twt_file_name_;

    bool ps_seismic_;

    bool nmo_corr_;
    double v_w_, z_w_;

    double top_time_window_, bot_time_window_, top_depth_window_, bot_depth_window_;
    bool time_window_specified_, depth_window_specified_;
};

#endif
