#ifndef SEISMIC_PARAMETERS_HPP
#define SEISMIC_PARAMETERS_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include <nrlib/surface/regularsurface.hpp>
#include <nrlib/volume/volume.hpp>
#include "modelsettings.hpp"
#include "seismic_output.hpp"

class Wavelet;
class ModelSettings;
class SeismicGeometry;
class SeismicOutput;

namespace NRLib {
    class EclipseGrid;
    class EclipseGeometry;
    class SegyGeometry;
    class StormContGrid;
}


class SeismicParameters {
  public:
    SeismicParameters(ModelSettings *model_settings);

    ~SeismicParameters() {};


    inline NRLib::StormContGrid &GetVpGrid()                            const { return *vpgrid_; };
    inline NRLib::StormContGrid &GetVsGrid()                            const { return *vsgrid_; };
    inline NRLib::StormContGrid &GetRhoGrid()                           const { return *rhogrid_; };
    inline NRLib::StormContGrid &GetZGrid()                             const { return *zgrid_; };
    inline NRLib::StormContGrid &GetTwtGrid()                           const { return *twtgrid_; };
    inline NRLib::StormContGrid &GetTwtSSGrid()                         const { return *twtssgrid_; };
    inline NRLib::StormContGrid &GetTwtPPGrid()                         const { return *twtppgrid_; };
    inline NRLib::StormContGrid &GetTwtShiftGrid()                      const { return *twt_timeshift_; };
    inline NRLib::StormContGrid &GetVrmsGrid()                          const { return *vrmsgrid_; };
    inline std::vector<NRLib::StormContGrid> &GetRGrids()               const { return *rgridvec_; };
    inline std::vector<NRLib::StormContGrid> &GetExtraParametersGrids() const { return *extra_parameter_grid_; };
    inline NRLib::EclipseGrid &GetEclipseGrid()                         const { return *eclipse_grid_; };

    inline NRLib::RegularSurface<double> &GetTopTime()                        { return top_time_; };
    inline NRLib::RegularSurface<double> &GetBottomTime()                     { return bot_time_; };
    inline NRLib::RegularSurface<double> &GetTopEclipse()                     { return topeclipse_; };
    inline NRLib::RegularSurface<double> &GetBottomEclipse()                  { return boteclipse_; };

    inline ModelSettings*       GetModelSettings()                      const { return model_settings_;};
    inline SeismicOutput*       GetSeismicOutput()                      const { return seismic_output_; };
    inline SeismicGeometry*     GetSeismicGeometry()                    const { return seismic_geometry_; };
    inline Wavelet*             GetWavelet()                            const { return wavelet_; };
    inline double               GetWaveletScale()                       const { return wavelet_scale_; };
    inline NRLib::SegyGeometry* GetSegyGeometry()                       const { return segy_geometry_; };

    inline size_t               GetTopK()                               const { return top_k_;      }
    inline size_t               GetBottomK()                            const { return bottom_k_;   }
    inline std::vector<double> &GetThetaVec()                                 { return theta_vec_;  }
    inline std::vector<double> &GetOffsetVec()                                { return offset_vec_; }
    inline float                GetMissingVal()                         const { return missing_;    }

    //inline double               GetTheta0()                             const { return theta_0_;  }
    //inline double               GetDTheta()                             const { return dtheta_;   }
    //inline size_t               GetNTheta()                             const { return ntheta_;   }
    //inline double               GetOffset0()                            const { return offset_0_; }
    //inline double               GetDOffset()                            const { return doffset_;  }
    //inline size_t               GetNOffset()                            const { return noffset_;  }

    inline bool GetTimeOutput()           const;
    inline bool GetDepthOutput()          const;
    inline bool GetTimeshiftOutput()      const;
    inline bool GetStackOutput()          const;
    inline bool GetSegyOutput()           const;
    inline bool GetTimeStormOutput()      const;
    inline bool GetDepthStormOutput()     const;
    inline bool GetTimeshiftStormOutput() const;
    inline bool GetStormOutput()          const;

    void SetSegyGeometry(const NRLib::SegyGeometry &geometry);

    void FindLoopIndeces(int               &n_xl,
                         int               &il_min,
                         int               &il_max,
                         int               &il_step,
                         int               &xl_min,
                         int               &xl_max,
                         int               &xl_step,
                         bool              &segy);

    void FindMaxTwtIndex(size_t &i_max, 
                         size_t &j_max, 
                         double &max_value);

    void GenerateTwt0AndZ0(std::vector<double> &twt_0,
                           std::vector<double> &z_0,
                           std::vector<double> &twts_0,
                           size_t              &time_samples_stretch,
                           bool                 ps_seis);

    std::vector<double> GenerateTwt0ForNMO(size_t &time_stretch_samples, 
                                           bool    ps_seis);

    std::vector<double> GenerateZ0ForNMO();

    std::vector<double> GenerateTWT0Shift(double twt_0_min,
                                          size_t n_samples);

    static void FindPSNMOThetaAndOffset(NRLib::Grid2D<double>     &thetagrid,
                                        NRLib::Grid2D<double>     &offset_down_grid,
                                        NRLib::Grid2D<double>     &offset_up_grid,
                                        const std::vector<double> &twt_pp_vec,
                                        const std::vector<double> &twt_ps_vec,
                                        const std::vector<double> &vrms_pp_vec,
                                        const std::vector<double> &vrms_ss_vec,
                                        const std::vector<double> &offset,
                                        bool                      save_theta = true);

    static double FindSinThetaPSWithNewtonsMethod(double start_value,
                                                  double offset,
                                                  double dU,
                                                  double dD,
                                                  double vr,
                                                  double tol,
                                                  size_t &n_it);

    void FindVrms(std::vector<double>       &vrms_vec,
                  std::vector<double>       &vrms_vec_reg,
                  const std::vector<double> &twt_vec,
                  const std::vector<double> &twt_vec_reg,
                  const std::vector<double> &v_vec,
                  double                     const_v,
                  double                     twt_wavelet_exstrapol,
                  size_t                     i,
                  size_t                     j,
                  bool                       include_regular) const;

    void  FindNMOReflections(NRLib::Grid2D<double>       &r_vec,
                             const NRLib::Grid2D<double> &theta_vec,
                             size_t                       i,
                             size_t                       j);
    
    void  FindReflections(NRLib::Grid2D<double>       &r_vec,
                          const std::vector<double>   &theta_vec,
                          size_t                       i,
                          size_t                       j);

    static void PrintElapsedTime(time_t start_time, std::string work);

    void DeleteEclipseGrid();
    void DeleteElasticParameterGrids();
    void DeleteExtraParameterGrids();
    void DeleteZandRandTWTGrids();
    void DeleteVrmsGrid();
    void DeleteWavelet();
    void DeleteGeometryAndOutput();

  private:
    void SetupWavelet();
    void ReadEclipseGrid();
    void FindGeometry();
    void FindSurfaceGeometry();
    void CalculateAngleSpan();
    void CalculateOffsetSpan();
    void CreateGrids();

    ModelSettings *model_settings_;
    SeismicGeometry *seismic_geometry_;
    SeismicOutput *seismic_output_;

    size_t ntheta_;
    double theta_0_;
    double dtheta_;
    double theta_max_;
    std::vector<double> theta_vec_;

    size_t noffset_;
    double offset_0_;
    double doffset_;
    double offset_max_;
    std::vector<double> offset_vec_;

    Wavelet *wavelet_;
    double   wavelet_scale_;

    NRLib::EclipseGrid *eclipse_grid_;

    size_t top_k_;
    size_t bottom_k_;

    NRLib::RegularSurface<double> top_time_;
    NRLib::RegularSurface<double> bot_time_;
    NRLib::RegularSurface<double> topeclipse_;
    NRLib::RegularSurface<double> boteclipse_;

    NRLib::SegyGeometry *segy_geometry_;

    NRLib::StormContGrid *zgrid_;
    NRLib::StormContGrid *vpgrid_;
    NRLib::StormContGrid *vsgrid_;
    NRLib::StormContGrid *rhogrid_;
    NRLib::StormContGrid *twtgrid_;
    NRLib::StormContGrid *twtssgrid_;
    NRLib::StormContGrid *twtppgrid_;
    NRLib::StormContGrid *twt_timeshift_;
    NRLib::StormContGrid *vrmsgrid_;
    std::vector<NRLib::StormContGrid> *rgridvec_;
    std::vector<NRLib::StormContGrid> *extra_parameter_grid_;   


    std::vector<double> twt_0_;
    std::vector<double> z_0_;

    float missing_;
};

bool SeismicParameters::GetTimeOutput() const {
  return (model_settings_->GetOutputTimeSegy()
    || model_settings_->GetOutputSeismicStackTimeSegy()
    || model_settings_->GetOutputSeismicTime()
    || model_settings_->GetOutputSeismicStackTimeStorm()
    || model_settings_->GetOutputPrenmoTimeSegy());
}

bool SeismicParameters::GetDepthOutput() const {
  return (model_settings_->GetOutputDepthSegy()
    || model_settings_->GetOutputSeismicStackDepthSegy()
    || model_settings_->GetOutputSeismicDepth()
    || model_settings_->GetOutputSeismicStackDepthStorm());
}

bool SeismicParameters::GetTimeshiftOutput() const {
  return (model_settings_->GetOutputTimeshiftSegy()
    || model_settings_->GetOutputSeismicStackTimeShiftSegy()
    || model_settings_->GetOutputSeismicTimeshift()
    || model_settings_->GetOutputSeismicStackTimeShiftStorm());
}
bool SeismicParameters::GetStackOutput() const {
  return (model_settings_->GetOutputSeismicStackTimeSegy()
    || model_settings_->GetOutputSeismicStackDepthSegy()
    || model_settings_->GetOutputSeismicStackTimeShiftSegy()
    || model_settings_->GetOutputSeismicStackTimeStorm()
    || model_settings_->GetOutputSeismicStackDepthStorm()
    || model_settings_->GetOutputSeismicStackTimeShiftStorm());
}
bool SeismicParameters::GetSegyOutput() const {
  return (model_settings_->GetOutputTimeSegy()
    || model_settings_->GetOutputSeismicStackTimeSegy()
    || model_settings_->GetOutputDepthSegy()
    || model_settings_->GetOutputSeismicStackDepthSegy()
    || model_settings_->GetOutputTimeshiftSegy()
    || model_settings_->GetOutputSeismicStackTimeShiftSegy()
    || model_settings_->GetOutputPrenmoTimeSegy());
}
bool SeismicParameters::GetTimeStormOutput() const {
  return (model_settings_->GetOutputSeismicStackTimeStorm() || model_settings_->GetOutputSeismicTime());
}
bool SeismicParameters::GetDepthStormOutput() const {
  return (model_settings_->GetOutputSeismicStackDepthStorm() || model_settings_->GetOutputSeismicDepth());
}
bool SeismicParameters::GetTimeshiftStormOutput() const {
  return (model_settings_->GetOutputSeismicStackTimeShiftStorm() || model_settings_->GetOutputSeismicTimeshift());
}

bool SeismicParameters::GetStormOutput() const {
  return (GetTimeStormOutput() || GetDepthStormOutput() || GetTimeshiftStormOutput());
}


#endif
