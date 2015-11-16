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


    NRLib::StormContGrid &vpGrid()                             { return *vpgrid_; };
    NRLib::StormContGrid &vsGrid()                             { return *vsgrid_; };
    NRLib::StormContGrid &rhoGrid()                            { return *rhogrid_; };
    NRLib::StormContGrid &zGrid()                              { return *zgrid_; };
    NRLib::StormContGrid &twtGrid()                            { return *twtgrid_; };
    NRLib::StormContGrid &twtSSGrid()                          { return *twtssgrid_; };
    NRLib::StormContGrid &twtPPGrid()                          { return *twtppgrid_; };
    NRLib::StormContGrid &twtShiftGrid()                       { return *twt_timeshift_; };
    NRLib::StormContGrid &vrmsGrid()                           { return *vrmsgrid_; };
    std::vector<NRLib::StormContGrid> &rGrids()                { return *rgridvec_; };
    std::vector<NRLib::StormContGrid> &extraParametersGrids()  { return *extra_parameter_grid_; };
    NRLib::EclipseGrid &eclipseGrid()                          { return *eclipse_grid_; };

    NRLib::RegularSurface<double> &topTime()                   { return top_time_; };
    NRLib::RegularSurface<double> &bottomTime()                { return bot_time_; };
    NRLib::RegularSurface<double> &topEclipse()                { return topeclipse_; };
    NRLib::RegularSurface<double> &bottomEclipse()             { return boteclipse_; };
    
    ModelSettings*       modelSettings()                       { return model_settings_;};
    SeismicOutput*       seismicOutput()                       { return seismic_output_; };
    SeismicGeometry*     seismicGeometry()                     { return seismic_geometry_; };
    Wavelet*             wavelet()                             { return wavelet_; };
    double               waveletScale() const                  { return wavelet_scale_; };
    NRLib::SegyGeometry* segyGeometry()                        { return segy_geometry_; };

    size_t topK()                                              { return top_k_; }
    size_t bottomK()                                           { return bottom_k_; }
    double theta0()                                            { return theta_0_; }
    double dTheta()                                            { return dtheta_; }
    size_t nTheta()                                            { return ntheta_; }
    std::vector<double> & GetThetaVec()                        { return theta_vec_; }
    double offset0()                                           { return offset_0_; }
    double dOffset()                                           { return doffset_; }
    size_t nOffset()                                           { return noffset_; }
    std::vector<double> & GetOffsetVec()                       { return offset_vec_; }
    float GetMissingVal()                                      { return missing_; }

    bool GetTimeOutput();
    bool GetDepthOutput();
    bool GetTimeshiftOutput();
    bool GetStackOutput();
    bool GetSegyOutput();
    bool GetTimeStormOutput();
    bool GetDepthStormOutput();
    bool GetTimeshiftStormOutput();
    bool GetStormOutput();

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

#endif
