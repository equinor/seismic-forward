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

    ~SeismicParameters() {
    };


    NRLib::StormContGrid &vpGrid()                             { return *vpgrid_; };
    NRLib::StormContGrid &vsGrid()                             { return *vsgrid_; };
    NRLib::StormContGrid &rhoGrid()                            { return *rhogrid_; };
    NRLib::StormContGrid &zGrid()                              { return *zgrid_; };
    NRLib::StormContGrid &twtGrid()                            { return *twtgrid_; };
    NRLib::StormContGrid &twtShiftGrid()                       { return *twt_timeshift_; };
    std::vector<NRLib::StormContGrid> &rGrids()                { return *rgridvec_; };
    std::vector<NRLib::StormContGrid> &extraParametersGrids()  { return *extra_parameter_grid_; };
    NRLib::EclipseGrid &eclipseGrid()                          { return *eclipse_grid_; };

    NRLib::StormContGrid &vrmsGrid()                           { return *vrmsgrid_; };


    void deleteEclipseGrid();
    void deleteElasticParameterGrids();
    void deleteExtraParameterGrids();
    void deleteZandRandTWTGrids();
    void deleteVrmsGrid();
    void deleteWavelet();
    void deleteGeometryAndOutput();

    size_t topK()                    { return top_k_; }
    size_t bottomK()                 { return bottom_k_; }
    double theta0()                  { return theta_0_; }
    double dTheta()                  { return dtheta_; }
    size_t nTheta()                  { return ntheta_; }
    std::vector<double> & GetThetaVec() { return theta_vec_; }
    double offset0()                 { return offset_0_; }
    double dOffset()                 { return doffset_; }
    size_t nOffset()                 { return noffset_; }
    std::vector<double> & GetOffsetVec() { return offset_vec_; }
    

    void   findLoopIndeces(int               &n_xl,
                           int               &il_min,
                           int               &il_max,
                           int               &il_step,
                           int               &xl_min,
                           int               &xl_max,
                           int               &xl_step,
                           bool              &segy);

    std::vector<double> GenerateTwt0ForNMO(size_t & time_stretch_samples);
    std::vector<double> GenerateZ0ForNMO();

    std::vector<double> GenerateTWT0Shift(double twt_0_min,
                                          size_t n_samples);

    void getSeisLimits(std::vector<double>  twt_0,
                       std::vector<double>  vrms_vec,
                       std::vector<double>  offset_vec,
                       std::vector<size_t> &n_min,
                       std::vector<size_t> &n_max);


   void   findVrmsPos(std::vector<double>       &vrms_vec,
                       std::vector<double>       &vrms_vec_reg,
                       const std::vector<double> &twt_0,
                       size_t                    i,
                       size_t                    j,
                       bool                      include_regular = true);

   void   findNMOReflections(NRLib::Grid2D<double>       &r_vec,
                             const NRLib::Grid2D<double> &theta_vec,
                             const std::vector<double>   &offset_vec,
                             size_t                       i,
                             size_t                       j);

   void   findReflections(NRLib::Grid2D<double>       &r_vec,
                          const std::vector<double>   &theta_vec,
                          size_t                       i,
                          size_t                       j);

    NRLib::RegularSurface<double> &topTime()       { return top_time_; };
    NRLib::RegularSurface<double> &bottomTime()    { return bot_time_; };

    NRLib::RegularSurface<double> &topEclipse()    { return topeclipse_; };
    NRLib::RegularSurface<double> &bottomEclipse() { return boteclipse_; };


    ModelSettings*       modelSettings()   { return model_settings_;};

    SeismicOutput*       seismicOutput()   { return seismic_output_; };
    SeismicGeometry*     seismicGeometry() { return seismic_geometry_; };


    Wavelet*             wavelet()         { return wavelet_; };
    double               waveletScale()    { return wavelet_scale_; };
    NRLib::SegyGeometry* segyGeometry()    { return segy_geometry_; };

    void setSegyGeometry(const NRLib::SegyGeometry &geometry);

    bool GetTimeOutput();
    bool GetDepthOutput();
    bool GetTimeshiftOutput();
    bool GetStackOutput();
    bool GetSegyOutput();
    bool GetTimeStormOutput();
    bool GetDepthStormOutput();
    bool GetTimeshiftStormOutput();
    bool GetStormOutput();

  private:
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
    double wavelet_scale_;

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
    std::vector<NRLib::StormContGrid> *rgridvec_;
    std::vector<NRLib::StormContGrid> *extra_parameter_grid_;

    NRLib::StormContGrid *vrmsgrid_;

    NRLib::StormContGrid *twt_timeshift_;


    std::vector<double> twt_0_;
    std::vector<double> z_0_;

    void setupWavelet();

    void readEclipseGrid();

    void findGeometry();

    void findSurfaceGeometry();

    void calculateAngleSpan();
    void calculateOffsetSpan();

    void createGrids();
};

#endif
