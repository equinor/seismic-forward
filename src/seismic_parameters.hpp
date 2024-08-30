#ifndef SEISMIC_PARAMETERS_HPP
#define SEISMIC_PARAMETERS_HPP

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/volume/volume.hpp"

#include "utils/trace.hpp"

#include "seismic_output.hpp"
#include "modelsettings.hpp"

#include <tbb/concurrent_queue.h>

#include <stdio.h>
#include <string>
#include <vector>

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

class SeismicParameters
{
  public:
    SeismicParameters(ModelSettings *model_settings);

    ~SeismicParameters() {};

    inline ModelSettings                      * GetModelSettings()        const { return model_settings_       ;}

    inline NRLib::StormContGrid               & GetVpGrid()               const { return *vpgrid_              ;}
    inline NRLib::StormContGrid               & GetVsGrid()               const { return *vsgrid_              ;}
    inline NRLib::StormContGrid               & GetRhoGrid()              const { return *rhogrid_             ;}
    inline NRLib::StormContGrid               & GetZGrid()                const { return *zgrid_               ;}
    inline NRLib::StormContGrid               & GetTwtGrid()              const { return *twtgrid_             ;}
    inline NRLib::StormContGrid               & GetTwtSSGrid()            const { return *twtssgrid_           ;}
    inline NRLib::StormContGrid               & GetTwtPPGrid()            const { return *twtppgrid_           ;}
    inline NRLib::StormContGrid               & GetTwtShiftGrid()         const { return *twt_timeshift_       ;}
    inline NRLib::StormContGrid               & GetVrmsGrid()             const { return *vrmsgrid_            ;}
    inline std::vector<NRLib::StormContGrid>  & GetRGrids()               const { return *rgridvec_            ;}
    inline std::vector<NRLib::StormContGrid*>   GetExtraParametersGrids() const { return  extra_parameter_grid_;}
    inline NRLib::EclipseGrid                 & GetEclipseGrid()          const { return *eclipse_grid_        ;}

    inline NRLib::RegularSurface<double>      & GetTopTime()                    { return  top_time_            ;}
    inline NRLib::RegularSurface<double>      & GetBottomTime()                 { return  bot_time_            ;}
    inline NRLib::RegularSurface<double>      & GetTopEclipse()                 { return  topeclipse_          ;}
    inline NRLib::RegularSurface<double>      & GetBottomEclipse()              { return  boteclipse_          ;}

    inline SeismicOutput                      * GetSeismicOutput()        const { return  seismic_output_      ;}
    inline SeismicGeometry                    * GetSeismicGeometry()      const { return  seismic_geometry_    ;}
    inline Wavelet                            * GetWavelet()              const { return  wavelet_             ;}
    inline double                               GetWaveletScale()         const { return  wavelet_scale_       ;}
    inline NRLib::SegyGeometry                * GetSegyGeometry()         const { return  segy_geometry_       ;}

    inline size_t                               GetTopK()                 const { return  top_k_               ;}
    inline size_t                               GetBottomK()              const { return  bottom_k_            ;}
    inline std::vector<double>                & GetThetaVec()                   { return  theta_vec_           ;}
    inline std::vector<double>                & GetOffsetVec()                  { return  offset_vec_          ;}
    inline float                                GetMissingVal()           const { return  missing_             ;}

    void SetSegyGeometry(const NRLib::SegyGeometry * geometry);

    void FindLoopIndeces(int  & n_xl,
                         int  & il_min,
                         int  & il_max,
                         int  & il_step,
                         int  & xl_min,
                         int  & xl_max,
                         int  & xl_step,
                         bool & segy);

    void FindMaxTwtIndex(size_t & i_max,
                         size_t & j_max,
                         double & max_value);

    void GenerateTwt0AndZ0(ModelSettings       * model_settings,
                           std::vector<double> & twt_0,
                           std::vector<double> & z_0,
                           std::vector<double> & twts_0,
                           size_t              & time_samples_stretch,
                           bool                  ps_seis);

    void GenerateTwt0ForNMO(std::vector<double> & twt_0,
                            size_t              & nt_stretch,
                            double              & stretch_factor,
                            const bool            ps_seis,
                            const double          twt_wavelet,
                            const size_t          nt,
                            const double          dt,
                            const double          t0,
                            const size_t          nzrefl);

    void GenerateZ0ForNMO(std::vector<double> & z_0,
                          const double          stretch_factor);

    std::vector<double> GenerateTWT0Shift(double twt_0_min,
                                          size_t n_samples);

    static void FindPSNMOThetaAndOffset(NRLib::Grid2D<double>     & thetagrid,
                                        NRLib::Grid2D<double>     & offset_down_grid,
                                        NRLib::Grid2D<double>     & offset_up_grid,
                                        const std::vector<double> & twt_pp_vec,
                                        const std::vector<double> & twt_ps_vec,
                                        const std::vector<double> & vrms_pp_vec,
                                        const std::vector<double> & vrms_ss_vec,
                                        const std::vector<double> & offset,
                                        bool                        save_theta = true);

    static double FindSinThetaPSWithNewtonsMethod(double   start_value,
                                                  double   offset,
                                                  double   dU,
                                                  double   dD,
                                                  double   vr,
                                                  double   tol,
                                                  size_t & n_it);

    void FindVrms(std::vector<double>       & vrms_vec,
                  const std::vector<double> & twt_vec,
                  const std::vector<double> & v_vec,
                  double                      z_res) const;

    void FindVrmsRegular(const std::vector<double> & vrms_vec,
                         std::vector<double>       & vrms_vec_reg,
                         const std::vector<double> & twt_vec,
                         const std::vector<double> & twt_vec_reg,
                         const std::vector<double> & v_vec,
                         double                      const_v,
                         double                      twt_wavelet_exstrapol) const;

    void  FindReflections(NRLib::Grid2D<double>       & r_vec,
                          const NRLib::Grid2D<double> & theta_vec,
                          size_t                        i,
                          size_t                        j);

    static void PrintElapsedTime(time_t      start_time,
                                 std::string work);

    tbb::concurrent_queue<Trace*> FindTracesInForward(size_t & n_traces);


    static void MonitorInitialize(size_t   n_traces,
                                  float  & monitor_size,
                                  float  & next_monitor);

    static void Monitor(size_t   trace,
                        float    monitor_size,
                        float  & next_monitor);

    static std::vector<double> LinInterp1D(const std::vector<double> & x_in,
                                           const std::vector<double> & y_in,
                                           const std::vector<double> & x_out);


    static std::vector<double> SplineInterp1D(const std::vector<double> & x_in,
                                              const std::vector<double> & y_in,
                                              const std::vector<double> & x_out,
                                              double                      extrap_value);


    static void FindExtrapolationRegion(NRLib::Grid2D<bool> & extrapolate,
                                        SeismicGeometry     & seismic_geometry,
                                        double                xmin,
                                        double                ymin,
                                        double                etdx,
                                        double                etdy);

    void DeleteEclipseGrid();
    void DeleteElasticParameterGrids();
    void DeleteExtraParameterGrids();
    void DeleteZandRandTWTGrids();
    void DeleteVrmsGrid();
    void DeleteWavelet();
    void DeleteGeometryAndOutput();

private:

  void SetupWavelet(Wavelet           *& wavelet,
                    const bool           use_zero_time_from_header,
                    const bool           write_wavelet,
                    const bool           use_ricker,
                    const double         peakF,
                    const double         length,
                    const double         length_factor,
                    const double         dt,
                    const std::string  & file_name,
                    const std::string  & file_format,
                    const std::string  & prefix);

  void ReadEclipseGrid(NRLib::EclipseGrid             *& eclipse_grid,
                       const std::string               & filename,
                       const std::vector<std::string>  & names,
                       const std::vector<std::string>  & extra_parameter_names);

  void FindGeometry(SeismicGeometry              *& seismic_geometry,
                    NRLib::SegyGeometry          *& segy_geometry,
                    const NRLib::EclipseGeometry  & eclipse_geometry,
                    ModelSettings                 * model_settings);

  void FindTopAndBaseSurfaces(NRLib::RegularSurface<double> & top_time,
                              NRLib::RegularSurface<double> & bot_time,
                              NRLib::RegularSurface<double> & topeclipse,
                              NRLib::RegularSurface<double> & boteclipse,
                              size_t                        & top_k,
                              size_t                        & bot_k,
                              SeismicGeometry               * seismic_geometry,
                              const NRLib::EclipseGeometry  & eclipse_geometry,
                              Wavelet                       * wavelet,
                              ModelSettings                 * model_settings);

  void CreateGrids(SeismicGeometry * seismic_geometry,
                   ModelSettings   * model_settings);

  ModelSettings                      * model_settings_;
  SeismicGeometry                    * seismic_geometry_;
  SeismicOutput                      * seismic_output_;

  std::vector<double>                  theta_vec_;
  std::vector<double>                  offset_vec_;

  Wavelet                            * wavelet_;
  double                               wavelet_scale_;

  size_t                               top_k_;
  size_t                               bottom_k_;

  NRLib::EclipseGrid                 * eclipse_grid_;
  NRLib::SegyGeometry                * segy_geometry_;

  NRLib::RegularSurface<double>        top_time_;
  NRLib::RegularSurface<double>        bot_time_;
  NRLib::RegularSurface<double>        topeclipse_;
  NRLib::RegularSurface<double>        boteclipse_;

  NRLib::StormContGrid               * zgrid_;
  NRLib::StormContGrid               * vpgrid_;
  NRLib::StormContGrid               * vsgrid_;
  NRLib::StormContGrid               * rhogrid_;
  NRLib::StormContGrid               * twtgrid_;
  NRLib::StormContGrid               * twtssgrid_;
  NRLib::StormContGrid               * twtppgrid_;
  NRLib::StormContGrid               * twt_timeshift_;
  NRLib::StormContGrid               * vrmsgrid_;
  std::vector<NRLib::StormContGrid>  * rgridvec_;
  std::vector<NRLib::StormContGrid*>   extra_parameter_grid_;

  float                                missing_;
};

#endif
