#ifndef SEISMIC_FORWARD_HPP
#define SEISMIC_FORWARD_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include "modelsettings.hpp"
#include <seismic_parameters.hpp>


#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/surface/regularsurface.hpp>
//#include "nrlib/geometry/interpolation.hpp"


class SeismicForward {
  public:
    static void seismicForward(SeismicParameters &seismic_parameters);
    static void printTime();
  private:

    static void makeNMOSeismic(SeismicParameters &seismic_parameters);
    static void makeSeismic(SeismicParameters &seismic_parameters);
    static void makeSeismicOLD(SeismicParameters &seismic_parameters);

    static void generateSeismic(std::vector<NRLib::StormContGrid> &rgridvec,
                                NRLib::StormContGrid &twtgrid,
                                NRLib::StormContGrid &zgrid,
                                NRLib::StormContGrid &twt_timeshift,
                                std::vector<NRLib::StormContGrid> &timegridvec,
                                std::vector<NRLib::StormContGrid> &depthgridvec,
                                std::vector<NRLib::StormContGrid> &timeshiftgridvec,
                                Wavelet *wavelet,
                                double dt,
                                NRLib::RegularSurface<double> &bot,
                                NRLib::RegularSurface<double> &toptime,
                                double t0, double dz, double z0,
                                std::vector<double> &constvp,
                                double waveletScale,
                                bool time_output,
                                bool depth_output,
                                bool timeshift_output);

    static void generateSeismicTrace(SeismicParameters             &seismic_parameters,
                                     const std::vector<double>     &twt_vec,
                                     const std::vector<double>     &twt_0,
                                     const std::vector<double>     &theta_vec,
                                     NRLib::Grid2D<double>         &timegrid_pos,
                                     size_t                         i,
                                     size_t                         j,
                                     unsigned long                  seed,
                                     const NRLib::StormContGrid    &zgrid,
                                     std::vector<NRLib::StormContGrid> &rgridvec);


    static void generateNMOSeismicTrace(SeismicParameters             &seismic_parameters,
                                        const std::vector<double>     &twt_vec,
                                        const std::vector<double>     &twt_0,
                                        const std::vector<double>     &offset_vec,
                                        NRLib::Grid2D<double>         &timegrid_pos,
                                        NRLib::Grid2D<double>         &nmo_timegrid_pos,
                                        NRLib::Grid2D<double>         &twtx_reg,
                                        size_t                         i,
                                        size_t                         j,
                                        unsigned long                  seed,
                                        const NRLib::StormContGrid    &zgrid,
                                        std::vector<NRLib::StormContGrid> &rgridvec,
                                        size_t                        &max_sample);

    static void seisConvolutionNMO(NRLib::Grid2D<double>               &timegrid_pos,
                                   NRLib::Grid2D<double>               &refl_pos,
                                   NRLib::Grid2D<double>               &twtx,
                                   const NRLib::StormContGrid          &zgrid,
                                   const NRLib::RegularSurface<double> &toptime,
                                   Wavelet                             *wavelet,
                                   double                               waveletScale,
                                   const std::vector<double>           &offset,
                                   double                               t0,
                                   double                               dt,
                                   size_t                               i,
                                   size_t                               j,
                                   const std::vector<size_t>           &n_min,
                                   const std::vector<size_t>           &n_max);

    static void seisConvolution(NRLib::Grid2D<double>               &timegrid_pos,
                                NRLib::Grid2D<double>               &refl_pos,
                                const std::vector<double>           &twt,
                                const NRLib::StormContGrid          &zgrid,
                                const NRLib::RegularSurface<double> &toptime,
                                Wavelet                             *wavelet,
                                double                               waveletScale,
                                const std::vector<double>           &theta_vec,
                                double                               t0,
                                double                               dt,
                                size_t                               i,
                                size_t                               j,
                                size_t           n_min,
                                  size_t           n_max);

    static void generateSeismicOnFile(std::vector<NRLib::StormContGrid> &rgridvec,
                                      NRLib::StormContGrid &twtgrid,
                                      NRLib::StormContGrid &zgrid,
                                      NRLib::StormContGrid &twt_timeshift,
                                      Wavelet *wavelet,
                                      double dt,
                                      int nt, int nz, int nx, int ny,
                                      NRLib::RegularSurface<double> &bot,
                                      NRLib::RegularSurface<double> &toptime,
                                      double t0, double dz, double z0,
                                      std::vector<double> &constvp,
                                      double waveletScale,
                                      bool time_output,
                                      bool depth_output,
                                      bool timeshift_output);

    static double findTFromZ(double z,
                             std::vector<double> &zvec,
                             std::vector<double> &tvec);

    static void   convertSeis(const std::vector<double>   &twt_vec,
                              const std::vector<double>   &twt_0,
                              const std::vector<double>   &zgrid_vec,
                              const std::vector<double>   &z_0,
                              const NRLib::Grid2D<double> &seismic,
                              NRLib::Grid2D<double>       &conv_seismic,
                              const size_t                &max_sample);


    static void   NMOCorrect(const std::vector<double>   &t_in,
                             const NRLib::Grid2D<double> &data_in,
                             const NRLib::Grid2D<double> &t_out,
                             NRLib::Grid2D<double>       &data_out,
                             const std::vector<size_t>   &n_min,
                             const std::vector<size_t>   &n_max,
                             size_t                      &max_sample);

    static void   findNMOTheta(NRLib::Grid2D<double>     &thetagrid,
                               const std::vector<double> &twt_vec,
                               const std::vector<double> &vrms_vec,
                               const std::vector<double> &offset);


    static void   findTWTx(NRLib::Grid2D<double>     &twtx_grid,
                           const std::vector<double> &twt_vec,
                           const std::vector<double> &vrms_vec,
                           const std::vector<double> &offset);



    static bool   generateTraceOk(SeismicParameters &seismic_parameters,
                                  size_t i,
                                  size_t j);

    static void   extrapolZandTwtVec(std::vector<double>        &zgrid_vec_extrapol,
                                     std::vector<double>        &twt_vec_extrapol,
                                     const std::vector<double>  &twt_vec,
                                     const NRLib::StormContGrid &zgrid,
                                     double                      z_bot,
                                     double                      vp_bot,
                                     double                      vs_bot,
                                     size_t                      i,
                                     size_t                      j,
                                     bool                        ps_seis);

    static std::vector<double>   linInterp1D(const std::vector<double> &x_in,
                                             const std::vector<double> &y_in,
                                             const std::vector<double> &x_out);


    static std::vector<double>   splineInterp1D(const std::vector<double> &x_in,
                                                const std::vector<double> &y_in,
                                                const std::vector<double> &x_out,
                                                double                     extrap_value);
    static void monitorInitialize(size_t nx,
                                  size_t ny,
                                  float &monitor_size,
                                  float &next_monitor);

    static void monitor(size_t n_xl,
                        size_t il_steps,
                        size_t xl_steps,
                        float monitor_size,
                        float &next_monitor);



};

#endif
