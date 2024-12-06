#ifndef SEISMIC_FORWARD_HPP
#define SEISMIC_FORWARD_HPP

#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/random/randomgenerator.hpp"
#include "nrlib/random/normal.hpp"

#include "utils/gen_seis_trace_params.hpp"
#include "utils/output.hpp"
#include "utils/trace.hpp"

#include "seismic_parameters.hpp"
#include "modelsettings.hpp"

#include <tbb/concurrent_queue.h>

#include <algorithm>
#include <stdio.h>
#include <string>
#include <vector>

class SeismicForward {
  public:

    static void DoSeismicForward(SeismicParameters & seismic_parameters,
                                 ModelSettings     * model_settings);

  private:

    static void MakeReflections(NRLib::Grid2D<double>             & refl_pos,
                                std::vector<NRLib::StormContGrid> & rgridvec,
                                SeismicParameters                 & seismic_parameters,
                                NRLib::Grid2D<double>             & theta,
                                bool                                output_refl,
                                bool                                add_noise,
                                double                              std,
                                unsigned long                       seed,
                                size_t                              nx,
                                size_t                              i,
                                size_t                              j);

    static void MakeNMOSeismic(SeismicParameters & seismic_parameters,
                               ModelSettings     * model_settings);

    static void MakeSeismic(SeismicParameters & seismic_parameters,
                            ModelSettings     * model_settings);

    static void GenerateSeismicTraces(SeismicParameters         & seismic_parameters,
                                      const std::vector<double> & twt_0,
                                      const std::vector<double> & z_0,
                                      const std::vector<double> & twts_0,
                                      const std::vector<double> & theta_vec,
                                      const std::vector<double> & offset_vec,
                                      Output                    * output,
                                      Trace                     * trace,
                                      ResultTrace               *& result_trace);

    static void GenerateNMOSeismicTraces(SeismicParameters         & seismic_parameters,
                                         const std::vector<double> & twt_0,
                                         const std::vector<double> & z_0,
                                         const std::vector<double> & twts_0,
                                         const std::vector<double> & theta_vec,
                                         const std::vector<double> & offset_vec,
                                         const size_t                time_samples_stretch,
                                         Output                    * nmo_output,
                                         Trace                     * trace,
                                         ResultTrace               *& result_trace);

    static bool GenerateTraceOk(SeismicParameters & seismic_parameters,
                                ModelSettings     * model_settings,
                                size_t              i,
                                size_t              j);


    static void SeisConvolutionNMO(NRLib::Grid2D<double>               & timegrid_pos,
                                   NRLib::Grid2D<double>               & refl_pos,
                                   NRLib::Grid2D<double>               & twtx,
                                   const NRLib::StormContGrid          & zgrid,
                                   const NRLib::RegularSurface<double> & toptime,
                                   Wavelet                             * wavelet,
                                   double                                waveletScale,
                                   const std::vector<double>           & offset,
                                   double                                t0,
                                   double                                dt,
                                   size_t                                i,
                                   size_t                                j,
                                   const std::vector<size_t>           & n_min,
                                   const std::vector<size_t>           & n_max);

    static void SeisConvolution(NRLib::Grid2D<double>               & timegrid_pos,
                                NRLib::Grid2D<double>               & refl_pos,
                                const std::vector<double>           & twt,
                                const NRLib::StormContGrid          & zgrid,
                                const NRLib::RegularSurface<double> & toptime,
                                Wavelet                             * wavelet,
                                double                                waveletScale,
                                const std::vector<double>           & theta_vec,
                                double                                t0,
                                double                                dt,
                                size_t                                i,
                                size_t                                j,
                                size_t                                n_min,
                                size_t                                n_max);

    static void ConvertShiftAndStack(NRLib::Grid2D<double>      & nmo_timegrid_stack_pos,
                                     NRLib::Grid2D<double>      & nmo_timeshiftgrid_stack_pos,
                                     NRLib::Grid2D<double>      & nmo_depthgrid_stack_pos,
                                     NRLib::Grid2D<double>      & nmo_depthgrid_pos,
                                     NRLib::Grid2D<double>      & nmo_timeshiftgrid_pos,
                                     NRLib::Grid2D<double>      & nmo_timegrid_pos,
                                     const NRLib::StormContGrid & twt_timeshift,
                                     const std::vector<double>  & twt_vec,
                                     const std::vector<double>  & twts_0,
                                     const std::vector<double>  & twt_0,
                                     const NRLib::StormContGrid & zgrid,
                                     const std::vector<double>  & z_0,
                                     const double                 vp_bot,
                                     const double                 vs_bot,
                                     const double                 z_wavelet_bot_stretched,
                                     const double                 twt_wavelet_stretched,
                                     const double                 sd,
                                     const int                    tshift,
                                     const int                    zshift,
                                     const bool                   depth_conversion,
                                     const bool                   time_shift,
                                     const bool                   add_white_noise,
                                     const bool                   equal_noise,
                                     const bool                   ps_seis,
                                     const bool                   stack_time,
                                     const bool                   stack_timeshift,
                                     const bool                   stack_depth,
                                     const size_t                 nzrefl,
                                     const size_t                 noff,
                                     const unsigned long          seed,
                                     const size_t                 n_traces,
                                     const size_t                 nt_non_nmo,
                                     const size_t                 nz_non_nmo,
                                     const size_t                 i,
                                     const size_t                 j);

    static void GenerateWhiteNoiseAllOffsets(std::vector<std::vector<double> > & noise,
                                             bool                                equal_noise,
                                             unsigned long                       seed,
                                             double                              sd,
                                             size_t                              n_traces,
                                             size_t                              nt);

    static std::vector<double> GenerateWhiteNoise(unsigned long seed,
                                                  double        sd,
                                                  size_t        n);

    static void StackOffsets(NRLib::Grid2D<double>       & stack,
                             const NRLib::Grid2D<double> & offset);

    static void ConvertSeis(const std::vector<double>   & twt_vec,
                            const std::vector<double>   & twt_0,
                            const std::vector<double>   & zgrid_vec,
                            const std::vector<double>   & z_0,
                            const NRLib::Grid2D<double> & seismic,
                            NRLib::Grid2D<double>       & conv_seismic);

    static void NMOCorrect(const std::vector<double>   & t_in,
                           const NRLib::Grid2D<double> & data_in,
                           const NRLib::Grid2D<double> & t_out,
                           NRLib::Grid2D<double>       & data_out,
                           const std::vector<size_t>   & n_min,
                           const std::vector<size_t>   & n_max);

    static void FindNMOTheta(NRLib::Grid2D<double>     & thetagrid,
                             const std::vector<double> & twt_vec,
                             std::vector<double>       & vrms_vec,
                             const std::vector<double> & offset);

    static void ResampleTwtPS(std::vector<double>       & twt_pp_reg,
                              std::vector<double>       & twt_ss_reg,
                              const std::vector<double> & twt_pp,
                              const std::vector<double> & twt_ss,
                              const std::vector<double> & twt,
                              const std::vector<double> & twt_0,
                              double                      twt_pp_above,
                              double                      twt_ss_above,
                              double                      twt_above,
                              double                      twt_pp_below,
                              double                      twt_ss_below,
                              double                      twt_below);

    static void ResampleOffsetPS(const std::vector<double>   & twt,
                                 const NRLib::Grid2D<double> & offset_pp,
                                 const NRLib::Grid2D<double> & offset_pp_above,
                                 const NRLib::Grid2D<double> & offset_pp_below,
                                 const std::vector<double>   & twt_0,
                                 const std::vector<double>   & offset_tot_vec,
                                 NRLib::Grid2D<double>       & offset_pp_reg,
                                 NRLib::Grid2D<double>       & offset_ss_reg,
                                 double                        twt_above,
                                 double                        twt_below);

    static void ResampleTWTx(const NRLib::Grid2D<double> & twtx_grid,
                             const NRLib::Grid2D<double> & twtx_below,
                             NRLib::Grid2D<double>       & twtx_grid_reg,
                             const std::vector<double>   & twt_vec,
                             const double                & twt_below,
                             const std::vector<double>   & twt_0);

    static void FindTWTxPS(NRLib::Grid2D<double>       & twtx_grid,
                           const std::vector<double>   & twt_ss_vec,
                           const std::vector<double>   & twt_pp_vec,
                           const std::vector<double>   & vrms_pp_vec,
                           const std::vector<double>   & vrms_ss_vec,
                           const NRLib::Grid2D<double> & offset_ss,
                           const NRLib::Grid2D<double> & offset_pp,
                           bool                          offset_without_stretch);

    static void FindTWTx(NRLib::Grid2D<double>     & twtx_grid,
                         const std::vector<double> & twt_vec,
                         const std::vector<double> & vrms_vec,
                         const std::vector<double> & offset,
                         bool                        offset_without_stretch);

    static void FindSeisLimits(const NRLib::Grid2D<double> & twtx_grid,
                               const std::vector<double>   & twt_0,
                               std::vector<size_t>         & n_min,
                               std::vector<size_t>         & n_max,
                               double                        twt_wave);

    static void ExtrapolTwtAndTimeShift(std::vector<double>        & timeshift_extrapol,
                                        std::vector<double>        & twt_extrapol,
                                        const NRLib::StormContGrid & twt_timeshift,
                                        const std::vector<double>  & twt,
                                        const double                 twt_wavelet_extrapol,
                                        const size_t                 i,
                                        const size_t                 j);

    static void ExtrapolZandTwtVec(std::vector<double>        & zgrid_vec_extrapol,
                                   std::vector<double>        & twt_vec_extrapol,
                                   const std::vector<double>  & twt_vec,
                                   const NRLib::StormContGrid & zgrid,
                                   double                       vp_bot,
                                   double                       vs_bot,
                                   double                       z_wavelet_bot,
                                   size_t                       i,
                                   size_t                       j,
                                   bool                         ps_seis);

    static void MonitorInitialize(size_t   n_traces,
                                  float  & monitor_size,
                                  float  & next_monitor);

    static void Monitor(size_t   trace,
                        float    monitor_size,
                        float  & next_monitor);

    static void PrintSeisType(bool                  nmo,
                              bool                  ps_seis,
                              std::vector<double> & off_theta_vec,
                              bool                  offset_without_stretch);

    static void PrintTime();

    static std::vector<double> LinInterp1D(const std::vector<double> & x_in,
                                           const std::vector<double> & y_in,
                                           const std::vector<double> & x_out);

    static std::vector<double> SplineInterp1D(const std::vector<double> & x_in,
                                              const std::vector<double> & y_in,
                                              const std::vector<double> & x_out,
                                              double                      extrap_value);

};

#endif
