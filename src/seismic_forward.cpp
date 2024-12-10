
#include "nrlib/geometry/interpolation.hpp"

#include "nrlib/math/constants.hpp"

#include "utils/output.hpp"

#include "seismic_geometry.hpp"
#include "seismic_output.hpp"
#include "seismic_forward.hpp"
#include "wavelet.hpp"

#include <thread>
#include <ctime>
#include <map>

//---------------------------------------------------------------------------
void SeismicForward::DoSeismicForward(SeismicParameters & seismic_parameters,
                                      ModelSettings     * model_settings)
//---------------------------------------------------------------------------
{
  if (model_settings->GetNMOCorr()) {
    MakeNMOSeismic(seismic_parameters,
                   model_settings);
  }
  else {
    MakeSeismic(seismic_parameters,
                model_settings);
  }
}

//----------------------------------------------------------------------
void SeismicForward::MakeSeismic(SeismicParameters & seismic_parameters,
                                 ModelSettings     * model_settings)
//----------------------------------------------------------------------
{
  bool                  ps_seis   = model_settings->GetPSSeismic();
  std::vector<double> & theta_vec = seismic_parameters.GetThetaVec();
  std::vector<double>   twt_0;
  std::vector<double>   z_0;
  std::vector<double>   twts_0;
  size_t                dummy;

  seismic_parameters.GenerateTwt0AndZ0(model_settings,
                                       twt_0,
                                       z_0,
                                       twts_0,
                                       dummy,
                                       ps_seis);

  Output output(seismic_parameters,
                model_settings,
                twt_0,
                z_0,
                twts_0,
                theta_vec,
                twt_0.size());

  time_t t1 = time(0);

  size_t                n_traces;
  std::vector<Trace*>   seismic_traces = seismic_parameters.FindTracesInForward(n_traces);
  size_t                nzrefl         = seismic_parameters.GetSeismicGeometry()->zreflectorcount();
  std::vector<double>   dummy_vec;
  float                 monitor_size;
  float                 next_monitor;

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d traces to be generated.\n", n_traces);

  PrintSeisType(false, ps_seis, theta_vec, false);

  MonitorInitialize(n_traces, monitor_size, next_monitor);
  for (size_t k = 0; k < n_traces; ++k) {

    Trace       * trace = seismic_traces[k];

    ResultTrace result_trace(model_settings,
                             trace,
                             nzrefl,
                             twt_0.size(),
                             z_0.size(),
                             twts_0.size(),
                             twt_0.size(),
                             theta_vec.size());

    size_t i = trace->GetI();
    size_t j = trace->GetJ();

    if (GenerateTraceOk(seismic_parameters, model_settings, i, j)) {
      GenerateSeismicTraces(seismic_parameters,
                            model_settings,
                            twt_0,
                            z_0,
                            twts_0,
                            theta_vec,
                            dummy_vec,
                            &output,
                            trace,
                            &result_trace);
    }
    else {
      result_trace.SetIsEmpty(true);
    }

    output.AddTrace(&result_trace,
                    model_settings,
                    seismic_parameters.GetSeismicOutput());

    Monitor(i, monitor_size, next_monitor);
    delete trace;
  }

  std::cout << "\n";
  seismic_parameters.PrintElapsedTime(t1, "generating seismic");

  //write storm grid if requested
  output.WriteStatisticsForSeismic(model_settings);
  output.WriteSeismicStorm(model_settings,
                           seismic_parameters.GetSeismicOutput(),
                           seismic_parameters.GetRGrids());

  seismic_parameters.DeleteZandRandTWTGrids();
  seismic_parameters.DeleteElasticParameterGrids();
  seismic_parameters.DeleteWavelet();
  seismic_parameters.DeleteGeometryAndOutput();
}

//---------------------------------------------------------------------
void SeismicForward::MakeNMOSeismic(SeismicParameters & seismic_parameters,
                                    ModelSettings     * model_settings)
//---------------------------------------------------------------------
{
  bool                  ps_seis                = model_settings->GetPSSeismic();
  bool                  offset_without_stretch = model_settings->GetOffsetWithoutStretch();
  std::vector<double> & offset_vec             = seismic_parameters.GetOffsetVec();
  size_t                time_samples_stretch   = seismic_parameters.GetSeismicGeometry()->nt();
  std::vector<double>   twt_0;
  std::vector<double>   z_0;
  std::vector<double>   twts_0;

  seismic_parameters.GenerateTwt0AndZ0(model_settings,
                                       twt_0,
                                       z_0,
                                       twts_0,
                                       time_samples_stretch,
                                       ps_seis);

  Output output(seismic_parameters,
                model_settings,
                twt_0,
                z_0,
                twts_0,
                offset_vec,
                time_samples_stretch);

  time_t t1 = time(0);

  size_t                n_traces;
  std::vector<Trace*>   seismic_traces = seismic_parameters.FindTracesInForward(n_traces);
  size_t                nzrefl         = seismic_parameters.GetSeismicGeometry()->zreflectorcount();
  std::vector<double>   dummy_vec;
  float                 monitor_size;
  float                 next_monitor;


  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d traces to be generated.\n", n_traces);

  PrintSeisType(true,
                ps_seis,
                offset_vec,
                offset_without_stretch);

  MonitorInitialize(n_traces, monitor_size, next_monitor);

  for (size_t k = 0; k < n_traces; ++k) {

    Trace * trace = seismic_traces[k];

    ResultTrace result_trace(model_settings,
                             trace,
                             nzrefl,
                             twt_0.size(),
                             z_0.size(),
                             twts_0.size(),
                             time_samples_stretch,
                             offset_vec.size());

    size_t i = trace->GetI();
    size_t j = trace->GetJ();

    if (GenerateTraceOk(seismic_parameters, model_settings, i, j)) {
      GenerateNMOSeismicTraces(seismic_parameters,
                               model_settings,
                               twt_0,
                               z_0,
                               twts_0,
                               dummy_vec,
                               offset_vec,
                               time_samples_stretch,
                               &output,
                               trace,
                               &result_trace);
    }
    else {
      result_trace.SetIsEmpty(true);
    }

    output.AddTrace(&result_trace,
                    model_settings,
                    seismic_parameters.GetSeismicOutput());

    Monitor(i, monitor_size, next_monitor);
    delete trace;
  }

  std::cout << "\n";
  seismic_parameters.PrintElapsedTime(t1, "generating seismic");

  //write storm grid if requested
  output.WriteStatisticsForSeismic(model_settings);
  output.WriteSeismicStorm(model_settings,
                           seismic_parameters.GetSeismicOutput(),
                           seismic_parameters.GetRGrids());

  seismic_parameters.DeleteZandRandTWTGrids();
  seismic_parameters.DeleteElasticParameterGrids();
  seismic_parameters.DeleteWavelet();
  seismic_parameters.DeleteGeometryAndOutput();
}

//--------------------------------------------------------------------------------
void SeismicForward::GenerateNMOSeismicTraces(SeismicParameters         & seismic_parameters,
                                              ModelSettings             * model_settings,
                                              const std::vector<double> & twt_0,
                                              const std::vector<double> & z_0,
                                              const std::vector<double> & twts_0,
                                              const std::vector<double> & dummy_vec,
                                              const std::vector<double> & offset_vec,
                                              const size_t                time_samples_stretch,
                                              Output                    * nmo_output,
                                              Trace                     * trace,
                                              ResultTrace               * result_trace)
//--------------------------------------------------------------------------------
{
  size_t i = trace->GetI();
  size_t j = trace->GetJ();

  NRLib::Grid2D<double>             & timegrid_pos                = result_trace->GetPreNMOTimeTrace();
  NRLib::Grid2D<double>             & nmo_timegrid_pos            = result_trace->GetTimeTrace();
  NRLib::Grid2D<double>             & nmo_timegrid_stack_pos      = result_trace->GetTimeStackTrace();
  NRLib::Grid2D<double>             & nmo_depthgrid_pos           = result_trace->GetDepthTrace();
  NRLib::Grid2D<double>             & nmo_depthgrid_stack_pos     = result_trace->GetDepthStackTrace();
  NRLib::Grid2D<double>             & nmo_timeshiftgrid_pos       = result_trace->GetTimeShiftTrace();
  NRLib::Grid2D<double>             & nmo_timeshiftgrid_stack_pos = result_trace->GetTimeShiftStackTrace();
  NRLib::Grid2D<double>             & twtx_reg                    = result_trace->GetTWTxReg();
  NRLib::Grid2D<double>             & theta_pos                   = result_trace->GetTheta();
  NRLib::Grid2D<double>             & refl_pos                    = result_trace->GetRefl();
  NRLib::Grid2D<double>             & twtx                        = result_trace->GetTWTx();
  NRLib::Grid2D<double>               dummygrid;

  std::vector<NRLib::StormContGrid> & rgridvec                    = seismic_parameters.GetRGrids();
  NRLib::RegularSurface<double>     & toptime                     = seismic_parameters.GetTopTime();
  NRLib::StormContGrid              & zgrid                       = seismic_parameters.GetZGrid();
  NRLib::StormContGrid              & twtgrid                     = seismic_parameters.GetTwtGrid();
  NRLib::StormContGrid              & vpgrid                      = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid              & twt_timeshift               = seismic_parameters.GetTwtShiftGrid();
  size_t                              nx                          = seismic_parameters.GetSeismicGeometry()->nx();
  size_t                              ny                          = seismic_parameters.GetSeismicGeometry()->ny();
  double                              dz                          = seismic_parameters.GetSeismicGeometry()->dz();
  double                              nt_non_nmo                  = seismic_parameters.GetSeismicGeometry()->nt();
  double                              nz_non_nmo                  = seismic_parameters.GetSeismicGeometry()->nz();
  double                              dt                          = seismic_parameters.GetSeismicGeometry()->dt();
  double                              t0                          = seismic_parameters.GetSeismicGeometry()->t0();
  double                              z0                          = seismic_parameters.GetSeismicGeometry()->z0();
  size_t                              nzrefl                      = seismic_parameters.GetSeismicGeometry()->zreflectorcount();

  double                              wavelet_scale               = seismic_parameters.GetWaveletScale();
  Wavelet                           * wavelet                     = seismic_parameters.GetWavelet();
  double                              twt_wavelet                 = wavelet->GetTwtLength();

  std::vector<double>                 constvp                     = model_settings->GetConstVp();
  std::vector<double>                 constvs                     = model_settings->GetConstVs();
  bool                                ps_seis                     = model_settings->GetPSSeismic();
  double                              z_wavelet_bot               = model_settings->GetZWaveletBot();
  double                              z_extrapol_factor           = model_settings->GetZExtrapolFactor();
  bool                                offset_without_stretch      = model_settings->GetOffsetWithoutStretch();
  double                              z_w                         = model_settings->GetZw();
  double                              v_w                         = model_settings->GetVw();
  bool                                output_refl                 = model_settings->GetOutputReflections();
  bool                                add_noise                   = model_settings->GetAddNoiseToReflCoef();
  bool                                add_white_noise             = model_settings->GetAddWhiteNoise();
  bool                                equal_noise                 = model_settings->GetUseEqualNoiseForOffsets();
  double                              sd1                         = model_settings->GetStandardDeviation1();
  double                              sd2                         = model_settings->GetStandardDeviation2();
  unsigned long                       seed1                       = model_settings->GetSeed1();
  unsigned long                       seed2                       = model_settings->GetSeed2();

  double                              tmin                        = twt_0[0];
  int                                 noff                        = offset_vec.size();
  int                                 nt                          = twt_0.size();

  std::vector<size_t> n_min(noff);
  std::vector<size_t> n_max(noff);
  std::vector<double> twt_vec(nzrefl);
  std::vector<double> vp_vec(nzrefl);

  //get twt and vp at all layers from twtgrid and vpgrid
  for (size_t k = 0; k < nzrefl; ++k) {
    twt_vec[k] = twtgrid(i, j, k);
    vp_vec[k]  = vpgrid(i, j, k);
  }

  if (ps_seis) { //------------PS seismic------------
    //setup and get vectors and grid for calculation
    std::vector<double>     twt_ss_vec(nzrefl);
    std::vector<double>     twt_pp_vec(nzrefl);
    std::vector<double>     twt_ss_vec_reg(nt);
    std::vector<double>     twt_pp_vec_reg(nt);
    std::vector<double>     vs_vec(nzrefl);

    std::vector<double>     vrms_pp_vec(nzrefl);
    std::vector<double>     vrms_pp_vec_reg(nt);
    std::vector<double>     vrms_ss_vec(nzrefl);
    std::vector<double>     vrms_ss_vec_reg(nt);

    NRLib::Grid2D<double> & offset_pp     = result_trace->GetOffsetPP();
    NRLib::Grid2D<double> & offset_ss     = result_trace->GetOffsetSS();
    NRLib::Grid2D<double> & offset_pp_reg = result_trace->GetOffsetPPReg();
    NRLib::Grid2D<double> & offset_ss_reg = result_trace->GetOffsetSSReg();

    NRLib::StormContGrid  & vsgrid        = seismic_parameters.GetVsGrid();
    NRLib::StormContGrid  & twtssgrid     = seismic_parameters.GetTwtSSGrid();
    NRLib::StormContGrid  & twtppgrid     = seismic_parameters.GetTwtPPGrid();

    //get SS and PP twt, and vs at all layers from twtgrid and vsgrid
    for (size_t k = 0; k < nzrefl; ++k) {
      twt_ss_vec[k] = twtssgrid(i, j, k);
      twt_pp_vec[k] = twtppgrid(i, j, k);
      vs_vec[k]     = vsgrid(i, j, k);
    }

    //find twtx resample all
    double twt_ss_extrapol = 2000 / constvs[2] * z_wavelet_bot * z_extrapol_factor;
    double twt_pp_extrapol = 2000 / constvp[2] * z_wavelet_bot * z_extrapol_factor;
    double twt_ps_extrapol = 0.5 * (twt_ss_extrapol + twt_pp_extrapol);

    double twt_below       = twt_vec   [nzrefl - 1] + twt_ps_extrapol;
    double twt_pp_below    = twt_pp_vec[nzrefl - 1] + twt_pp_extrapol;
    double twt_ss_below    = twt_ss_vec[nzrefl - 1] + twt_ss_extrapol;

    double twt_w           = 2000 * z_w /v_w; //twt_above pp, ps, ss

    //resample twt_pp and twt_ss
    ResampleTwtPS(twt_pp_vec_reg, twt_ss_vec_reg, twt_pp_vec, twt_ss_vec, twt_vec, twt_0, twt_w, twt_w, twt_w, twt_pp_below, twt_ss_below, twt_below);

    //find vrms for pp and ss -  for each layer and regularly sampled
    seismic_parameters.FindVrms       (vrms_pp_vec, twt_pp_vec, vp_vec, zgrid(i, j, 0));
    seismic_parameters.FindVrms       (vrms_ss_vec, twt_ss_vec, vs_vec, zgrid(i, j, 0));
    seismic_parameters.FindVrmsRegular(vrms_pp_vec, vrms_pp_vec_reg, twt_pp_vec, twt_pp_vec_reg, vp_vec, constvp[2], twt_pp_extrapol);
    seismic_parameters.FindVrmsRegular(vrms_ss_vec, vrms_ss_vec_reg, twt_ss_vec, twt_ss_vec_reg, vs_vec, constvs[2], twt_ss_extrapol);

    //find theta and offset - for each layer for each offset:
    seismic_parameters.FindPSNMOThetaAndOffset(theta_pos, offset_pp, offset_ss, twt_pp_vec, twt_ss_vec, vrms_pp_vec, vrms_ss_vec, offset_vec);

    //find offset above and below reservoir - and resample regularly
    NRLib::Grid2D<double> offset_pp_above(1, noff), offset_ss_above(1, noff);
    NRLib::Grid2D<double> offset_pp_below(1, noff), offset_ss_below(1, noff);
    // - above reservoir
    std::vector<double>   twt_pp_one (1, twt_pp_vec_reg[0] ); //twt_w
    std::vector<double>   twt_ss_one (1, twt_ss_vec_reg[0] ); //twt_w
    std::vector<double>   vrms_pp_one(1, vrms_pp_vec_reg[0]); //v_w
    std::vector<double>   vrms_ss_one(1, vrms_ss_vec_reg[0]); //v_w
    seismic_parameters.FindPSNMOThetaAndOffset(dummygrid, offset_pp_above, offset_ss_above, twt_pp_one, twt_ss_one, vrms_pp_one, vrms_ss_one, offset_vec, false);
    // - below reservoir
    twt_pp_one [0] = twt_pp_below;
    twt_ss_one [0] = twt_ss_below;
    vrms_pp_one[0] = std::sqrt(1 / twt_pp_below * (constvp[2] * constvp[2] * (twt_pp_below - twt_pp_vec[nzrefl - 1]) + vrms_pp_vec[nzrefl - 1] * vrms_pp_vec[nzrefl - 1] * twt_pp_vec[nzrefl - 1]));
    vrms_ss_one[0] = std::sqrt(1 / twt_ss_below * (constvs[2] * constvs[2] * (twt_ss_below - twt_ss_vec[nzrefl - 1]) + vrms_ss_vec[nzrefl - 1] * vrms_ss_vec[nzrefl - 1] * twt_ss_vec[nzrefl - 1]));
    seismic_parameters.FindPSNMOThetaAndOffset(dummygrid, offset_pp_below, offset_ss_below, twt_pp_one, twt_ss_one, vrms_pp_one, vrms_ss_one, offset_vec, false);
    // - resample regularly
    ResampleOffsetPS(twt_vec, offset_pp, offset_pp_above, offset_pp_below, twt_0, offset_vec, offset_pp_reg, offset_ss_reg, tmin, twt_below);

    //find twtx - for each layer for each offset - and regularly sampled
    FindTWTxPS(twtx,     twt_ss_vec,     twt_pp_vec,     vrms_pp_vec,     vrms_ss_vec,     offset_ss,     offset_pp,     offset_without_stretch);
    FindTWTxPS(twtx_reg, twt_ss_vec_reg, twt_pp_vec_reg, vrms_pp_vec_reg, vrms_ss_vec_reg, offset_ss_reg, offset_pp_reg, offset_without_stretch);

    //find limits for where to generate seismic, for each offset
    FindSeisLimits(twtx, twt_0, n_min, n_max, twt_wavelet);

  }
  else { //------------PP seismic------------
    std::vector<double> vrms_vec(nzrefl);
    std::vector<double> vrms_vec_reg(twt_0);
    double twt_wavelet_extrapol = twt_wavelet * z_extrapol_factor;
    seismic_parameters.FindVrms       (vrms_vec, twt_vec, vp_vec, zgrid(i, j, 0));
    seismic_parameters.FindVrmsRegular(vrms_vec, vrms_vec_reg, twt_vec, twt_0, vp_vec, constvp[2], twt_wavelet_extrapol);

    FindNMOTheta  (theta_pos, twt_vec, vrms_vec    , offset_vec);                         // Find theta - for each layer for each offset:
    FindTWTx      (twtx     , twt_vec, vrms_vec    , offset_vec, offset_without_stretch); // Find twtx for each layer for each offset, and regularly in time:
    FindTWTx      (twtx_reg , twt_0  , vrms_vec_reg, offset_vec, offset_without_stretch);
    FindSeisLimits(twtx     , twt_0  , n_min       , n_max     , twt_wavelet);            // Find limits for where to generate seismic, for each offset
  }      //----------------------------------

  MakeReflections(refl_pos,             // Also add noise if requested
                  rgridvec,
                  seismic_parameters,
                  theta_pos,
                  output_refl,
                  add_noise,
                  sd2,
                  seed2,
                  nx,
                  i,
                  j);

  SeisConvolutionNMO(timegrid_pos,      // Generate seismic
                     refl_pos,
                     twtx,
                     zgrid,
                     toptime,
                     wavelet,
                     wavelet_scale,
                     offset_vec,
                     tmin,
                     dt,
                     i,
                     j,
                     n_min,
                     n_max);
  //
  // NMO correction:
  //
  if (offset_without_stretch) {
    nmo_timegrid_pos = NRLib::Grid2D<double>(timegrid_pos);
  }
  else {
    NMOCorrect(twt_0,
               timegrid_pos,
               twtx_reg,
               nmo_timegrid_pos, // output
               n_min,
               n_max);
  }

  bool depth_conversion = nmo_output->GetDepthSegyOk()          || nmo_output->GetDepthStackSegyOk()     || model_settings->GetDepthStormOutput();
  bool time_shift       = nmo_output->GetTimeshiftSegyOk()      || nmo_output->GetTimeshiftStackSegyOk() || model_settings->GetTimeshiftStormOutput();
  bool stack_time       = model_settings->GetStackOutput()      || model_settings->GetStormOutput();
  bool stack_timeshift  = nmo_output->GetTimeshiftStackSegyOk() || model_settings->GetTimeshiftStormOutput();
  bool stack_depth      = nmo_output->GetDepthStackSegyOk()     || model_settings->GetDepthStormOutput();

  int  tshift           = static_cast<int>(floor((t0 - twt_0[0]) / dt + 0.5));   //| Align noise with zero offset
  int  zshift           = static_cast<int>(floor((z0 - z_0[0]  ) / dz + 0.5));   //|

  ConvertShiftAndStack(nmo_timegrid_stack_pos,
                       nmo_timeshiftgrid_stack_pos,
                       nmo_depthgrid_stack_pos,
                       nmo_depthgrid_pos,
                       nmo_timeshiftgrid_pos,
                       nmo_timegrid_pos,
                       twt_timeshift,
                       twt_vec,
                       twts_0,
                       twt_0,
                       zgrid,
                       z_0,
                       constvp[2],
                       constvs[2],
                       z_wavelet_bot*z_extrapol_factor,
                       twt_wavelet*z_extrapol_factor,
                       sd1,
                       tshift,
                       zshift,
                       depth_conversion,
                       time_shift,
                       add_white_noise,
                       equal_noise,
                       ps_seis,
                       stack_time,
                       stack_timeshift,
                       stack_depth,
                       nzrefl,
                       noff,
                       seed1 + static_cast<long>(i + nx*j),
                       nx*ny, // Number of traces
                       nt_non_nmo,
                       nz_non_nmo,
                       i,
                       j);

}

//---------------------------------------------------------------------
void SeismicForward::GenerateSeismicTraces(SeismicParameters         & seismic_parameters,
                                           ModelSettings             * model_settings,
                                           const std::vector<double> & twt_0,
                                           const std::vector<double> & z_0,
                                           const std::vector<double> & twts_0,
                                           const std::vector<double> & theta_vec,
                                           const std::vector<double> & dummy_vec,
                                           Output                    * output,
                                           Trace                     * trace,
                                           ResultTrace               * result_trace)
//---------------------------------------------------------------------
{
  size_t i = trace->GetI();
  size_t j = trace->GetJ();

  NRLib::Grid2D<double>             & timegrid_pos            = result_trace->GetTimeTrace();
  NRLib::Grid2D<double>             & timegrid_stack_pos      = result_trace->GetTimeStackTrace();
  NRLib::Grid2D<double>             & depthgrid_pos           = result_trace->GetDepthTrace();
  NRLib::Grid2D<double>             & depthgrid_stack_pos     = result_trace->GetDepthStackTrace();
  NRLib::Grid2D<double>             & timeshiftgrid_pos       = result_trace->GetTimeShiftTrace();
  NRLib::Grid2D<double>             & timeshiftgrid_stack_pos = result_trace->GetTimeShiftStackTrace();

  std::vector<NRLib::StormContGrid> & rgridvec                = seismic_parameters.GetRGrids();
  NRLib::RegularSurface<double>     & toptime                 = seismic_parameters.GetTopTime();
  NRLib::StormContGrid              & zgrid                   = seismic_parameters.GetZGrid();
  NRLib::StormContGrid              & twtgrid                 = seismic_parameters.GetTwtGrid();
  NRLib::StormContGrid              & twt_timeshift           = seismic_parameters.GetTwtShiftGrid();

  size_t                              nx                      = seismic_parameters.GetSeismicGeometry()->nx();
  size_t                              ny                      = seismic_parameters.GetSeismicGeometry()->ny();
  double                              dz                      = seismic_parameters.GetSeismicGeometry()->dz();
  size_t                              nt                      = seismic_parameters.GetSeismicGeometry()->nt();
  size_t                              nz                      = seismic_parameters.GetSeismicGeometry()->nz();
  double                              dt                      = seismic_parameters.GetSeismicGeometry()->dt();
  double                              t0                      = seismic_parameters.GetSeismicGeometry()->t0();
  double                              z0                      = seismic_parameters.GetSeismicGeometry()->z0();
  size_t                              nzrefl                  = seismic_parameters.GetSeismicGeometry()->zreflectorcount();
  double                              wavelet_scale           = seismic_parameters.GetWaveletScale();
  Wavelet                           * wavelet                 = seismic_parameters.GetWavelet();
  double                              twt_wavelet             = wavelet->GetTwtLength();

  NRLib::Grid2D<double>               refl_pos(nzrefl, theta_vec.size());
  std::vector<double>                 twt_vec(nzrefl);

  std::vector<double>                 constvp                 = model_settings->GetConstVp();
  std::vector<double>                 constvs                 = model_settings->GetConstVs();
  bool                                ps_seis                 = model_settings->GetPSSeismic();
  double                              z_wavelet_bot           = model_settings->GetZWaveletBot();
  bool                                output_refl             = model_settings->GetOutputReflections();
  bool                                add_noise               = model_settings->GetAddNoiseToReflCoef();
  bool                                add_white_noise         = model_settings->GetAddWhiteNoise();
  bool                                equal_noise             = model_settings->GetUseEqualNoiseForOffsets();
  double                              sd1                     = model_settings->GetStandardDeviation1();
  double                              sd2                     = model_settings->GetStandardDeviation2();
  unsigned long                       seed1                   = model_settings->GetSeed1();
  unsigned long                       seed2                   = model_settings->GetSeed2();

  double                              tmin                    = twt_0[0];
  int                                 ntheta                  = theta_vec.size();

  size_t                              n_min                   = 0;
  size_t                              n_max                   = nt;

  //get twt at all layers from twtgrid
  for (size_t k = 0; k < nzrefl; ++k) {
    twt_vec[k]   = twtgrid(i,j,k);
  }

  //size_t kdim = seismic_parameters.GetBottomK() - seismic_parameters.GetTopK() + 3;
  NRLib::Grid2D<double> theta(nzrefl, ntheta);
  for (size_t k = 0 ; k < nzrefl ; ++k) {
    for (size_t t = 0 ; t < ntheta ; ++t) {
      theta(k, t) = theta_vec[t];
    }
  }

  MakeReflections(refl_pos,             // Also add noise if requested
                  rgridvec,
                  seismic_parameters,
                  theta,
                  output_refl,
                  add_noise,
                  sd2,
                  seed2,
                  nx,
                  i,
                  j);

  SeisConvolution(timegrid_pos,         // Generate seismic
                  refl_pos,
                  twt_vec,
                  zgrid,
                  toptime,
                  wavelet,
                  wavelet_scale,
                  theta_vec,
                  tmin,
                  dt,
                  i,
                  j,
                  n_min,
                  n_max);

  bool depth_conversion = output->GetDepthSegyOk()          || output->GetDepthStackSegyOk()     || model_settings->GetDepthStormOutput();
  bool time_shift       = output->GetTimeshiftSegyOk()      || output->GetTimeshiftStackSegyOk() || model_settings->GetTimeshiftStormOutput();
  bool stack_time       = model_settings->GetStackOutput()  || model_settings->GetStormOutput();
  bool stack_timeshift  = output->GetTimeshiftStackSegyOk() || model_settings->GetTimeshiftStormOutput();
  bool stack_depth      = output->GetDepthStackSegyOk()     || model_settings->GetDepthStormOutput();

  int  tshift           = static_cast<int>(floor((t0 - twt_0[0]) / dt + 0.5));   //| Align noise with zero offset
  int  zshift           = static_cast<int>(floor((z0 - z_0[0]  ) / dz + 0.5));   //|

  ConvertShiftAndStack(timegrid_stack_pos,
                       timeshiftgrid_stack_pos,
                       depthgrid_stack_pos,
                       depthgrid_pos,
                       timeshiftgrid_pos,
                       timegrid_pos,
                       twt_timeshift,
                       twt_vec,
                       twts_0,
                       twt_0,
                       zgrid,
                       z_0,
                       constvp[2],
                       constvs[2],
                       z_wavelet_bot,
                       twt_wavelet,
                       sd1,
                       tshift,
                       zshift,
                       depth_conversion,
                       time_shift,
                       add_white_noise,
                       equal_noise,
                       ps_seis,
                       stack_time,
                       stack_timeshift,
                       stack_depth,
                       nzrefl,
                       ntheta,
                       seed1 + static_cast<long>(i + nx*j),
                       nx*ny, // Number of traces
                       nt,
                       nz,
                       i,
                       j);

}

//-----------------------------------------------------------------------------------
void SeismicForward::MakeReflections(NRLib::Grid2D<double>             & refl,
                                     std::vector<NRLib::StormContGrid> & rgridvec,
                                     SeismicParameters                 & seismic_parameters,
                                     NRLib::Grid2D<double>             & theta,
                                     bool                                output_refl,
                                     bool                                add_noise,
                                     double                              std,
                                     unsigned long                       seed,
                                     size_t                              nx,
                                     size_t                              i,
                                     size_t                              j)
//-----------------------------------------------------------------------------------
{
  seismic_parameters.FindReflections(refl, theta, i, j);

  size_t m = refl.GetNI();
  size_t n = refl.GetNJ();

  if (output_refl) { // Keep reflections for zero offset if output on storm
    for (size_t k = 0 ; k < m ; ++k) {
      rgridvec[0](i, j, k) = static_cast<float>(refl(k, 0));
    }
  }
  if (add_noise) {
    std::vector<double> noise(m);
    for (int jj = 0 ; jj < n ; jj++) {
      noise = GenerateWhiteNoise(seed + static_cast<long>(i + nx*j), std, m); // Gives equal noise for equal angles
      for (int ii = 0 ; ii < m ; ii++) {
        refl(ii, jj) += noise[ii];
      }
    }
    if (output_refl) { // Keep reflections for zero offset if output on storm and white noise
      for (size_t k = 0 ; k < m ; ++k) {
        rgridvec[1](i, j, k) = static_cast<float>(refl(k, 0));
      }
    }
  }
}

//--------------------------------------------------------------------------------------------
void SeismicForward::ConvertShiftAndStack(NRLib::Grid2D<double>      & timegrid_stack,
                                          NRLib::Grid2D<double>      & timeshiftgrid_stack,
                                          NRLib::Grid2D<double>      & depthgrid_stack,
                                          NRLib::Grid2D<double>      & depthgrid,
                                          NRLib::Grid2D<double>      & timeshiftgrid,
                                          NRLib::Grid2D<double>      & timegrid,
                                          const NRLib::StormContGrid & twt_timeshift,
                                          const std::vector<double>  & twt_vec,
                                          const std::vector<double>  & twts_0,
                                          const std::vector<double>  & twt_0,
                                          const NRLib::StormContGrid & zgrid,
                                          const std::vector<double>  & z_0,
                                          const double                 vp_bot,
                                          const double                 vs_bot,
                                          const double                 z_wavelet_bot_stretched, // NMO => z_wavelet_bot * z_extrapol_factor
                                          const double                 twt_wavelet_stretched,   // NMO => twt_wavelet   * z_extrapol_factor
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
                                          const size_t                 j)
//--------------------------------------------------------------------------------------
{
  //
  // Depth conversion:
  //
  if (depth_conversion) {
    std::vector<double> zgrid_vec_extrapol(nzrefl + 2);
    std::vector<double> twt_vec_extrapol  (nzrefl + 2);
    ExtrapolZandTwtVec(zgrid_vec_extrapol,
                       twt_vec_extrapol,
                       twt_vec,
                       zgrid,
                       vp_bot,
                       vs_bot,
                       z_wavelet_bot_stretched,
                       i,
                       j,
                       ps_seis);
    ConvertSeis(twt_vec_extrapol,
                twt_0,
                zgrid_vec_extrapol,
                z_0,
                timegrid,
                depthgrid);
  }

  //
  // Time shift:
  //
  if (time_shift) {
    std::vector<double> timeshiftgrid_vec_extrapol(nzrefl + 2);
    std::vector<double> twt_vec_extrapol          (nzrefl + 2);
    ExtrapolTwtAndTimeShift(timeshiftgrid_vec_extrapol,
                            twt_vec_extrapol,
                            twt_timeshift,
                            twt_vec,
                            twt_wavelet_stretched,
                            i,
                            j);
    ConvertSeis(twt_vec_extrapol,
                twt_0,
                timeshiftgrid_vec_extrapol,
                twts_0,
                timegrid,
                timeshiftgrid);
  }

  //
  // Add white noise to seismic signal
  //
  if (add_white_noise) {
    std::vector<std::vector<double> > noise_time(noff);
    GenerateWhiteNoiseAllOffsets(noise_time, equal_noise, seed, sd, n_traces, nt_non_nmo);

    std::vector<std::vector<double> > noise_depth(noff);
    GenerateWhiteNoiseAllOffsets(noise_depth, equal_noise, seed, sd, n_traces, nz_non_nmo);

    for (int off = 0 ; off < noff ; off++) {
      for (int ii = 0 ; ii < nt_non_nmo ; ii++) {
        timegrid(tshift + ii, off) += noise_time[off][ii];           // Add noise to time
      }
      if (time_shift) {
        for (int ii = 0 ; ii < nt_non_nmo ; ii++) {
          timeshiftgrid(tshift + ii, off) += noise_time[off][ii];    // Add noise to time shift. Is tshift correct?
        }
      }
      if (depth_conversion) {
        for (int ii = 0 ; ii < nz_non_nmo ; ii++) {
          depthgrid(zshift + ii, off) += noise_depth[off][ii];       // Add noise to depth
        }
      }
    }
  }

  //
  // Stack time, time shift and depth
  //
  if (stack_time) {
    StackOffsets(timegrid_stack, timegrid);
  }
  if (stack_timeshift) {
    StackOffsets(timeshiftgrid_stack, timeshiftgrid);
  }
  if (stack_depth) {
    StackOffsets(depthgrid_stack, depthgrid);
  }
}

//---------------------------------------------------------------------
void SeismicForward::StackOffsets(NRLib::Grid2D<double>       & stack,
                                  const NRLib::Grid2D<double> & offset)
//---------------------------------------------------------------------
{
  size_t nk   = offset.GetNI();
  size_t noff = offset.GetNJ();
  float scale = static_cast<float>(1.0 / noff);
  for (size_t k = 0 ; k < nk ; ++k) {
    stack(k,0) = 0.0;
    for (size_t off = 0; off < noff; ++off) {
      stack(k,0) += scale * offset(k,off);
    }
  }
}

//------------------------------------------------------------------------------------------------
void SeismicForward::GenerateWhiteNoiseAllOffsets(std::vector<std::vector<double> > & noise,
                                                  bool                                equal_noise,
                                                  unsigned long                       seed,
                                                  double                              sd,
                                                  size_t                              n_traces,
                                                  size_t                              nt)
//------------------------------------------------------------------------------------------------
{
  noise[0] = GenerateWhiteNoise(seed, sd, nt);
  for (int off = 1 ; off < noise.size() ; off++) {
    if (equal_noise)
      noise[off] = noise[0];                                        // Equal noise for each offset
    else
      noise[off] = GenerateWhiteNoise(seed + n_traces*off, sd, nt); // Unique noise for each offset
  }
}

//-------------------------------------------------------------------------------
std::vector<double> SeismicForward::GenerateWhiteNoise(unsigned long seed,
                                                       double        sd,
                                                       size_t        n)
//-------------------------------------------------------------------------------
{
  std::vector<double> noise(n);
  NRLib::RandomGenerator rg;
  rg.Initialize(seed);
  NRLib::Normal normal(0, sd);

  for (size_t i = 0 ; i < n ; ++i) {
    noise[i] = static_cast<float>(normal.Draw(rg));
  }
  return noise;
}

//---------------------------------------------------------
void SeismicForward::MonitorInitialize(size_t n_traces,
                                       float &monitor_size,
                                       float &next_monitor)
//---------------------------------------------------------
{
  monitor_size = static_cast<float>(n_traces) * 0.02f;
  if (monitor_size < 1.0f)
    monitor_size = 1.0f;
  next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
}

//-----------------------------------------------
void SeismicForward::Monitor(size_t trace,
                             float monitor_size,
                             float &next_monitor)
//-----------------------------------------------
{
  if (trace +1 >= static_cast<size_t>(next_monitor)) {
    next_monitor += monitor_size;
    std::cout << "^";
    fflush(stdout);
    if (next_monitor > monitor_size*51){
      std::cout << "\n";
    }
  }
}

//------------------------------------------------------------------------------
void SeismicForward::PrintSeisType(bool                  nmo,
                                   bool                  ps_seis,
                                   std::vector<double> & off_theta_vec,
                                   bool                  offset_without_stretch)
//------------------------------------------------------------------------------
{
  if (nmo) {
    if (offset_without_stretch) {
      if (ps_seis) NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PS-seismic for offsets: ");
      else         NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PP-seismic for offsets: ");
    }
    else {
      if (ps_seis) NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PS-seismic with stretch for offsets: ");
      else         NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PP-seismic with stretch for offsets: ");
    }
    for (size_t i = 0; i < off_theta_vec.size(); ++i) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "%.1f ", off_theta_vec[i]);
    }
  }
  else {
    if (ps_seis)   NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic PS-seismic for angles: ");
    else           NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic PP-seismic for angles: ");
    for (size_t i = 0; i < off_theta_vec.size(); ++i) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "%.1f ", (off_theta_vec[i]/NRLib::Degree));
    }
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n");
}

void SeismicForward::PrintTime()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  std::cout << "Time: "
    << (now->tm_hour) << ':'
    << (now->tm_min) << ':'
    <<  now->tm_sec
    << "\n";
}

//--------------------------------------------------------------------------
bool SeismicForward::GenerateTraceOk(SeismicParameters & seismic_parameters,
                                     ModelSettings     * model_settings,
                                     size_t              i,
                                     size_t              j)
//--------------------------------------------------------------------------
{
  bool   generate_ok = false;
  double const_vp    = model_settings->GetConstVp()[1];
  double const_vs    = model_settings->GetConstVs()[1];
  double const_rho   = model_settings->GetConstRho()[1];

  NRLib::StormContGrid & vpgrid    = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid & vsgrid    = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid & rhogrid   = seismic_parameters.GetRhoGrid();
  NRLib::StormContGrid & twtgrid   = seismic_parameters.GetTwtGrid();

  if (twtgrid(i, j, 0) != -999) {
    size_t nk = vpgrid.GetNK();

    for (size_t k = 1; k < nk - 1 ; k++) {
      if (generate_ok == false) {
        if (vpgrid(i, j, k) != const_vp){
          generate_ok = true;
        }
        if (vsgrid(i, j, k) != const_vs){
          generate_ok = true;
        }
        if (rhogrid(i, j, k) != const_rho){
          generate_ok = true;
        }
      }
      else{
        break;
      }
    }
  }
  return generate_ok;
}

std::vector<double> SeismicForward::LinInterp1D(const std::vector<double> &x_in,
                                                const std::vector<double> &y_in,
                                                const std::vector<double> &x_out)

{
  std::vector<double> x_in_copy(x_in.size());
  std::vector<double> y_in_copy(y_in.size());
  x_in_copy[0] = x_in[0];
  y_in_copy[0] = y_in[0];
  size_t index = 1; //sjekk om denne skal være 0???
  for (size_t i = 1; i < x_in.size(); ++i){
    if (x_in[i] != x_in[i-1]){
      x_in_copy[index] = x_in[i];
      y_in_copy[index] = y_in[i];
      ++index;
    }
  }
  x_in_copy.resize(index);
  y_in_copy.resize(index);
  return NRLib::Interpolation::Interpolate1D(x_in_copy, y_in_copy, x_out, "linear");
}

std::vector<double> SeismicForward::SplineInterp1D(const std::vector<double> &x_in,
                                                   const std::vector<double> &y_in,
                                                   const std::vector<double> &x_out,
                                                   double                     extrap_value)

{
  std::vector<double> x_in_copy(x_in.size());
  std::vector<double> y_in_copy(y_in.size());
  x_in_copy[0] = x_in[0];
  y_in_copy[0] = y_in[0];
  size_t index = 1;
  for (size_t i = 1; i < x_in.size(); ++i){
    if (x_in[i] != x_in[i-1]){
      x_in_copy[index] = x_in[i];
      y_in_copy[index] = y_in[i];
      ++index;
    }
  }
  x_in_copy.resize(index);
  y_in_copy.resize(index);
  return NRLib::Interpolation::Interpolate1D(x_in_copy, y_in_copy, x_out, "spline", extrap_value);
}

//---------------------------------------------------------------------------------------------
void SeismicForward::ExtrapolTwtAndTimeShift(std::vector<double>        & timeshift_extrapol,
                                             std::vector<double>        & twt_extrapol,
                                             const NRLib::StormContGrid & twt_timeshift,
                                             const std::vector<double>  & twt,
                                             const double                 twt_wavelet_extrapol,
                                             const size_t                 i,
                                             const size_t                 j)
//---------------------------------------------------------------------------------------------
{
  size_t nzrefl = twt_extrapol.size() - 2;

  timeshift_extrapol[0] = 0.0;
  twt_extrapol      [0] = 0.0;
  for (size_t k = 1; k < nzrefl + 1 ; ++k) {
    timeshift_extrapol[k] = twt_timeshift(i, j, k - 1);
    twt_extrapol[k]       = twt[k - 1];
  }
  timeshift_extrapol[nzrefl + 1] = twt_timeshift(i,j,nzrefl - 1) + twt_wavelet_extrapol *twt_timeshift(i, j, nzrefl - 1)/ twt[nzrefl - 1];
  twt_extrapol      [nzrefl + 1] = twt[nzrefl - 1]               + twt_wavelet_extrapol;
}

//----------------------------------------------------------------------------------------
void SeismicForward::ExtrapolZandTwtVec(std::vector<double>        & zgrid_vec_extrapol,
                                        std::vector<double>        & twt_vec_extrapol,
                                        const std::vector<double>  & twt_vec,
                                        const NRLib::StormContGrid & zgrid,
                                        double                       vp_bot,
                                        double                       vs_bot,
                                        double                       z_wavelet_bot,
                                        size_t                       i,
                                        size_t                       j,
                                        bool                         ps_seis)
//----------------------------------------------------------------------------------------
{
  double vel_bot;
  if (ps_seis)
    vel_bot = 2 / (1 / vp_bot + 1 / vs_bot);
  else
    vel_bot = vp_bot;
  size_t nzrefl = zgrid_vec_extrapol.size() - 2;
  zgrid_vec_extrapol[0] = 0.0;
  twt_vec_extrapol  [0] = 0.0;
  for (size_t k = 0; k < nzrefl; ++k) {
    twt_vec_extrapol  [k+1] = twt_vec[k];
    zgrid_vec_extrapol[k+1] = zgrid(i, j, k);
  }
  zgrid_vec_extrapol[nzrefl+1] = zgrid_vec_extrapol[nzrefl] + z_wavelet_bot;
  twt_vec_extrapol  [nzrefl+1] = twt_vec_extrapol[nzrefl]   + 2000 * z_wavelet_bot / vel_bot;
}

void SeismicForward::ConvertSeis(const std::vector<double>   & twt_vec,
                                 const std::vector<double>   & twt_0,
                                 const std::vector<double>   & zgrid_vec,
                                 const std::vector<double>   & z_0,
                                 const NRLib::Grid2D<double> & seismic,
                                 NRLib::Grid2D<double>       & conv_seismic)
{
  size_t ni   = seismic.GetNI();
  size_t noff = seismic.GetNJ();
  size_t nk   = conv_seismic.GetNI();

  std::vector<double> seismic_vec(ni);
  std::vector<double> conv_seismic_vec(nk);

  std::vector<double> zt_reg  // nt               (y_out)
    = LinInterp1D(twt_vec,    // nzrefl + 1 or 2  (x_in )
                  zgrid_vec,  // nzrefl + 1 or 2  (y_in )
                  twt_0);     // nt               (x_out)

  zt_reg.resize(seismic.GetNI());

  for (size_t off = 0 ; off < noff ; off++) {
    for (size_t k = 0 ; k < ni ; k++) {
      seismic_vec[k] = seismic(k, off);
    }

    conv_seismic_vec                 // nk   | y_out
      = SplineInterp1D(zt_reg,       // ni   | x_in
                       seismic_vec,  // ni   | y_in
                       z_0,          // nk   | x_out
                       0.0);         //      | Extrapolation value

    for (size_t k = 0; k < nk; k++) {
      conv_seismic(k, off) = conv_seismic_vec[k];
    }
  }
}

void SeismicForward::NMOCorrect(const std::vector<double>   & t_in,
                                const NRLib::Grid2D<double> & data_in,
                                const NRLib::Grid2D<double> & t_out,
                                NRLib::Grid2D<double>       & data_out,
                                const std::vector<size_t>   & n_min,
                                const std::vector<size_t>   & n_max)
{
  size_t nt_in = data_in.GetNI();
  size_t noff  = data_in.GetNJ();
  std::vector<double> data_vec_in(nt_in);
  std::vector<double> data_vec_out(nt_in);
  std::vector<double> t_vec_in(nt_in);
  std::vector<double> t_vec_out(nt_in);

  for (size_t off = 0; off < noff; off++) {
    size_t n_min_max = (n_max[off]-n_min[off]+1);
    data_vec_in.resize(n_min_max);
    t_vec_in.resize(n_min_max);
    size_t index = 0;
    t_vec_out.resize(nt_in);

    //only interpolate FROM within min-max
    for (size_t k = n_min[off]; k <= n_max[off]; ++k) {
      data_vec_in[k - n_min[off]] = data_in(k,off);
      t_vec_in[k - n_min[off]]    = t_in[k];
    }

    // not necessary to interpolate AT values higher than max t
    // t_out not monotonously increasing, must check that we are inside
    bool inside = false;
    for (size_t k = 0; k < nt_in; k++) {
      t_vec_out[k]   = t_out(k, off);
      if (inside == false && (t_vec_out[k] > t_vec_in[0] && t_vec_out[k] < t_vec_in[n_min_max - 1])){
        inside = true;
      }
      if (inside == true && t_vec_out[k] > t_vec_in[n_min_max-1]){
        break;
      }
      ++index;
    }
    t_vec_out.resize(index);

    //interpolate
    data_vec_out = SplineInterp1D(t_vec_in, data_vec_in, t_vec_out, 0);
    if (index > data_out.GetNI()){
      std::cout << "DEBUG ERROR: stretch not properly accounted for.\n";
      index = data_out.GetNI();
    }

    for (size_t k = 0; k < index; k++) {
      data_out(k, off) = data_vec_out[k];
    }
    //fill in zeros at higher values than max t
    for (size_t k = index; k < data_out.GetNI(); k++){
      data_out(k, off) = 0.0;
    }
  }
}

void SeismicForward::FindNMOTheta(NRLib::Grid2D<double>       &thetagrid,
                                  const std::vector<double>   &twt_vec,
                                  std::vector<double>         &vrms_vec,
                                  const std::vector<double>   &offset){

  for (size_t off = 0; off < offset.size(); off++) {
    for (size_t k = 0; k < twt_vec.size(); k++) {
      double tmp = offset[off] / (vrms_vec[k]*twt_vec[k] / 1000);
      thetagrid(k, off) = atan(tmp);
    }
  }
}

void SeismicForward::ResampleTwtPS(std::vector<double> &twt_pp_reg,
                                   std::vector<double> &twt_ss_reg,
                                   const std::vector<double> &twt_pp,
                                   const std::vector<double> &twt_ss,
                                   const std::vector<double> &twt,
                                   const std::vector<double> &twt_0,
                                   double twt_pp_above,
                                   double twt_ss_above,
                                   double twt_above,
                                   double twt_pp_below,
                                   double twt_ss_below,
                                   double twt_below)
{
  size_t nzrefl = twt.size();
  std::vector<double> twt_vec(nzrefl + 2), data_vec(nzrefl + 2);
  twt_vec[0] = twt_above;
  data_vec[0] = twt_pp_above;
  for (size_t i = 0; i < nzrefl; ++i) {
    twt_vec [i + 1] = twt   [i];
    data_vec[i + 1] = twt_pp[i];
  }
  twt_vec [nzrefl + 1] = twt_below;
  data_vec[nzrefl + 1] = twt_pp_below;
  twt_pp_reg = LinInterp1D(twt_vec, data_vec, twt_0);
  data_vec[0] = twt_ss_above;
  for (size_t i = 0; i < nzrefl; ++i) {
    data_vec[i + 1] = twt_ss[i];
  }
  data_vec[nzrefl + 1] = twt_ss_below;
  twt_ss_reg = LinInterp1D(twt_vec, data_vec, twt_0);
}

void SeismicForward::ResampleOffsetPS(const std::vector<double>   &twt,
                                      const NRLib::Grid2D<double> &offset_pp,
                                      const NRLib::Grid2D<double> &offset_pp_above,
                                      const NRLib::Grid2D<double> &offset_pp_below,
                                      const std::vector<double>   &twt_0,
                                      const std::vector<double>   &offset_tot_vec,
                                      NRLib::Grid2D<double>       &offset_pp_reg,
                                      NRLib::Grid2D<double>       &offset_ss_reg,
                                      double                       twt_above,
                                      double                       twt_below)
{
  size_t nzrefl = twt.size();
  size_t nsampl = twt_0.size();
  std::vector<double> off_vec_pp_reg(twt_0.size());
  std::vector<double> off_vec_pp(nzrefl + 2), twt_vec(nzrefl + 2);
  twt_vec[0] = twt_above;
  for (size_t i = 0; i < nzrefl; ++i) {
    twt_vec[i + 1] = twt[i];
  }
  twt_vec[nzrefl + 1] = twt_below;

  for (size_t off = 0; off < offset_tot_vec.size(); ++off) {
    off_vec_pp[0] = offset_pp_above(0, off);
    for (size_t i = 0; i < nzrefl; ++i) {
      off_vec_pp[i+1] = offset_pp(i, off);
    }
    off_vec_pp[nzrefl + 1] = offset_pp_below(0, off);
    off_vec_pp_reg = LinInterp1D(twt_vec, off_vec_pp, twt_0);
    for (size_t i = 0; i < nsampl; ++i) {
      offset_pp_reg(i, off) = off_vec_pp_reg[i];
      offset_ss_reg(i, off) = offset_tot_vec[off] - off_vec_pp_reg[i];
    }
  }
}

void SeismicForward::ResampleTWTx(const NRLib::Grid2D<double> &twtx_grid,
                                  const NRLib::Grid2D<double> &twtx_below,
                                  NRLib::Grid2D<double>       &twtx_grid_reg,
                                  const std::vector<double>   &twt_vec,
                                  const double                &twt_below,
                                  const std::vector<double>   &twt_0)
{
  size_t nlayers = twtx_grid.GetNI() + 1;
  size_t noff    = twtx_grid.GetNJ();
  size_t nreg    = twt_0.size();
  std::vector<double> twtx_vec_reg(nreg);
  std::vector<double> twtx_vec(nlayers);
  std::vector<double> twt_vec_nlay(nlayers);

  for (size_t i = 0; i < nlayers - 1; ++i) {
    twt_vec_nlay[i] = twt_vec[i];
  }
  twt_vec_nlay[nlayers - 1] = twt_below;

  for (size_t off = 0; off < noff; ++off) {
    for (size_t i = 0; i < nlayers - 1; ++i) {
      twtx_vec[i] = twtx_grid(i, off);
    }
    twtx_vec[nlayers - 1] = twtx_below(0, off);

    twtx_vec_reg = LinInterp1D(twt_vec_nlay, twtx_vec, twt_0);
    for (size_t i = 0; i < twtx_grid_reg.GetNI(); ++i) {
      twtx_grid_reg(i, off) = twtx_vec_reg[i];
    }
  }
}

void SeismicForward::FindTWTxPS(NRLib::Grid2D<double>       &twtx_grid,
                                const std::vector<double>   &twt_ss_vec,
                                const std::vector<double>   &twt_pp_vec,
                                const std::vector<double>   &vrms_pp_vec,
                                const std::vector<double>   &vrms_ss_vec,
                                const NRLib::Grid2D<double> &offset_ss,
                                const NRLib::Grid2D<double> &offset_pp,
                                bool                         offset_without_stretch)
{
  double twtx_pp, twtx_ss;
  for (size_t off = 0; off < offset_ss.GetNJ(); ++off) {
    for (size_t k = 0; k < twt_pp_vec.size(); ++k) {
      if (offset_without_stretch == false) {
        twtx_pp = std::sqrt(twt_pp_vec[k] * twt_pp_vec[k] / 4 + 1000 * 1000 * (offset_pp(k, off) * offset_pp(k, off) / (vrms_pp_vec[k] * vrms_pp_vec[k])));
        twtx_ss = std::sqrt(twt_ss_vec[k] * twt_ss_vec[k] / 4 + 1000 * 1000 * (offset_ss(k, off) * offset_ss(k, off) / (vrms_ss_vec[k] * vrms_ss_vec[k])));
        twtx_grid(k, off) = twtx_pp + twtx_ss;
      }
      else {
        twtx_grid(k, off) = twt_pp_vec[k] + twt_ss_vec[k]; //without NMO stretch
      }
    }
  }
}

//--------------------------------------------------------------------------
void SeismicForward::FindSeisLimits(const NRLib::Grid2D<double> & twtx_grid,
                                    const std::vector<double>   & twt_0,
                                    std::vector<size_t>         & n_min,
                                    std::vector<size_t>         & n_max,
                                    double                        twt_wave)
//--------------------------------------------------------------------------
{
  size_t i_min, i_max;
  size_t nzrefl = twtx_grid.GetNI();
  for (size_t off = 0; off < n_min.size(); ++off) {
    double twtx_min = 10000.0;
    double twtx_max =     0.0;
    for (size_t i = 0; i < nzrefl; ++i) {
      if (twtx_grid(i, off) > twtx_max)
        twtx_max = twtx_grid(i, off);
      if (twtx_grid(i, off) < twtx_min)
        twtx_min = twtx_grid(i, off);
    }
    twtx_min -= twt_wave;
    twtx_max += twt_wave;
    i_min = 0;
    i_max = twt_0.size() - 1;
    for (size_t i = 0; i < twt_0.size(); ++i) {
      if (twt_0[i] > twtx_max) {
        i_max = i;
        break;
      }
    }
    for (int i = twt_0.size() - 1; i >= 0; --i) {
      if (twt_0[static_cast<size_t>(i)] < twtx_min) {
        i_min = static_cast<size_t>(i);
        break;
      }
    }
    n_min[off] = i_min;
    n_max[off] = i_max;
  }
}

//-------------------------------------------------------------------------------
void SeismicForward::FindTWTx(NRLib::Grid2D<double>     & twtx_grid,
                              const std::vector<double> & twt_vec,
                              const std::vector<double> & vrms_vec,
                              const std::vector<double> & offset,
                              bool                        offset_without_stretch)
//-------------------------------------------------------------------------------
{
  double twtx;
  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = 0; k < twt_vec.size(); ++k) {
      if (offset_without_stretch == false) {
        double t_off = 1000.0*offset[off]/vrms_vec[k];
        twtx = twt_vec[k] * twt_vec[k] + t_off*t_off;
        twtx_grid(k, off) = std::sqrt(twtx);
      }
      else {
        twtx_grid(k,off) = twt_vec[k];
      }
    }
  }
}

//-----------------------------------------------------------------------------------------
void SeismicForward::SeisConvolutionNMO(NRLib::Grid2D<double>               & timegrid_pos,
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
                                        const std::vector<size_t>           & n_max)
//-----------------------------------------------------------------------------------------
{
  double rickerLimit = wavelet->GetTwtLength();
  size_t nt          = timegrid_pos.GetNI();
  size_t nc          = refl_pos.GetNI();
  size_t noff        = offset.size();

  double x, y, z;
  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  double topt = toptime.GetZ(x, y);

  if (toptime.IsMissing(topt) == false) {

    for (size_t off = 0 ; off < noff ; off++) {
      double t = t0;

      for (size_t k = 0 ; k < nt ; k++) {
        if (k >= n_min[off] && k <= n_max[off]) {
          double seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            double twtx_kk = twtx(kk, off);
            if (fabs(twtx_kk - t) < rickerLimit) {
              double ricker = waveletScale * wavelet->FindWaveletPoint(twtx_kk - t);
              seis += refl_pos(kk, off) * ricker;
            }
          }
          timegrid_pos(k, off) = static_cast<float>(seis);
        }
        else {
          timegrid_pos(k, off) = 0.0;
        }
        t = t + dt;
      }
    }
  }
  else {
    for (size_t k = 0 ; k < nt ; k++){
      for (size_t off = 0 ; off < noff ; off++) {
        timegrid_pos(k, off) = 0.0;
      }
    }
  }
}

//--------------------------------------------------------------------------------------
void SeismicForward::SeisConvolution(NRLib::Grid2D<double>               & timegrid_pos,
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
                                     size_t                                n_max)
//--------------------------------------------------------------------------------------
{
  double rickerLimit = wavelet->GetTwtLength();
  size_t nt          = timegrid_pos.GetNI();
  size_t nc          = refl_pos.GetNI();

  double x, y, z;
  zgrid.FindCenterOfCell(i, j, 0, x, y, z);

  double topt = toptime.GetZ(x, y);

  if (toptime.IsMissing(topt) == false) {

    for (size_t theta = 0; theta < theta_vec.size(); theta++) {
      double t = t0;

      for (size_t k = 0; k < nt; k++) {
        if (k > n_min && k < n_max) {
          double seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            double twt_kk = twt[kk];
            if (fabs(twt_kk - t) < rickerLimit) {
              double ricker = waveletScale * wavelet->FindWaveletPoint(twt_kk - t);
              seis += refl_pos(kk, theta) * ricker;
            }
          }
          timegrid_pos(k, theta) = static_cast<float>(seis);
        }
        else {
          timegrid_pos(k, theta) = 0.0;
        }
        t = t + dt;
      }
    }
  }
  else {
    for (size_t k = 0; k < nt; k++){
      for (size_t theta = 0; theta < theta_vec.size(); theta++) {
        timegrid_pos(k, theta) = 0.0;
      }
    }
  }
}
