
#include "nrlib/geometry/interpolation.hpp"
#include "nrlib/math/constants.hpp"

#include "physics/zoeppritz.hpp"
#include "physics/zoeppritz_ps.hpp"
#include "physics/zoeppritz_pp.hpp"
#include "physics/wavelet.hpp"

#include "utils/output.hpp"

#include "seismic_geometry.hpp"
#include "seismic_output.hpp"
#include "seismic_forward.hpp"


#include "tbb/compat/thread"

#include <ctime>
#include <map>

//---------------------------------------------------------------------------
void SeismicForward::DoSeismicForward(SeismicParameters & seismic_parameters)
//---------------------------------------------------------------------------
{
  if (seismic_parameters.GetModelSettings()->GetNMOCorr()) {
    MakeNMOSeismic(seismic_parameters);
  }
  else {
    MakeSeismic(seismic_parameters);
  }
}

//----------------------------------------------------------------------
void SeismicForward::MakeSeismic(SeismicParameters & seismic_parameters)
//----------------------------------------------------------------------
{
  if (seismic_parameters.GetTimeOutput() || seismic_parameters.GetDepthOutput() || seismic_parameters.GetTimeshiftOutput()) {

    bool                  ps_seis   = seismic_parameters.GetModelSettings()->GetPSSeismic();
    std::vector<double> & theta_vec = seismic_parameters.GetThetaVec();
    std::vector<double>   twt_0;
    std::vector<double>   z_0;
    std::vector<double>   twts_0;
    size_t                dummy;

    seismic_parameters.GenerateTwt0AndZ0(twt_0,
                                         z_0,
                                         twts_0,
                                         dummy,
                                         ps_seis);

    Output seis_output(seismic_parameters,
                       twt_0,
                       z_0,
                       twts_0,
                       theta_vec,
                       twt_0.size());

    time_t t1 = time(0);

    size_t                        n_traces;
    tbb::concurrent_queue<Trace*> seismic_traces = seismic_parameters.FindTracesInForward(n_traces);
    size_t                        max_threads    = seismic_parameters.GetModelSettings()->GetMaxThreads();
    size_t                        queue_capacity = seismic_parameters.GetModelSettings()->GetTracesInMemory();
    unsigned int                  n              = std::thread::hardware_concurrency();
    size_t                        n_threads      = static_cast<size_t>(n);

    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d traces to be generated.\n", n_traces);

    if (n_threads > max_threads) {
      n_threads = max_threads;
    }
    if (n_threads > 1) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d of %d available threads are used, and queue capacity is %d traces.\n",
                                  n_threads, n, queue_capacity);
    }
    else {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d of %d available threads are used.\n", n_threads, n);
    }

    PrintSeisType(false, ps_seis, theta_vec, false);

    tbb::concurrent_queue<ResultTrace*>         empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*> result_queue;
    result_queue.set_capacity(queue_capacity);
    std::vector<std::thread*>                   worker_thread;

    std::vector<double> dummy_vec;
    GenSeisTraceParams parameters(seismic_parameters,
                                  twt_0,
                                  z_0,
                                  twts_0,
                                  theta_vec,
                                  dummy_vec,
                                  empty_queue,
                                  result_queue,
                                  seismic_traces,
                                  n_traces,
                                  0);

    if (n_threads > 1){
      for (size_t i = 0; i < n_threads - 1; ++i) {
        worker_thread.push_back(new std::thread(GenerateSeismicTracesQueue, &seis_output, &parameters));
      }

      std::thread write_thread(WriteSeismicTraces, &parameters, &seis_output);

      for (size_t i = 0; i < n_threads-1; ++i) {
        worker_thread[i]->join();
        delete worker_thread[i];
      }
      write_thread.join();
    }
    else {
      float monitor_size, next_monitor;
      MonitorInitialize(n_traces, monitor_size, next_monitor);
      for (size_t i = 0; i < n_traces; ++i) {
        Trace * trace;
        parameters.seismic_traces.try_pop(trace);
        GenerateSeismicTraces(&seis_output, &parameters, trace);
        ResultTrace *result_trace;
        parameters.result_queue.try_pop(result_trace);
        WriteTrace(result_trace, parameters.seismic_parameters, &seis_output);
        Monitor(i, monitor_size, next_monitor);
        delete result_trace;
        delete trace;
      }
    }

    ResultTrace *result_trace;
    while (empty_queue.try_pop(result_trace)){
      delete result_trace;
    }

    std::cout << "\n";
    seismic_parameters.PrintElapsedTime(t1, "generating seismic");

    //write storm grid if requested
    seis_output.WriteSeismicStorm(seismic_parameters);

    seismic_parameters.DeleteZandRandTWTGrids();
    seismic_parameters.DeleteElasticParameterGrids();
    seismic_parameters.DeleteWavelet();
    seismic_parameters.DeleteGeometryAndOutput();
  }
  else {
    std::cout << "No output of seismic is requested; hence no seismic is generated.\n";
  }
}

//--------------------------------------------------------------
void SeismicForward::MakeNMOSeismic(SeismicParameters & seispar)
//--------------------------------------------------------------
{
  if (seispar.GetTimeOutput() || seispar.GetDepthOutput() || seispar.GetTimeshiftOutput()) {
    bool                  ps_seis                = seispar.GetModelSettings()->GetPSSeismic();
    bool                  offset_without_stretch = seispar.GetModelSettings()->GetOffsetWithoutStretch();
    std::vector<double> & offset_vec             = seispar.GetOffsetVec();
    size_t                time_samples_stretch   = seispar.GetSeismicGeometry()->nt();
    std::vector<double>   twt_0;
    std::vector<double>   z_0;
    std::vector<double>   twts_0;

    seispar.GenerateTwt0AndZ0(twt_0,
                              z_0,
                              twts_0,
                              time_samples_stretch,
                              ps_seis);

    Output nmo_output(seispar,
                      twt_0,
                      z_0,
                      twts_0,
                      offset_vec,
                      time_samples_stretch);

    time_t t1 = time(0);

    size_t                        n_traces;
    tbb::concurrent_queue<Trace*> seismic_traces = seispar.FindTracesInForward(n_traces);
    size_t                        max_threads    = seispar.GetModelSettings()->GetMaxThreads();
    size_t                        queue_capacity = seispar.GetModelSettings()->GetTracesInMemory();
    unsigned int                  n              = std::thread::hardware_concurrency();
    size_t                        n_threads      = static_cast<size_t>(n);

    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d traces to be generated.\n", n_traces);

    if (n_threads > max_threads) {
      n_threads = max_threads;
    }
    if (n_threads > 1) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d of %d available threads are used, and queue capacity is %d traces.\n",
                                  n_threads, n, queue_capacity);
    }
    else {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n%d of %d available threads are used.\n", n_threads, n);
    }

    PrintSeisType(true, ps_seis, offset_vec, offset_without_stretch);

    tbb::concurrent_queue<ResultTrace*>         empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*> result_queue;
    result_queue.set_capacity(queue_capacity);
    std::vector<std::thread*>                   worker_thread;

    std::vector<double> dummy_vec;
    GenSeisTraceParams parameters(seispar,
                                  twt_0,
                                  z_0,
                                  twts_0,
                                  dummy_vec,
                                  offset_vec,
                                  empty_queue,
                                  result_queue,
                                  seismic_traces,
                                  n_traces,
                                  time_samples_stretch);

    if (n_threads > 1) {
      //loop over available threads, send work to each thread
      for (size_t i = 0; i < n_threads - 1; ++i) {
        worker_thread.push_back(new std::thread(GenerateNMOSeismicTracesQueue, &nmo_output, &parameters));
      }
      std::thread write_thread(WriteSeismicTraces, &parameters, &nmo_output);

      for (size_t i = 0; i < n_threads-1; ++i) {
        worker_thread[i]->join();
        delete worker_thread[i];
      }
      write_thread.join();
    }
    else {
      float monitor_size, next_monitor;
      MonitorInitialize(n_traces, monitor_size, next_monitor);
      for (size_t i = 0; i < n_traces; ++i) {
        Trace * trace;
        parameters.seismic_traces.try_pop(trace);
        GenerateNMOSeismicTraces(&nmo_output, &parameters, trace);
        ResultTrace *result_trace;
        parameters.result_queue.try_pop(result_trace);
        WriteTrace(result_trace,
                   parameters.seismic_parameters,
                   &nmo_output);
        Monitor(i, monitor_size, next_monitor);
        delete result_trace;
        delete trace;
      }
    }

    ResultTrace *result_trace;
    while (empty_queue.try_pop(result_trace)){
      delete result_trace;
    }

    std::cout << "\n";
    seispar.PrintElapsedTime(t1, "generating seismic");

    //write storm grid if requested
    nmo_output.WriteSeismicStorm(seispar);

    seispar.DeleteZandRandTWTGrids();
    seispar.DeleteElasticParameterGrids();
    seispar.DeleteWavelet();
    seispar.DeleteGeometryAndOutput();
  }
  else {
    std::cout << "No output of seismic is requested; hence no seismic is generated.\n";
  }
}

void SeismicForward::GenerateNMOSeismicTracesQueue(Output             * nmo_output,
                                                   GenSeisTraceParams * param)
{
  Trace *trace;
  while (param->seismic_traces.try_pop(trace)) {
    GenerateNMOSeismicTraces(nmo_output, param, trace);
  }
}

//----------------------------------------------------------------------------
void SeismicForward::GenerateNMOSeismicTraces(Output             * nmo_output,
                                              GenSeisTraceParams * param,
                                              Trace              * trace)
//----------------------------------------------------------------------------
{
  ResultTrace * result_trace;
  if (!param->empty_queue.try_pop(result_trace)){
    result_trace = new ResultTrace(param->seismic_parameters,
                                   param->twt_0,
                                   param->z_0,
                                   param->twts_0,
                                   param->time_samples_stretch,
                                   param->offset_vec);
  }

  result_trace->SetJobID(trace);

  size_t i = trace->GetI();
  size_t j = trace->GetJ();
  double x = trace->GetX();
  double y = trace->GetY();

  if (GenerateTraceOk(param->seismic_parameters, i, j)) {
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

    std::vector<NRLib::StormContGrid> & rgridvec                    = param->seismic_parameters.GetRGrids();
    NRLib::RegularSurface<double>     & toptime                     = param->seismic_parameters.GetTopTime();
    NRLib::StormContGrid              & zgrid                       = param->seismic_parameters.GetZGrid();
    NRLib::StormContGrid              & twtgrid                     = param->seismic_parameters.GetTwtGrid();
    NRLib::StormContGrid              & vpgrid                      = param->seismic_parameters.GetVpGrid();

    size_t                              nx                          = param->seismic_parameters.GetSeismicGeometry()->nx();
    double                              dt                          = param->seismic_parameters.GetSeismicGeometry()->dt();
    double                              tmin                        = param->twt_0[0];
    size_t                              nzrefl                      = param->seismic_parameters.GetSeismicGeometry()->zreflectorcount();
    std::vector<double>                 constvp                     = param->seismic_parameters.GetModelSettings()->GetConstVp();
    std::vector<double>                 constvs                     = param->seismic_parameters.GetModelSettings()->GetConstVs();
    unsigned long                       seed                        = param->seismic_parameters.GetModelSettings()->GetSeed();
    bool                                ps_seis                     = param->seismic_parameters.GetModelSettings()->GetPSSeismic();

    double                              wavelet_scale               = param->seismic_parameters.GetWaveletScale();
    Wavelet                           * wavelet                     = param->seismic_parameters.GetWavelet();
    double                              twt_wavelet                 = wavelet->GetTwtLength();
    double                              z_wavelet_bot               = param->seismic_parameters.GetModelSettings()->GetZWaveletBot();
    double                              z_extrapol_factor           = param->seismic_parameters.GetModelSettings()->GetZExtrapolFactor();
    bool                                offset_without_stretch      = param->seismic_parameters.GetModelSettings()->GetOffsetWithoutStretch();

    std::vector<size_t> n_min(param->offset_vec.size());
    std::vector<size_t> n_max(param->offset_vec.size());
    std::vector<double> twt_vec(nzrefl);
    std::vector<double> vp_vec(nzrefl);

    //get twt and vp at all layers from twtgrid and vpgrid
    for (size_t k = 0; k < nzrefl; ++k) {
      twt_vec[k] = twtgrid(i, j, k);
      vp_vec[k]  = vpgrid(i, j, k);
    }

    if (ps_seis) { //------------PS seismic------------
      //setup and get vectors and grid for calculation
      std::vector<double> twt_ss_vec(nzrefl);
      std::vector<double> twt_pp_vec(nzrefl);
      std::vector<double> twt_ss_vec_reg(param->twt_0.size());
      std::vector<double> twt_pp_vec_reg(param->twt_0.size());
      std::vector<double> vs_vec(nzrefl);

      std::vector<double> vrms_pp_vec(nzrefl);
      std::vector<double> vrms_pp_vec_reg(param->twt_0.size());
      std::vector<double> vrms_ss_vec(nzrefl);
      std::vector<double> vrms_ss_vec_reg(param->twt_0.size());

      NRLib::Grid2D<double> & offset_pp     = result_trace->GetOffsetPP();
      NRLib::Grid2D<double> & offset_ss     = result_trace->GetOffsetSS();
      NRLib::Grid2D<double> & offset_pp_reg = result_trace->GetOffsetPPReg();
      NRLib::Grid2D<double> & offset_ss_reg = result_trace->GetOffsetSSReg();

      NRLib::StormContGrid  & vsgrid        = param->seismic_parameters.GetVsGrid();
      NRLib::StormContGrid  & twtssgrid     = param->seismic_parameters.GetTwtSSGrid();
      NRLib::StormContGrid  & twtppgrid     = param->seismic_parameters.GetTwtPPGrid();

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

      double z_w             = param->seismic_parameters.GetModelSettings()->GetZw();
      double v_w             = param->seismic_parameters.GetModelSettings()->GetVw();
      double twt_w           = 2000 * z_w /v_w; //twt_above pp, ps, ss

      //resample twt_pp and twt_ss
      ResampleTwtPS(twt_pp_vec_reg, twt_ss_vec_reg, twt_pp_vec, twt_ss_vec, twt_vec, param->twt_0, twt_w, twt_w, twt_w, twt_pp_below, twt_ss_below, twt_below);

      //find vrms for pp and ss -  for each layer and regularly sampled
      param->seismic_parameters.FindVrms(vrms_pp_vec, vrms_pp_vec_reg, twt_pp_vec, twt_pp_vec_reg, vp_vec, constvp[2], twt_pp_extrapol, i, j, true);
      param->seismic_parameters.FindVrms(vrms_ss_vec, vrms_ss_vec_reg, twt_ss_vec, twt_ss_vec_reg, vs_vec, constvs[2], twt_ss_extrapol, i, j, true);

      //find theta and offset - for each layer for each offset:
      param->seismic_parameters.FindPSNMOThetaAndOffset(theta_pos, offset_pp, offset_ss, twt_pp_vec, twt_ss_vec, vrms_pp_vec, vrms_ss_vec, param->offset_vec);

      //find offset above and below reservoir - and resample regularly
      NRLib::Grid2D<double> offset_pp_above(1, param->offset_vec.size()), offset_ss_above(1, param->offset_vec.size());
      NRLib::Grid2D<double> offset_pp_below(1, param->offset_vec.size()), offset_ss_below(1, param->offset_vec.size());
      // - above reservoir
      std::vector<double>   twt_pp_one (1, twt_pp_vec_reg[0]);  //twt_w
      std::vector<double>   twt_ss_one (1, twt_ss_vec_reg[0]);  //twt_w
      std::vector<double>   vrms_pp_one(1, vrms_pp_vec_reg[0]); //v_w
      std::vector<double>   vrms_ss_one(1, vrms_ss_vec_reg[0]); //v_w
      param->seismic_parameters.FindPSNMOThetaAndOffset(dummygrid, offset_pp_above, offset_ss_above, twt_pp_one, twt_ss_one, vrms_pp_one, vrms_ss_one, param->offset_vec, false);
      // - below reservoir
      twt_pp_one [0] = twt_pp_below;
      twt_ss_one [0] = twt_ss_below;
      vrms_pp_one[0] = std::sqrt(1 / twt_pp_below * (constvp[2] * constvp[2] * (twt_pp_below - twt_pp_vec[nzrefl - 1]) + vrms_pp_vec[nzrefl - 1] * vrms_pp_vec[nzrefl - 1] * twt_pp_vec[nzrefl - 1]));
      vrms_ss_one[0] = std::sqrt(1 / twt_ss_below * (constvs[2] * constvs[2] * (twt_ss_below - twt_ss_vec[nzrefl - 1]) + vrms_ss_vec[nzrefl - 1] * vrms_ss_vec[nzrefl - 1] * twt_ss_vec[nzrefl - 1]));
      param->seismic_parameters.FindPSNMOThetaAndOffset(dummygrid, offset_pp_below, offset_ss_below, twt_pp_one, twt_ss_one, vrms_pp_one, vrms_ss_one, param->offset_vec, false);
      // - resample regularly
      ResampleOffsetPS(twt_vec, offset_pp, offset_pp_above, offset_pp_below, param->twt_0, param->offset_vec, offset_pp_reg, offset_ss_reg, param->twt_0[0], twt_below);

      //find twtx - for each layer for each offset - and regularly sampled
      FindTWTxPS(twtx,     twt_ss_vec,     twt_pp_vec,     vrms_pp_vec,     vrms_ss_vec,     offset_ss,     offset_pp,     offset_without_stretch);
      FindTWTxPS(twtx_reg, twt_ss_vec_reg, twt_pp_vec_reg, vrms_pp_vec_reg, vrms_ss_vec_reg, offset_ss_reg, offset_pp_reg, offset_without_stretch);

      //find limits for where to generate seismic, for each offset
      FindSeisLimits(twtx, param->twt_0, n_min, n_max, twt_wavelet);

    }
    else { //------------PP seismic------------
      std::vector<double> vrms_vec(nzrefl);
      std::vector<double> vrms_vec_reg(param->twt_0);
      double twt_wavelet_extrapol = twt_wavelet * z_extrapol_factor;
      param->seismic_parameters.FindVrms(vrms_vec, vrms_vec_reg, twt_vec, param->twt_0, vp_vec, constvp[2], twt_wavelet_extrapol, i, j, true);

      //find theta - for each layer for each offset:
      FindNMOTheta(theta_pos, twt_vec, vrms_vec, param->offset_vec);

      //find twtx for each layer for each offset, and regularly in time:
      FindTWTx(twtx,     twt_vec,      vrms_vec,     param->offset_vec, offset_without_stretch);
      FindTWTx(twtx_reg, param->twt_0, vrms_vec_reg, param->offset_vec, offset_without_stretch);

      //find limits for where to generate seismic, for each offset
      FindSeisLimits(twtx, param->twt_0, n_min, n_max, twt_wavelet);
    }//------------------------

    //find reflection coeff - for each layer for each offset:
    param->seismic_parameters.FindNMOReflections(refl_pos, theta_pos, i, j);

    //keep reflections for zero offset if output on storm
    if (param->seismic_parameters.GetModelSettings()->GetOutputReflections()){
      for (size_t k = 0; k < nzrefl; ++k) {
        rgridvec[0](i,j,k) = float(refl_pos(k,0));
      }
    }
    //add noise to reflections
    if (param->seismic_parameters.GetModelSettings()->GetWhiteNoise()) {
      double deviation = param->seismic_parameters.GetModelSettings()->GetStandardDeviation();
      param->seismic_parameters.AddNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos);
      //keep reflections for zero offset if output on storm and white noise
      if (param->seismic_parameters.GetModelSettings()->GetOutputReflections()) {
        for (size_t k = 0; k < nzrefl; ++k) {
          rgridvec[1](i,j,k) = float(refl_pos(k,0));
        }
      }
    }

    //generate seismic
    SeisConvolutionNMO(timegrid_pos,
                       refl_pos,
                       twtx,
                       zgrid,
                       toptime,
                       wavelet,
                       wavelet_scale,
                       param->offset_vec,
                       tmin,
                       dt,
                       i,
                       j,
                       n_min,
                       n_max);

    //NMO correction:
    size_t max_sample;
    if (offset_without_stretch) {
      max_sample       = nmo_timegrid_pos.GetNI();
      nmo_timegrid_pos = NRLib::Grid2D<double>(timegrid_pos);
    }
    else {
      NMOCorrect(param->twt_0, timegrid_pos, twtx_reg, nmo_timegrid_pos, n_min, n_max, max_sample);
    }

    //stacking of offsets:
    if (param->seismic_parameters.GetStackOutput() || param->seismic_parameters.GetStormOutput()) {
      float noffset_inv = static_cast<float>(1.0 / param->offset_vec.size());
      for (size_t k = 0; k < nmo_timegrid_stack_pos.GetNI(); ++k) {
        nmo_timegrid_stack_pos(k,0) = 0.0;
        for (size_t off = 0; off < param->offset_vec.size(); ++off) {
          nmo_timegrid_stack_pos(k,0) += noffset_inv * nmo_timegrid_pos(k,off);
        }
      }
    }

    //depth conversion:
    if (nmo_output->GetDepthSegyOk() || nmo_output->GetDepthStackSegyOk() || param->seismic_parameters.GetDepthStormOutput()) {
      std::vector<double> zgrid_vec_extrapol(nzrefl+2);
      std::vector<double> twt_vec_extrapol(nzrefl+2);
      ExtrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, constvp[2], constvs[2], z_wavelet_bot*z_extrapol_factor, i, j, ps_seis);
      if (nmo_output->GetDepthSegyOk()) {
        ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, nmo_timegrid_pos, nmo_depthgrid_pos, max_sample);
      }
      if (nmo_output->GetDepthStackSegyOk() || param->seismic_parameters.GetDepthStormOutput()){
        ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, nmo_timegrid_stack_pos, nmo_depthgrid_stack_pos, max_sample);
      }
    }

    //timeshift:
    if (nmo_output->GetTimeshiftSegyOk() || nmo_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()) {
      std::vector<double> timeshiftgrid_vec_extrapol(nzrefl+2);
      std::vector<double> twt_vec_extrapol          (nzrefl+2);
      timeshiftgrid_vec_extrapol[0] = 0;
      twt_vec_extrapol          [0] = 0;
      NRLib::StormContGrid &twt_timeshift = param->seismic_parameters.GetTwtShiftGrid();
      for (size_t k = 0; k < nzrefl; ++k) {
        timeshiftgrid_vec_extrapol[k+1] = twt_timeshift(i,j,k);
        twt_vec_extrapol[k+1]           = twt_vec[k];
      }
      timeshiftgrid_vec_extrapol[nzrefl + 1] = twt_timeshift(i,j,nzrefl-1) + z_extrapol_factor * twt_wavelet *twt_timeshift(i, j, nzrefl - 1)/ twt_vec[nzrefl - 1];
      twt_vec_extrapol          [nzrefl + 1] = twt_vec[nzrefl -1]          + z_extrapol_factor * twt_wavelet;
      if (nmo_output->GetTimeshiftSegyOk()) {
        ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, nmo_timegrid_pos, nmo_timeshiftgrid_pos, max_sample);
      }
      if (nmo_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()){
        ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, nmo_timegrid_stack_pos, nmo_timeshiftgrid_stack_pos, max_sample);
      }
    }

    if (false) {
      //if (i == 0 && j == 0) {
      std::vector<double> z_vector(zgrid.GetNK());
      for (size_t ii = 0; ii < zgrid.GetNK(); ++ii) {
        z_vector[ii] = zgrid(i, j, ii);
      }
      param->seismic_parameters.GetSeismicOutput()->PrintVector(z_vector, "z_vector.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintVector(twt_vec, "twt_vec.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(twtx, "twtx.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(twtx_reg, "twtx_reg.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintVector(param->twt_0, "twt_0.txt");

      param->seismic_parameters.GetSeismicOutput()->PrintVector(vp_vec, "vp_vec.txt");

      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(refl_pos, "refl_pos.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(theta_pos, "theta_pos.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintVector(param->offset_vec, "offset_vec.txt");

      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(timegrid_pos, "timegrid_pos.txt");
      param->seismic_parameters.GetSeismicOutput()->PrintMatrix(nmo_timegrid_pos, "nmo_timegrid_pos.txt");
    }
    result_trace->SetIsEmpty(false);
    param->result_queue.push(result_trace);
  }
  else {
    result_trace->SetIsEmpty(true);
    param->result_queue.push(result_trace);
  }
}

void SeismicForward::GenerateSeismicTracesQueue(Output             *seis_output,
                                                GenSeisTraceParams *param)
{
  Trace *trace;
  while (param->seismic_traces.try_pop(trace)) {
    GenerateSeismicTraces(seis_output, param, trace);
  }
}

void SeismicForward::GenerateSeismicTraces(Output             *seis_output,
                                           GenSeisTraceParams *param,
                                           Trace              *trace)
{
  ResultTrace *result_trace;
  if (!param->empty_queue.try_pop(result_trace)){
    result_trace = new ResultTrace(param->seismic_parameters, param->twt_0, param->z_0, param->twts_0, param->twt_0.size(), param->theta_vec);
  }

  result_trace->SetJobID(trace);

  size_t i = trace->GetI();
  size_t j = trace->GetJ();
  double x = trace->GetX();
  double y = trace->GetY();

  if (GenerateTraceOk(param->seismic_parameters, i, j)) {
    NRLib::Grid2D<double> &timegrid_pos            = result_trace->GetTimeTrace();
    NRLib::Grid2D<double> &timegrid_stack_pos      = result_trace->GetTimeStackTrace();
    NRLib::Grid2D<double> &depthgrid_pos           = result_trace->GetDepthTrace();
    NRLib::Grid2D<double> &depthgrid_stack_pos     = result_trace->GetDepthStackTrace();
    NRLib::Grid2D<double> &timeshiftgrid_pos       = result_trace->GetTimeShiftTrace();
    NRLib::Grid2D<double> &timeshiftgrid_stack_pos = result_trace->GetTimeShiftStackTrace();

    std::vector<NRLib::StormContGrid> &rgridvec   = param->seismic_parameters.GetRGrids();
    size_t nx                   = param->seismic_parameters.GetSeismicGeometry()->nx();
    size_t nt                   = param->seismic_parameters.GetSeismicGeometry()->nt();
    double dt                   = param->seismic_parameters.GetSeismicGeometry()->dt();
    double tmin                 = param->twt_0[0];
    size_t nzrefl               = param->seismic_parameters.GetSeismicGeometry()->zreflectorcount();
    std::vector<double> constvp = param->seismic_parameters.GetModelSettings()->GetConstVp();
    std::vector<double> constvs = param->seismic_parameters.GetModelSettings()->GetConstVs();
    unsigned long seed          = param->seismic_parameters.GetModelSettings()->GetSeed();
    bool ps_seis                = param->seismic_parameters.GetModelSettings()->GetPSSeismic();

    size_t n_min = 0;
    size_t n_max = nt;

    NRLib::Grid2D<double> refl_pos(nzrefl, param->theta_vec.size());
    std::vector<double>   twt_vec (nzrefl);

    double wavelet_scale        = param->seismic_parameters.GetWaveletScale();
    Wavelet * wavelet           = param->seismic_parameters.GetWavelet();
    double z_wavelet_bot        = param->seismic_parameters.GetModelSettings()->GetZWaveletBot();

    NRLib::RegularSurface<double> &toptime        = param->seismic_parameters.GetTopTime();
    NRLib::StormContGrid          &zgrid          = param->seismic_parameters.GetZGrid();
    NRLib::StormContGrid          &twtgrid        = param->seismic_parameters.GetTwtGrid();

    //get twt at all layers from twtgrid
    for (size_t k = 0; k < nzrefl; ++k) {
      twt_vec[k]   = twtgrid(i,j,k);
    }
    param->seismic_parameters.FindReflections(refl_pos, param->theta_vec, i, j);

    //keep reflections for zero offset if output on storm
    if (param->seismic_parameters.GetModelSettings()->GetOutputReflections()){
      for (size_t k = 0; k < nzrefl; ++k) {
        rgridvec[0](i,j,k) = float(refl_pos(k,0));
      }
    }
    //add noise to reflections
    if (param->seismic_parameters.GetModelSettings()->GetWhiteNoise()) {
      double deviation = param->seismic_parameters.GetModelSettings()->GetStandardDeviation();
      param->seismic_parameters.AddNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos); //nb, make unique seed when i and j loop is made
      //keep reflections for zero offset if output on storm and white noise
      if (param->seismic_parameters.GetModelSettings()->GetOutputReflections()) {
        for (size_t k = 0; k < nzrefl; ++k) {
          rgridvec[1](i,j,k) = float(refl_pos(k,0));
        }
      }
    }

    //generate seismic
    SeisConvolution(timegrid_pos,
                    refl_pos,
                    twt_vec,
                    zgrid,
                    toptime,
                    wavelet,
                    wavelet_scale,
                    param->theta_vec,
                    tmin,
                    dt,
                    i,
                    j,
                    n_min,
                    n_max);

    //stacking of angles:
    if (param->seismic_parameters.GetStackOutput() || param->seismic_parameters.GetStormOutput()) {
      float ntheta_inv = static_cast<float>(1.0 / param->theta_vec.size());
      for (size_t k = 0; k < param->twt_0.size(); ++k) {
        timegrid_stack_pos(k,0) = 0.0;
        for (size_t off = 0; off < param->theta_vec.size(); ++off) {
          timegrid_stack_pos(k,0) += ntheta_inv * timegrid_pos(k,off);
        }
      }
    }

    //depth conversion:
    if (seis_output->GetDepthSegyOk() || seis_output->GetDepthStackSegyOk() || param->seismic_parameters.GetDepthStormOutput()) {
      std::vector<double> zgrid_vec_extrapol(nzrefl+2);
      std::vector<double> twt_vec_extrapol(nzrefl+2);
      ExtrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, constvp[2], constvs[2], z_wavelet_bot, i, j, ps_seis);
      if (seis_output->GetDepthSegyOk()) {
        ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, timegrid_pos, depthgrid_pos, timegrid_pos.GetNI());
      }
      if (seis_output->GetDepthStackSegyOk() || param->seismic_parameters.GetDepthStormOutput()){
        ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, timegrid_stack_pos, depthgrid_stack_pos, timegrid_stack_pos.GetNI());
      }
    }
    //timeshift:
    if (seis_output->GetTimeshiftSegyOk() || seis_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()) {
      std::vector<double> timeshiftgrid_vec_extrapol(nzrefl+1);
      std::vector<double> twt_vec_extrapol(nzrefl+1);
      timeshiftgrid_vec_extrapol[0] = 0;
      twt_vec_extrapol[0]           = 0;
      NRLib::StormContGrid &twt_timeshift = param->seismic_parameters.GetTwtShiftGrid();
      for (size_t k = 0; k < nzrefl; ++k) {
        timeshiftgrid_vec_extrapol[k+1] = twt_timeshift(i,j,k);
        twt_vec_extrapol[k+1]           = twt_vec[k];
      }
      if (seis_output->GetTimeshiftSegyOk()) {
        ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, timegrid_pos, timeshiftgrid_pos, timegrid_pos.GetNI());
      }
      if (seis_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()){
        ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, timegrid_stack_pos, timeshiftgrid_stack_pos, timegrid_stack_pos.GetNI());
      }
    }
    result_trace->SetIsEmpty(false);
    param->result_queue.push(result_trace);
  }
  else {
    result_trace->SetIsEmpty(true);
    param->result_queue.push(result_trace);
  }
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
      if (ps_seis)
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic PS-seismic for offsets: ");
      else
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic PP-seismic for offsets: ");
    }
    else {
      if (ps_seis)
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PS-seismic for offsets: ");
      else
        NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PP-seismic for offsets: ");
    }
    for (size_t i = 0; i < off_theta_vec.size(); ++i) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "%.1f ", off_theta_vec[i]);
    }
  }
  else {
    if (ps_seis)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PS-seismic for angles: ");
    else
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGenerating synthetic NMO PP-seismic for angles: ");
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

bool SeismicForward::GenerateTraceOk(SeismicParameters &seismic_parameters,
                                     size_t i,
                                     size_t j)
{
  bool generate_ok = false;
  double const_vp  = seismic_parameters.GetModelSettings()->GetConstVp()[1];
  double const_vs  = seismic_parameters.GetModelSettings()->GetConstVs()[1];
  double const_rho = seismic_parameters.GetModelSettings()->GetConstRho()[1];

  NRLib::StormContGrid &vpgrid    = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid &vsgrid    = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid &rhogrid   = seismic_parameters.GetRhoGrid();
  NRLib::StormContGrid &twtgrid   = seismic_parameters.GetTwtGrid();
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

void SeismicForward::ExtrapolZandTwtVec(std::vector<double>        &zgrid_vec_extrapol,
                                        std::vector<double>        &twt_vec_extrapol,
                                        const std::vector<double>  &twt_vec,
                                        const NRLib::StormContGrid &zgrid,
                                        double                      vp_bot,
                                        double                      vs_bot,
                                        double                      z_wavelet_bot,
                                        size_t                      i,
                                        size_t                      j,
                                        bool                        ps_seis)
{
  double vel_bot;
  if (ps_seis)
    vel_bot = 2 / (1 / vp_bot + 1 / vs_bot);
  else
    vel_bot = vp_bot;
  size_t nzrefl = zgrid_vec_extrapol.size() - 2;
  zgrid_vec_extrapol[0] = 0;
  twt_vec_extrapol  [0] = 0;
  for (size_t k = 0; k < nzrefl; ++k) {
    twt_vec_extrapol  [k+1] = twt_vec[k];
    zgrid_vec_extrapol[k+1] = zgrid(i, j, k);
  }
  zgrid_vec_extrapol[nzrefl+1] = zgrid_vec_extrapol[nzrefl] + z_wavelet_bot;
  twt_vec_extrapol  [nzrefl+1] = twt_vec_extrapol[nzrefl]   + 2000 * z_wavelet_bot / vel_bot;
}

void SeismicForward::ConvertSeis(const std::vector<double>   &twt_vec,
                                 const std::vector<double>   &twt_0,
                                 const std::vector<double>   &zgrid_vec,
                                 const std::vector<double>   &z_0,
                                 const NRLib::Grid2D<double> &seismic,
                                 NRLib::Grid2D<double>       &conv_seismic,
                                 const size_t                &max_sample)
{
  size_t nk = conv_seismic.GetNI();
  std::vector<double> seismic_vec(max_sample);
  std::vector<double> conv_seismic_vec(nk);

  std::vector<double> zt_reg = LinInterp1D(twt_vec, zgrid_vec, twt_0);
  zt_reg.resize(max_sample);

  for (size_t off = 0; off < seismic.GetNJ(); off++) {
    for (size_t k = 0; k < max_sample; k++) {
      seismic_vec[k] = seismic(k, off);
    }
    conv_seismic_vec = SplineInterp1D(zt_reg, seismic_vec, z_0, 0);
    for (size_t k = 0; k < nk; k++) {
      conv_seismic(k, off) = conv_seismic_vec[k];
    }
  }
}

void SeismicForward::NMOCorrect(const std::vector<double>   &t_in,
                                const NRLib::Grid2D<double> &data_in,
                                const NRLib::Grid2D<double> &t_out,
                                NRLib::Grid2D<double>       &data_out,
                                const std::vector<size_t>   &n_min,
                                const std::vector<size_t>   &n_max,
                                size_t                      &max_sample)
{
  max_sample = 0;
  size_t nt_in = data_in.GetNI();
  size_t noff  = data_in.GetNJ();
  std::vector<double> data_vec_in(nt_in), t_vec_in(nt_in), t_vec_out(nt_in), data_vec_out(nt_in);
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

    //not necessary to interpolate AT values higher than max t
    //t_out not monotonously increasing, must check that we are inside
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
    if (index > max_sample) {
      max_sample = index;
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

void SeismicForward::FindSeisLimits(const NRLib::Grid2D<double> &twtx_grid,
                                    const std::vector<double>   &twt_0,
                                    std::vector<size_t>         &n_min,
                                    std::vector<size_t>         &n_max,
                                    double                       twt_wave)
{
  size_t i_min, i_max;
  size_t nzrefl = twtx_grid.GetNI();
  for (size_t off = 0; off < n_min.size(); ++off) {
    double twtx_min = 10000, twtx_max = 0;
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

void SeismicForward::FindTWTx(NRLib::Grid2D<double>     &twtx_grid,
                              const std::vector<double> &twt_vec,
                              const std::vector<double> &vrms_vec,
                              const std::vector<double> &offset,
                              bool                       offset_without_stretch)
{
  double twtx;
  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = 0; k < twt_vec.size(); ++k) {
      if (offset_without_stretch == false) {
        twtx = twt_vec[k] * twt_vec[k] + 1000 * 1000 * (offset[off] * offset[off] / (vrms_vec[k] * vrms_vec[k]));
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
  size_t nt = timegrid_pos.GetNI();
  size_t nc = refl_pos.GetNI();
  double seis, twtx_kk;
  double x, y, z;
  double topt        = 0.0;
  double rickerLimit = wavelet->GetTwtLength();

  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  topt = toptime.GetZ(x, y);
  if (toptime.IsMissing(topt) == false) {
    for (size_t off = 0; off < offset.size(); off++) {
      double t = t0;

      for (size_t k = 0; k < nt; k++) {
        if ((k > n_min[off] || k == n_min[off]) && k < (n_max[off] + 1)) {
          seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            twtx_kk = twtx(kk, off);
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
    for (size_t k = 0; k < nt; k++){
      for (size_t off = 0; off < offset.size(); off++) {
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

  size_t nt = timegrid_pos.GetNI();
  size_t nc = refl_pos.GetNI();
  double seis, twt_kk;
  double x, y, z;
  double topt        = 0.0;
  double rickerLimit = wavelet->GetTwtLength();

  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  topt = toptime.GetZ(x, y);
  if (toptime.IsMissing(topt) == false) {
    for (size_t theta = 0; theta < theta_vec.size(); theta++) {
      double t = t0;

      for (size_t k = 0; k < nt; k++) {
        if (k > n_min && k < n_max) {
          seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            twt_kk = twt[kk];
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

void SeismicForward::WriteSeismicTraces(GenSeisTraceParams *param,
                                        Output             *seis_output)
{
  float monitor_size, next_monitor;
  MonitorInitialize(param->n_traces, monitor_size, next_monitor);
  size_t trace = 0;
  std::map<size_t, ResultTrace*> finished_jobs;
  std::map<size_t, ResultTrace*>::iterator it;
  while (trace < param->n_traces) {
    ResultTrace *result_trace;
    if (param->result_queue.try_pop(result_trace)) {
      finished_jobs.insert(std::pair<size_t, ResultTrace*>(result_trace->GetJobNumber(), result_trace));
    }
    it = finished_jobs.find(trace);
    while (it != finished_jobs.end()) {
      WriteTrace(finished_jobs[trace], param->seismic_parameters, seis_output);
      param->empty_queue.push(finished_jobs[trace]);
      finished_jobs.erase(it);
      Monitor(trace, monitor_size, next_monitor);
      ++trace;
      it = finished_jobs.find(trace);
    }
  }
}

void SeismicForward::WriteTrace(ResultTrace       *result_trace,
                                SeismicParameters &seismic_parameters,
                                Output            *seis_output)
{
  if (result_trace->GetIsEmpty()) {
    seis_output->AddZeroTrace(seismic_parameters,
      result_trace->GetX(),
      result_trace->GetY(),
      result_trace->GetI(),
      result_trace->GetJ());
  }
  else {
    seis_output->AddTrace(seismic_parameters,
      result_trace->GetTimeTrace(),
      result_trace->GetPreNMOTimeTrace(),
      result_trace->GetTimeStackTrace(),
      result_trace->GetDepthTrace(),
      result_trace->GetDepthStackTrace(),
      result_trace->GetTimeShiftTrace(),
      result_trace->GetTimeShiftStackTrace(),
      result_trace->GetTWTxReg(),
      result_trace->GetX(),
      result_trace->GetY(),
      result_trace->GetI(),
      result_trace->GetJ());
  }
}
