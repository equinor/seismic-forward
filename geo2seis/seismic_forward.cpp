#include <seismic_forward.hpp>
#include <seismic_regridding.hpp>

#include <physics/wavelet.hpp>

#include <seismic_geometry.hpp>
#include <seismic_output.hpp>
//#include <seismic_output.hpp>
#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>
#include "nrlib/geometry/interpolation.hpp"
#include <ctime>
//#include <thread>
//#include "tbb/concurrent_queue.h"
#include "tbb/compat/thread"

void SeismicForward::DoSeismicForward(SeismicParameters &seismic_parameters) {

  if (seismic_parameters.modelSettings()->GetNMOCorr()) {
    MakeNMOSeismic(seismic_parameters);
  }
  else {
    MakeSeismic(seismic_parameters);
  }
}

void SeismicForward::MakeSeismic(SeismicParameters &seismic_parameters) 
{
  if (seismic_parameters.GetTimeOutput() || seismic_parameters.GetDepthOutput() || seismic_parameters.GetTimeshiftOutput()) {

    bool ps_seis                   = seismic_parameters.modelSettings()->GetPSSeismic();
    std::vector<double> & theta_vec  = seismic_parameters.GetThetaVec();
    std::vector<double> twt_0, z_0, twts_0;
    size_t dummy;
    seismic_parameters.GenerateTwt0AndZ0(twt_0, z_0, twts_0, dummy);

    //prepare segy and storm files:
    Output seis_output(seismic_parameters, twt_0, z_0, twts_0, theta_vec, twt_0.size());

    time_t t1 = time(0);   // get time now

    ///parallelisation

    size_t n_traces;
    tbb::concurrent_queue<Trace*> seismic_traces = FindTracesInForward(seismic_parameters, n_traces);

    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";

    PrintSeisType(false, ps_seis, theta_vec);

    size_t n_threads = 3;
    size_t queue_capacity = 100;//50 * n_threads;

    tbb::concurrent_queue<ResultTrace*> empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*> result_queue;
    result_queue.set_capacity(queue_capacity);

    std::vector<std::thread> worker_thread;
    std::vector<std::thread> write_thread;

    std::vector<double> dummy_vec;
    Output * output = &seis_output;
    GenSeisTraceParams parameters_tmp(seismic_parameters,
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
    GenSeisTraceParams * parameters = &parameters_tmp;
    for (size_t i = 0; i < n_threads; ++i) {
      worker_thread.push_back(std::thread(GenerateSeismicTraces, output, parameters));
    }

    write_thread.push_back(std::thread(WriteSeismicTraces, parameters, output));

    for (size_t i = 0; i < n_threads; ++i) {
      worker_thread[i].join();
    }
    write_thread[0].join();


    ResultTrace *result_trace;
    while (empty_queue.try_pop(result_trace)){
      delete result_trace;
      ResultTrace *result_trace;
    }
    //delete result_trace;
    //delete output;
    //delete parameters;

    PrintElapsedTime(t1);

    //write storm grid if requested
    seis_output.WriteSeismicStorm(seismic_parameters);

    seismic_parameters.deleteZandRandTWTGrids();
    seismic_parameters.deleteElasticParameterGrids();
    seismic_parameters.deleteWavelet();
    seismic_parameters.deleteGeometryAndOutput();
  }
}

void SeismicForward::MakeNMOSeismic(SeismicParameters &seismic_parameters)
{
  if (seismic_parameters.GetTimeOutput() || seismic_parameters.GetDepthOutput() || seismic_parameters.GetTimeshiftOutput()) {

    bool ps_seis                = seismic_parameters.modelSettings()->GetPSSeismic();
    std::vector<double> & offset_vec = seismic_parameters.GetOffsetVec();
    std::vector<double> twt_0, z_0, twts_0;
    size_t time_samples_stretch;
    seismic_parameters.GenerateTwt0AndZ0(twt_0, z_0, twts_0, time_samples_stretch);

    //prepare segy and storm files
    Output nmo_output(seismic_parameters, twt_0, z_0, twts_0, offset_vec, time_samples_stretch);

    time_t t1 = time(0);   // get time now

    ///parallelisation
    size_t n_traces;
    tbb::concurrent_queue<Trace*> seismic_traces = FindTracesInForward(seismic_parameters, n_traces);

    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";

    PrintSeisType(true, ps_seis, offset_vec);

    size_t n_threads = 3;
    size_t queue_capacity = 100;//50 * n_threads;

    tbb::concurrent_queue<ResultTrace*> empty_queue;
    tbb::concurrent_bounded_queue<ResultTrace*> result_queue;
    result_queue.set_capacity(queue_capacity);

    std::vector<std::thread> worker_thread;
    std::vector<std::thread> write_thread;

    Output * output = &nmo_output;
    std::vector<double> dummy_vec;
    GenSeisTraceParams parameters_tmp(seismic_parameters,
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
    GenSeisTraceParams *parameters = & parameters_tmp;


    //loop over available threads, send work to each thread
    for (size_t i = 0; i < n_threads; ++i) {
      worker_thread.push_back(std::thread(GenerateNMOSeismicTraces, output, parameters));
    }
    write_thread.push_back(std::thread(WriteSeismicTraces, parameters, output));

    for (size_t i = 0; i < n_threads; ++i) {
      worker_thread[i].join();
    }
    write_thread[0].join();

    ResultTrace *result_trace;
    while (empty_queue.try_pop(result_trace)){
      delete result_trace;
      ResultTrace *result_trace;
    }
    //delete result_trace;
    //delete output;
    //delete parameters;

    PrintElapsedTime(t1);

    //write storm grid if requested
    nmo_output.WriteSeismicStorm(seismic_parameters);

    seismic_parameters.deleteZandRandTWTGrids();
    seismic_parameters.deleteElasticParameterGrids();
    seismic_parameters.deleteWavelet();
    seismic_parameters.deleteGeometryAndOutput();
  }
}


tbb::concurrent_queue<Trace*> SeismicForward::FindTracesInForward(SeismicParameters &seismic_parameters,
                                                                  size_t            &n_traces)
{
  tbb::concurrent_queue<Trace*> traces;
  int n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step;
  bool ilxl_loop = false;
  //find index min and max and whether loop over i,j or il,xl:
  seismic_parameters.findLoopIndeces(n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step, ilxl_loop);
  NRLib::SegyGeometry *geometry = seismic_parameters.segyGeometry();
  int il_steps = 0;
  int xl_steps = 0;
  //----------------------LOOP OVER I,J OR IL,XL---------------------------------
  size_t job_number = 0;
  for (int il = il_min; il <= il_max; il += il_step) { 
    //for (int il = 1350; il < 1452; il += il_step) {
    ++il_steps;
    xl_steps = 0;
    for (int xl = xl_min; xl <= xl_max; xl +=xl_step) {
    //for (int xl = 1280; xl < 1281; xl += xl_step) {
      ++xl_steps;
      size_t i, j;
      double x, y;
      if (ilxl_loop) { //loop over il,xl, find corresponding x,y,i,j
        geometry->FindXYFromILXL(il, xl, x, y);
        geometry->FindIndex(x, y, i, j);
      }
      else { //loop over i,j, no segy output
        i = il;
        j = xl;
        x = 0;
        y = 0;
      }
      Trace * trace = new Trace(job_number, x, y, i, j);
      traces.push(trace);
      ++job_number;
    }
  }
  n_traces = job_number;
  std::cout << n_traces <<  " traces to be generated.\n";
  return traces;
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
    if (param->result_queue.try_pop(result_trace)){
      finished_jobs.insert(std::pair<size_t, ResultTrace*> (result_trace->GetJobNumber(), result_trace));
    }
    //else {
    //  what?
    //}
    it = finished_jobs.find(trace);
    while(it != finished_jobs.end()) {
      WriteTrace(finished_jobs[trace], param->seismic_parameters, seis_output);
      param->empty_queue.push(finished_jobs[trace]);
      finished_jobs.erase(it);

      Monitor(trace, monitor_size, next_monitor);
      //std::cout << "write " << trace << "\n";
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
                         result_trace->GetTWTx(),
                         result_trace->GetX(),
                         result_trace->GetY(),
                         result_trace->GetI(),
                         result_trace->GetJ());
  }
}

void SeismicForward::GenerateNMOSeismicTraces(Output             *nmo_output,
                                              GenSeisTraceParams *param)
{
  Trace *trace;
  while (param->seismic_traces.try_pop(trace)){

    ResultTrace *result_trace;
    if (!param->empty_queue.try_pop(result_trace)){
      result_trace = new ResultTrace(param->seismic_parameters, param->twt_0, param->z_0, param->twts_0, param->time_samples_stretch, param->offset_vec);
    }

    result_trace->SetJobID(trace);

    size_t i = trace->GetI();
    size_t j = trace->GetJ();
    double x = trace->GetX();
    double y = trace->GetY();


    if (GenerateTraceOk(param->seismic_parameters, i, j)) {
      NRLib::Grid2D<double> &timegrid_pos                = result_trace->GetPreNMOTimeTrace();
      NRLib::Grid2D<double> &nmo_timegrid_pos            = result_trace->GetTimeTrace();
      NRLib::Grid2D<double> &nmo_timegrid_stack_pos      = result_trace->GetTimeStackTrace();
      NRLib::Grid2D<double> &nmo_depthgrid_pos           = result_trace->GetDepthTrace();
      NRLib::Grid2D<double> &nmo_depthgrid_stack_pos     = result_trace->GetDepthStackTrace();
      NRLib::Grid2D<double> &nmo_timeshiftgrid_pos       = result_trace->GetTimeShiftTrace();
      NRLib::Grid2D<double> &nmo_timeshiftgrid_stack_pos = result_trace->GetTimeShiftStackTrace();
      NRLib::Grid2D<double> &twtx_reg                    = result_trace->GetTWTx();

      std::vector<NRLib::StormContGrid> &rgridvec        = param->seismic_parameters.rGrids();
      size_t nx     = param->seismic_parameters.seismicGeometry()->nx();
      double dt     = param->seismic_parameters.seismicGeometry()->dt();
      double tmin   = param->twt_0[0] - 0.5*dt;
      size_t nzrefl = param->seismic_parameters.seismicGeometry()->zreflectorcount();
      std::vector<double> constvp = param->seismic_parameters.modelSettings()->GetConstVp();
      std::vector<double> constvs = param->seismic_parameters.modelSettings()->GetConstVs();
      unsigned long seed          = param->seismic_parameters.modelSettings()->GetSeed();
      bool ps_seis                = param->seismic_parameters.modelSettings()->GetPSSeismic();

      std::vector<size_t> n_min(param->offset_vec.size());
      std::vector<size_t> n_max(param->offset_vec.size());
      std::vector<double> vrms_vec    (nzrefl);
      std::vector<double> vrms_vec_reg(param->twt_0.size());

          //setup vectors for twt, vrms, theta, refl and seismic
      std::vector<double>   twt_vec  (nzrefl);
      NRLib::Grid2D<double> theta_pos(nzrefl, param->offset_vec.size()); 
      NRLib::Grid2D<double> refl_pos (nzrefl, param->offset_vec.size());
      NRLib::Grid2D<double> twtx     (nzrefl, param->offset_vec.size());

      double wavelet_scale        = param->seismic_parameters.waveletScale();
      Wavelet * wavelet           = param->seismic_parameters.wavelet();
    
      NRLib::RegularSurface<double>     &bottom_eclipse = param->seismic_parameters.bottomEclipse();
      NRLib::RegularSurface<double>     &toptime        = param->seismic_parameters.topTime();
      NRLib::StormContGrid              &zgrid          = param->seismic_parameters.zGrid();
      NRLib::StormContGrid              &twtgrid        = param->seismic_parameters.twtGrid();

      //get twt at all layers from twtgrid
      for (size_t k = 0; k < nzrefl; ++k) {
        twt_vec[k]   = twtgrid(i,j,k);
      }

      //calculate vrms per reflection and regularly sampled  
      param->seismic_parameters.findVrmsPos(vrms_vec, vrms_vec_reg, param->twt_0, i, j);

      //find min and max sample for seismic - for each offset.
      param->seismic_parameters.getSeisLimits(param->twt_0, vrms_vec, param->offset_vec, n_min, n_max);

      //find theta - for each layer for each offset:
      FindNMOTheta(theta_pos, twt_vec, vrms_vec, param->offset_vec);

      //find reflection coeff - for each layer for each offset:
      param->seismic_parameters.findNMOReflections(refl_pos, theta_pos, param->offset_vec, i, j);

      //keep reflections for zero offset if output on storm
      if (param->seismic_parameters.modelSettings()->GetOutputReflections()){
        for (size_t k = 0; k < nzrefl; ++k) {
          rgridvec[0](i,j,k) = float(refl_pos(k,0));
        }
      }
      //add noise to reflections
      if (param->seismic_parameters.modelSettings()->GetWhiteNoise()) {
        double deviation = param->seismic_parameters.modelSettings()->GetStandardDeviation();
        SeismicRegridding::AddNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos); //nb, make unique seed when i and j loop is made
        //keep reflections for zero offset if output on storm and white noise
        if (param->seismic_parameters.modelSettings()->GetOutputReflections()) {
          for (size_t k = 0; k < nzrefl; ++k) {
            rgridvec[1](i,j,k) = float(refl_pos(k,0));
          }
        }
      }

      //find twtx - for each layer for each offset:
      FindTWTx(twtx, twt_vec, vrms_vec, param->offset_vec);

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

      //sample twtx regularly - for each layer for each offset:
      FindTWTx(twtx_reg, param->twt_0, vrms_vec_reg, param->offset_vec);

      //NMO correction:
      size_t max_sample;
      NMOCorrect(param->twt_0, timegrid_pos, twtx_reg, nmo_timegrid_pos, n_min, n_max, max_sample);

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
        ExtrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, bottom_eclipse.GetZ(x,y), constvp[2], constvs[2], i, j, ps_seis);
        if (nmo_output->GetDepthSegyOk()) {
          ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, nmo_timegrid_pos, nmo_depthgrid_pos, max_sample);
        }
        if (nmo_output->GetDepthStackSegyOk() || param->seismic_parameters.GetDepthStormOutput()){
          ConvertSeis(twt_vec_extrapol, param->twt_0, zgrid_vec_extrapol, param->z_0, nmo_timegrid_stack_pos, nmo_depthgrid_stack_pos, max_sample);
        }
      }

      //timeshift:
      if (nmo_output->GetTimeshiftSegyOk() || nmo_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()) {
        std::vector<double> timeshiftgrid_vec_extrapol(nzrefl+1);
        std::vector<double> twt_vec_extrapol(nzrefl+1);
        timeshiftgrid_vec_extrapol[0] = 0;
        twt_vec_extrapol[0]           = 0;
        NRLib::StormContGrid &twt_timeshift = param->seismic_parameters.twtShiftGrid();
        for (size_t k = 0; k < nzrefl; ++k) {
          timeshiftgrid_vec_extrapol[k+1] = twt_timeshift(i,j,k);
          twt_vec_extrapol[k+1]           = twt_vec[k];
        }
        if (nmo_output->GetTimeshiftSegyOk()) {
          ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, nmo_timegrid_pos, nmo_timeshiftgrid_pos, max_sample);
        }
        if (nmo_output->GetTimeshiftStackSegyOk() || param->seismic_parameters.GetTimeshiftStormOutput()){
          ConvertSeis(twt_vec_extrapol, param->twt_0, timeshiftgrid_vec_extrapol, param->twts_0, nmo_timegrid_stack_pos, nmo_timeshiftgrid_stack_pos, max_sample);
        }

      }
      result_trace->SetIsEmpty(false);
      param->result_queue.push(result_trace);
      //std::cout << "id = " << trace->GetJobNumber() << "\n";
    }
    else {
      result_trace->SetIsEmpty(true);
      param->result_queue.push(result_trace);
    }
  }
}

void SeismicForward::GenerateSeismicTraces(Output             *seis_output,
                                           GenSeisTraceParams *param)
{

  Trace *trace;
  while (param->seismic_traces.try_pop(trace)){

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

      std::vector<NRLib::StormContGrid> &rgridvec   = param->seismic_parameters.rGrids();
      size_t nx                   = param->seismic_parameters.seismicGeometry()->nx();
      size_t nt                   = param->seismic_parameters.seismicGeometry()->nt();
      double dt                   = param->seismic_parameters.seismicGeometry()->dt();
      double tmin                 = param->twt_0[0] - 0.5*dt;
      size_t nzrefl               = param->seismic_parameters.seismicGeometry()->zreflectorcount();
      std::vector<double> constvp = param->seismic_parameters.modelSettings()->GetConstVp();
      std::vector<double> constvs = param->seismic_parameters.modelSettings()->GetConstVs();
      unsigned long seed          = param->seismic_parameters.modelSettings()->GetSeed();
      bool ps_seis                = param->seismic_parameters.modelSettings()->GetPSSeismic();

      size_t n_min = 0;
      size_t n_max = nt;

      NRLib::Grid2D<double> refl_pos(nzrefl, param->theta_vec.size());
      std::vector<double>   twt_vec (nzrefl);

      double wavelet_scale        = param->seismic_parameters.waveletScale();
      Wavelet * wavelet           = param->seismic_parameters.wavelet();

      NRLib::RegularSurface<double> &bottom_eclipse = param->seismic_parameters.bottomEclipse();
      NRLib::RegularSurface<double> &toptime        = param->seismic_parameters.topTime();
      NRLib::StormContGrid          &zgrid          = param->seismic_parameters.zGrid();
      NRLib::StormContGrid          &twtgrid        = param->seismic_parameters.twtGrid();

      //get twt at all layers from twtgrid
      for (size_t k = 0; k < nzrefl; ++k) {
        twt_vec[k]   = twtgrid(i,j,k);
      }
      param->seismic_parameters.findReflections(refl_pos, param->theta_vec, i, j);

      //keep reflections for zero offset if output on storm
      if (param->seismic_parameters.modelSettings()->GetOutputReflections()){
        for (size_t k = 0; k < nzrefl; ++k) {
          rgridvec[0](i,j,k) = float(refl_pos(k,0));
        }
      }
      //add noise to reflections
      if (param->seismic_parameters.modelSettings()->GetWhiteNoise()) {
        double deviation = param->seismic_parameters.modelSettings()->GetStandardDeviation();
        SeismicRegridding::AddNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos); //nb, make unique seed when i and j loop is made
        //keep reflections for zero offset if output on storm and white noise
        if (param->seismic_parameters.modelSettings()->GetOutputReflections()) {
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
        ExtrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, bottom_eclipse.GetZ(x,y), constvp[2], constvs[2], i, j, ps_seis);
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
        NRLib::StormContGrid &twt_timeshift = param->seismic_parameters.twtShiftGrid();
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
      //std::cout << "id = " << trace->GetJobNumber() << "\n";
    }
    else {
      result_trace->SetIsEmpty(true);
      param->result_queue.push(result_trace);
    }
  }
}

void SeismicForward::MonitorInitialize(size_t n_traces,
                                       float &monitor_size,
                                       float &next_monitor)
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

void SeismicForward::Monitor(size_t trace,
                             float monitor_size,
                             float &next_monitor)
{
  //int size_loop = n_xl * il_steps + xl_steps + 1;
  if (trace >= static_cast<size_t>(next_monitor)) {
    next_monitor += monitor_size;
    std::cout << "^";
    fflush(stdout);
    if (next_monitor > monitor_size*51){
      std::cout << "\n";
    }
  }
}

void SeismicForward::PrintSeisType(bool                 nmo,
                                   bool                 ps_seis,
                                   std::vector<double> &off_theta_vec)
{
  if (nmo) {
    if (ps_seis)
      std::cout << "Generating synthetic NMO PS-seismic for offsets: ";
    else
      std::cout << "Generating synthetic NMO PP-seismic for offsets: ";
  }
  else {
    if (ps_seis)
      std::cout << "Generating synthetic PS-seismic for angles: ";
    else
      std::cout << "Generating synthetic PP-seismic for angles: ";
  }
  for (size_t i = 0; i < off_theta_vec.size(); ++i){
    std::cout << off_theta_vec[i] << " ";
  }
  std::cout << "\n";
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

void SeismicForward::PrintElapsedTime(time_t start_time)
{
  time_t end_time       = time(0);   // get time now
  size_t seconds        = static_cast<size_t>(difftime(end_time,start_time));

  size_t hours          = static_cast<int>(seconds/3600);
  seconds               = seconds % 3600;
  size_t minutes        = static_cast<int>(seconds/60);
  seconds               = seconds % 60;

  std::string hours_s   = NRLib::ToString(hours);
  std::string zeros     = std::string(2 - hours_s.length(), '0');
  hours_s               = zeros + hours_s;
  std::string minutes_s = NRLib::ToString(minutes);
  zeros                 = std::string(2 - minutes_s.length(), '0');
  minutes_s             = zeros + minutes_s;
  std::string seconds_s = NRLib::ToString(seconds);
  zeros                 = std::string(2 - seconds_s.length(), '0');
  seconds_s             = zeros + seconds_s;

  std::cout << "\nTotal time generating seismic: "
    << hours_s << ':'
    << minutes_s << ':'
    << seconds_s
    << "\n";
}

bool SeismicForward::GenerateTraceOk(SeismicParameters &seismic_parameters,
                                     size_t i,
                                     size_t j)
{
  bool generate_ok = false;
  double const_vp  = seismic_parameters.modelSettings()->GetConstVp()[1];
  double const_vs  = seismic_parameters.modelSettings()->GetConstVs()[1];
  double const_rho = seismic_parameters.modelSettings()->GetConstRho()[1];

  NRLib::StormContGrid &vpgrid    = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid    = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid   = seismic_parameters.rhoGrid();
  NRLib::StormContGrid &twtgrid   = seismic_parameters.twtGrid();
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
                                        double                      z_bot,
                                        double                      vp_bot,
                                        double                      vs_bot,
                                        size_t                      i,
                                        size_t                      j,
                                        bool                        ps_seis)
{
  double vel_bot;
  if (ps_seis)
    vel_bot = 0.5*(vp_bot + vs_bot);
  else
    vel_bot = vp_bot;
  size_t nzrefl = zgrid_vec_extrapol.size() - 2;
  zgrid_vec_extrapol[0] = 0;
  twt_vec_extrapol[0]   = 0;
  for (size_t k = 0; k < nzrefl; ++k) {
    twt_vec_extrapol[k+1]   = twt_vec[k];
    zgrid_vec_extrapol[k+1] = zgrid(i, j, k);
  }
  zgrid_vec_extrapol[nzrefl+1] = z_bot;
  twt_vec_extrapol[nzrefl+1]   = twt_vec_extrapol[nzrefl] + 2000* (zgrid_vec_extrapol[nzrefl+1] - zgrid_vec_extrapol[nzrefl]) / vel_bot;
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
      ++index;
      if (inside == true && t_vec_out[k] > t_vec_in[n_min_max-1]){
        break;
      }
    }
    t_vec_out.resize(index);
    data_vec_out = SplineInterp1D(t_vec_in, data_vec_in, t_vec_out, 0);
    if (index > data_out.GetNI()){
      std::cout << "ERROR: stretch not properly accounted for\n";
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


void SeismicForward::FindNMOTheta(NRLib::Grid2D<double>     &thetagrid,
                                  const std::vector<double> &twt_vec,
                                  const std::vector<double> &vrms_vec,
                                  const std::vector<double> &offset){

  for (size_t off = 0; off < offset.size(); off++) {
    for (size_t k = 0; k < twt_vec.size(); k++) {
      double tmp = offset[off] / (vrms_vec[k]*twt_vec[k] / 1000);
      thetagrid(k, off) = atan(tmp); 
    }
  }
}


void SeismicForward::FindTWTx(NRLib::Grid2D<double>     &twtx_grid,
                              const std::vector<double> &twt_vec,
                              const std::vector<double> &vrms_vec,
                              const std::vector<double> &offset){
  double twtx;
  for (size_t off = 0; off < offset.size(); ++off) {
    for (size_t k = 0; k < twt_vec.size(); ++k) {
      twtx = twt_vec[k]*twt_vec[k] + 1000*1000*(offset[off]*offset[off]/(vrms_vec[k]*vrms_vec[k]));
      twtx_grid(k,off) = std::sqrt(twtx);
      //twtx_grid(k,off) = twt_vec[k]; //NBNB test for not NMO corr
    }
  }
}


void SeismicForward::SeisConvolutionNMO(NRLib::Grid2D<double>               &timegrid_pos,
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
                                        const std::vector<size_t>           &n_max)
{

  size_t nt = timegrid_pos.GetNI();
  size_t nc = refl_pos.GetNI();
  double seis;
  double x, y, z;
  double topt        = 0.0;
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  size_t count_out = 0;
  size_t count_in = 0;
  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  topt = toptime.GetZ(x, y);
  if (toptime.IsMissing(topt) == false) {
    for (size_t off = 0; off < offset.size(); off++) {
      double t = t0 + 0.5 * dt;

      for (size_t k = 0; k < nt; k++) {
        if (k > n_min[off] && k < n_max[off]) {
          ++count_in;
          seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            if (fabs(twtx(kk, off) - t) < rickerLimit) {
              double ricker = waveletScale * wavelet->FindWaveletPoint(twtx(kk, off) - t);
              seis += refl_pos(kk, off) * ricker;
            }
          }
          timegrid_pos(k, off) = static_cast<float>(seis);
        }
        else {
          timegrid_pos(k, off) = 0.0;
          ++count_out;
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
  //std::cout << "in = " << count_in << ", out = " << count_out << "\n";
}

void SeismicForward::SeisConvolution(NRLib::Grid2D<double>               &timegrid_pos,
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
                                     size_t           n_max)
{

  size_t nt = timegrid_pos.GetNI();
  size_t nc = refl_pos.GetNI();
  double seis;
  double x, y, z;
  double topt        = 0.0;
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  size_t count_out = 0;
  size_t count_in = 0;
  zgrid.FindCenterOfCell(i, j, 0, x, y, z);
  topt = toptime.GetZ(x, y);
  if (toptime.IsMissing(topt) == false) {
    for (size_t theta = 0; theta < theta_vec.size(); theta++) {
      double t = t0 + 0.5 * dt;

      for (size_t k = 0; k < nt; k++) {
        if (k > n_min && k < n_max) {
          ++count_in;
          seis = 0.0;
          for (size_t kk = 0; kk < nc; kk++) {
            if (fabs(twt[kk] - t) < rickerLimit) {
              double ricker = waveletScale * wavelet->FindWaveletPoint(twt[kk] - t);
              seis += refl_pos(kk, theta) * ricker;
            }
          }
          timegrid_pos(k, theta) = static_cast<float>(seis);
        }
        else {
          timegrid_pos(k, theta) = 0.0;
          ++count_out;
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
  //std::cout << "in = " << count_in << ", out = " << count_out << "\n";
}


