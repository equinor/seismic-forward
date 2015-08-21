#include <seismic_forward.hpp>
#include <seismic_regridding.hpp>

#include <physics/wavelet.hpp>

#include <seismic_geometry.hpp>
#include <seismic_output.hpp>
#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>
#include "nrlib/geometry/interpolation.hpp"
#include "utils/nmo_output.hpp"
#include "utils/seis_output.hpp"
#include <ctime>



void SeismicForward::seismicForward(SeismicParameters &seismic_parameters) {

  if (seismic_parameters.modelSettings()->GetNMOCorr()) {
    makeNMOSeismic(seismic_parameters);
  }
  else if (seismic_parameters.modelSettings()->GetOldModelling()) {
    makeSeismicOLD(seismic_parameters);
  }
  else {
    makeSeismic(seismic_parameters);
  }
}

void SeismicForward::makeSeismic(SeismicParameters &seismic_parameters) 
{
  if (seismic_parameters.GetTimeOutput() || seismic_parameters.GetDepthOutput() || seismic_parameters.GetTimeshiftOutput()) {

    ModelSettings * model_settings = seismic_parameters.modelSettings();
    size_t nx                   = seismic_parameters.seismicGeometry()->nx();
    size_t ny                   = seismic_parameters.seismicGeometry()->ny();
    size_t nz                   = seismic_parameters.seismicGeometry()->nz();
    size_t nt                   = seismic_parameters.seismicGeometry()->nt();
    size_t nzrefl               = seismic_parameters.seismicGeometry()->zreflectorcount();
    std::vector<double> constvp = seismic_parameters.modelSettings()->GetConstVp();
    std::vector<double> constvs = seismic_parameters.modelSettings()->GetConstVs();
    unsigned long seed          = seismic_parameters.modelSettings()->GetSeed();
    bool ps_seis                = seismic_parameters.modelSettings()->GetPSSeismic();

    NRLib::StormContGrid              &zgrid          = seismic_parameters.zGrid();
    NRLib::StormContGrid              &twtgrid        = seismic_parameters.twtGrid();
    NRLib::RegularSurface<double>     &bottom_eclipse = seismic_parameters.bottomEclipse();
    std::vector<NRLib::StormContGrid> &rgridvec       = seismic_parameters.rGrids();

    //find twt_0 and z_0
    double tmin = seismic_parameters.seismicGeometry()->t0();
    double dt   = seismic_parameters.seismicGeometry()->dt();
    std::vector<double> twt_0(nt);
    for (size_t i = 0; i < nt; ++i){
      twt_0[i] = tmin + (0.5 + i)*dt;
    }
    double zmin = seismic_parameters.seismicGeometry()->z0();
    double dz   = seismic_parameters.seismicGeometry()->dz();
    std::vector<double> z_0(nz);
    for (size_t i = 0; i < nz; ++i){
      z_0[i] = zmin + (0.5 + i)*dz;
    }
    std::vector<double> & theta_vec  = seismic_parameters.GetThetaVec(); 
    
    //setup vectors for twt, vrms, theta, refl and seismic
    std::vector<double> twt_vec     (nzrefl);   

    NRLib::Grid2D<double> timegrid_pos (twt_0.size(), theta_vec.size());
    NRLib::Grid2D<double> timegrid_stack_pos;
    NRLib::Grid2D<double> timeshiftgrid_pos;
    NRLib::Grid2D<double> timeshiftgrid_stack_pos;
    NRLib::Grid2D<double> depthgrid_pos;
    NRLib::Grid2D<double> depthgrid_stack_pos;

    std::vector<double> twts_0;
    if (model_settings->GetTwtFileName() != "") {
      twts_0 = seismic_parameters.GenerateTWT0Shift(twt_0[0], twt_0.size());
    }

    if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
      timegrid_stack_pos.Resize(twt_0.size(), 1);
    }
    if (seismic_parameters.GetTimeshiftOutput()) {
      timeshiftgrid_pos.      Resize(twts_0.size(), theta_vec.size());
      timeshiftgrid_stack_pos.Resize(twts_0.size(), 1);
    }
    if (seismic_parameters.GetDepthOutput()){
      depthgrid_pos.      Resize(z_0.size(), theta_vec.size());
      depthgrid_stack_pos.Resize(z_0.size(), 1);
    }

    //prepare segy and storm files:      
    SeisOutput seis_output(seismic_parameters, twt_0, z_0, twts_0);

    printTime();
    std::cout << "Generating synthetic seismic for angles: ";
    for (size_t i = 0; i < theta_vec.size(); ++i){
      std::cout << theta_vec[i]/NRLib::Degree << " ";
    }
    std::cout << "\n";

    float monitor_size, next_monitor;
    monitorInitialize(nx, ny, monitor_size, next_monitor);
    int n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step;
    bool ilxl_loop = false;
    //find index min and max and whether loop over i,j or il,xl:
    seismic_parameters.findLoopIndeces(n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step, ilxl_loop);
    NRLib::SegyGeometry *geometry = seismic_parameters.segyGeometry();
    int il_steps = 0;
    int xl_steps = 0;
    //----------------------LOOP OVER I,J OR IL,XL---------------------------------
    for (int il = il_min; il <= il_max; il += il_step) {
    //for (int il = 1350; il < 1351; il += il_step) {
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
        //----------------------BEGIN GEN SEIS WITH NMO FOR I,J---------------------------------
        if (generateTraceOk(seismic_parameters, i, j)) {
          //get twt at all layers from twtgrid
          for (size_t k = 0; k < nzrefl; ++k) {
            twt_vec[k]   = twtgrid(i,j,k);
          }
          generateSeismicTrace(seismic_parameters,
                               twt_vec,
                               twt_0,
                               theta_vec,
                               timegrid_pos,
                               i,j,
                               seed,
                               zgrid,
                               rgridvec);

          //stacking of angles:
          if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
            float ntheta_inv = static_cast<float>(1.0 / theta_vec.size());
            for (size_t k = 0; k < twt_0.size(); ++k) {
              timegrid_stack_pos(k,0) = 0.0;
              for (size_t off = 0; off < theta_vec.size(); ++off) {
                timegrid_stack_pos(k,0) += ntheta_inv * timegrid_pos(k,off);
              }
            }
          }

          //depth conversion:
          if (seis_output.GetDepthSegyOk() || seis_output.GetDepthStackSegyOk() || seismic_parameters.GetDepthStormOutput()) {
            std::vector<double> zgrid_vec_extrapol(nzrefl+2);
            std::vector<double> twt_vec_extrapol(nzrefl+2);
            extrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, bottom_eclipse.GetZ(x,y), constvp[2], constvs[2], i, j, ps_seis);
            if (seis_output.GetDepthSegyOk()) {
              convertSeis(twt_vec_extrapol, twt_0, zgrid_vec_extrapol, z_0, timegrid_pos, depthgrid_pos, timegrid_pos.GetNI());
            }
            if (seis_output.GetDepthStackSegyOk() || seismic_parameters.GetDepthStormOutput()){
              convertSeis(twt_vec_extrapol, twt_0, zgrid_vec_extrapol, z_0, timegrid_stack_pos, depthgrid_stack_pos, timegrid_stack_pos.GetNI());
            }
          }

          //timeshift:
          if (seis_output.GetTimeshiftSegyOk() || seis_output.GetTimeshiftStackSegyOk() || seismic_parameters.GetTimeshiftStormOutput()) {
            std::vector<double> timeshiftgrid_vec_extrapol(nzrefl+1);
            std::vector<double> twt_vec_extrapol(nzrefl+1);
            timeshiftgrid_vec_extrapol[0] = 0;
            twt_vec_extrapol[0]           = 0;
            NRLib::StormContGrid &twt_timeshift = seismic_parameters.twtShiftGrid();
            for (size_t k = 0; k < nzrefl; ++k) {
              timeshiftgrid_vec_extrapol[k+1] = twt_timeshift(i,j,k);
              twt_vec_extrapol[k+1]           = twt_vec[k];
            }
            if (seis_output.GetTimeshiftSegyOk()) {
              convertSeis(twt_vec_extrapol, twt_0, timeshiftgrid_vec_extrapol, twts_0, timegrid_pos, timeshiftgrid_pos, timegrid_pos.GetNI());
            }
            if (seis_output.GetTimeshiftStackSegyOk() || seismic_parameters.GetTimeshiftStormOutput()){
              convertSeis(twt_vec_extrapol, twt_0, timeshiftgrid_vec_extrapol, twts_0, timegrid_stack_pos, timeshiftgrid_stack_pos, timegrid_stack_pos.GetNI());
            } 
          }

          //print trace to segy and store in storm grid
          seis_output.AddTrace(seismic_parameters,
                               timegrid_pos,
                               timegrid_stack_pos,
                               depthgrid_pos,
                               depthgrid_stack_pos,
                               timeshiftgrid_pos,
                               timeshiftgrid_stack_pos,
                               x, y, i, j);
          }
        //-----------------OUTSIDE ECLIPSE GRID, WHERE NO SEISMIC GENERATED, ZERO TRACE-----------------
        else {
          //print zero trace to segy and store zero trace in storm grid
          seis_output.AddZeroTrace(seismic_parameters, x, y, i, j);
        }
        monitor(n_xl, il_steps, xl_steps, monitor_size, next_monitor);

      }

    }
    printTime();

    //write storm grid if requested
    seis_output.WriteSeismicStorm(seismic_parameters);
    seismic_parameters.deleteZandRandTWTGrids();
    seismic_parameters.deleteElasticParameterGrids();
    seismic_parameters.deleteWavelet();
    seismic_parameters.deleteGeometryAndOutput();
  }
}

void SeismicForward::makeSeismicOLD(SeismicParameters &seismic_parameters) 
{

  ModelSettings * model_settings = seismic_parameters.modelSettings();
  size_t nx     = seismic_parameters.seismicGeometry()->nx();
  size_t ny     = seismic_parameters.seismicGeometry()->ny();
  size_t nz     = seismic_parameters.seismicGeometry()->nz();
  size_t nt     = seismic_parameters.seismicGeometry()->nt();
  double dt     = seismic_parameters.seismicGeometry()->dt();
  size_t nzrefl = seismic_parameters.seismicGeometry()->zreflectorcount();
  size_t ntheta = seismic_parameters.nTheta();
  double tmin   = seismic_parameters.seismicGeometry()->t0();
  double dz     = seismic_parameters.seismicGeometry()->dz();
  double d1     = seismic_parameters.seismicGeometry()->z0();
  
  double wavelet_scale   = seismic_parameters.waveletScale();
  Wavelet * wavelet      = seismic_parameters.wavelet();  

  NRLib::StormContGrid              &zgrid          = seismic_parameters.zGrid();
  NRLib::StormContGrid              &twtgrid        = seismic_parameters.twtGrid();
  std::vector<NRLib::StormContGrid> &rgridvec       = seismic_parameters.rGrids();
  NRLib::RegularSurface<double>     &toptime        = seismic_parameters.topTime();
  NRLib::RegularSurface<double>     &bottom_eclipse = seismic_parameters.bottomEclipse();

  NRLib::Volume volume   = seismic_parameters.seismicGeometry()->createDepthVolume();
  NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();

  std::vector<double> constvp = seismic_parameters.modelSettings()->GetConstVp();

  bool time_output            = seismic_parameters.GetTimeOutput();
  bool depth_output           = seismic_parameters.GetDepthOutput();
  bool timeshift_output       = seismic_parameters.GetTimeshiftOutput();

  NRLib::StormContGrid *twt_timeshift;
  if (model_settings->GetTwtFileName() != "") {
    twt_timeshift = new NRLib::StormContGrid(model_settings->GetTwtFileName());
    if (twtgrid.GetNI() != nx || twtgrid.GetNJ() != ny || twtgrid.GetNK() != nzrefl) {
      printf("TWT from file has wrong dimension. Aborting. \n");
      exit(1);
    }
  } else {
    twt_timeshift = NULL;
  }

  printTime();

  double memory_use     = static_cast<double>(4.0 * (nx * ny * nzrefl * (2 + ntheta) + nx * ny * nz * ntheta * depth_output + nx * ny * nt * ntheta * time_output + (0.5 * nx * ny * nz))); // last term is a buffer
  double memory_limit   = model_settings->GetMemoryLimit();

  ///-----------------------------------------------------------------------------------------------
  ///----------------------------MEMORY USE < MEMORY LIMIT------------------------------------------
  ///-----------------------------------------------------------------------------------------------
  if (memory_use < memory_limit) {
    // ---- Create grids for seismic data
    std::vector<NRLib::StormContGrid> timegridvec(ntheta);
    std::vector<NRLib::StormContGrid> depthgridvec(ntheta);
    std::vector<NRLib::StormContGrid> timeshiftgridvec(ntheta);

    NRLib::StormContGrid timegrid(0, 0, 0);
    if (time_output) {
      timegrid = NRLib::StormContGrid(volume_t, nx, ny, nt);
    }
    NRLib::StormContGrid timeshiftgrid(0, 0, 0);
    if (timeshift_output) {
      timeshiftgrid = NRLib::StormContGrid(volume_t, nx, ny, nt);
    }
    NRLib::StormContGrid depthgrid(0, 0, 0);
    if (depth_output) {
      depthgrid = NRLib::StormContGrid(volume, nx, ny, nz);
    }

    for (size_t i = 0; i < ntheta; i++) {
      timegridvec[i] = timegrid;
      depthgridvec[i] = depthgrid;
      timeshiftgridvec[i] = timeshiftgrid;
    }

    //-----------------Generate seismic----------------------------
    if (model_settings->GetPSSeismic()) {
      printf("\nGenerating PS seismic:\n");
    }

    generateSeismic(rgridvec,
      twtgrid,
      zgrid,
      *twt_timeshift,
      timegridvec,
      depthgridvec,
      timeshiftgridvec,
      wavelet,
      dt,
      bottom_eclipse,
      toptime,
      tmin, dz, d1,
      constvp,
      wavelet_scale,
      time_output,
      depth_output,
      timeshift_output);


    seismic_parameters.deleteZandRandTWTGrids();
    if (twt_timeshift != NULL) {
      delete twt_timeshift;
    }

    printTime();

    if (model_settings->GetOutputTimeSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicTimeSegy(seismic_parameters, timegridvec);
    }

    if (model_settings->GetOutputSeismicTime()) {
      seismic_parameters.seismicOutput()->writeSeismicTimeStorm(seismic_parameters, timegridvec);
    }

    if (model_settings->GetOutputTimeshiftSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicTimeshiftSegy(seismic_parameters, timeshiftgridvec);
    }

    if (model_settings->GetOutputSeismicTimeshift()) {
      seismic_parameters.seismicOutput()->writeSeismicTimeshiftStorm(seismic_parameters, timeshiftgridvec);
    }

    if (model_settings->GetOutputDepthSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicDepthSegy(seismic_parameters, depthgridvec);
    }

    if (model_settings->GetOutputSeismicDepth()) {
      seismic_parameters.seismicOutput()->writeSeismicDepthStorm(seismic_parameters, depthgridvec);
    }

    if (model_settings->GetOutputSeismicStackTimeStorm() || model_settings->GetOutputSeismicStackTimeSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicStackTime(seismic_parameters, timegridvec);
    }

    if (model_settings->GetOutputSeismicStackTimeShiftStorm() || model_settings->GetOutputSeismicStackTimeShiftSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicStackTimeshift(seismic_parameters, timeshiftgridvec);
    }

    if (model_settings->GetOutputSeismicStackDepthStorm() || model_settings->GetOutputSeismicStackDepthSegy()) {
      seismic_parameters.seismicOutput()->writeSeismicStackDepth(seismic_parameters, depthgridvec);
    }
  }
  
  ///-----------------------------------------------------------------------------------------------
  ///----------------------------END(MEMORY USE < MEMORY LIMIT)-------------------------------------
  ///-----------------------------------------------------------------------------------------------

  ///-----------------------------------------------------------------------------------------------
  ///-------------------------------MEMORY USE TOO BIG----------------------------------------------
  ///-----------------------------------------------------------------------------------------------
  else {

    generateSeismicOnFile(rgridvec,
      twtgrid,
      zgrid,
      *twt_timeshift,
      wavelet,
      dt,
      nt, nz, nx, ny,
      bottom_eclipse,
      toptime,
      tmin, dz, d1,
      constvp,
      wavelet_scale,
      time_output,
      depth_output,
      timeshift_output);



    seismic_parameters.deleteZandRandTWTGrids();
    if (twt_timeshift != NULL) {
      delete twt_timeshift;
    }
    seismic_parameters.seismicOutput()->writeSeismicTimeSeismicOnFile(seismic_parameters, time_output);
    seismic_parameters.seismicOutput()->writeSeismicDepthSeismicOnFile(seismic_parameters, depth_output);
  }
  //-----------------------------------------------------------------------------------------------
  //-------------------------------END( MEMORY USE TOO BIG)----------------------------------------
  //-----------------------------------------------------------------------------------------------
  seismic_parameters.deleteWavelet();
  seismic_parameters.deleteGeometryAndOutput();
}



void SeismicForward::makeNMOSeismic(SeismicParameters &seismic_parameters)
{
  if (seismic_parameters.GetTimeOutput() || seismic_parameters.GetDepthOutput() || seismic_parameters.GetTimeshiftOutput()) {
    ModelSettings * model_settings = seismic_parameters.modelSettings();
    size_t nx                   = seismic_parameters.seismicGeometry()->nx();
    size_t ny                   = seismic_parameters.seismicGeometry()->ny();
    size_t nz                   = seismic_parameters.seismicGeometry()->nz();
    size_t nt                   = seismic_parameters.seismicGeometry()->nt();
    size_t nzrefl               = seismic_parameters.seismicGeometry()->zreflectorcount();
    std::vector<double> constvp = seismic_parameters.modelSettings()->GetConstVp();
    std::vector<double> constvs = seismic_parameters.modelSettings()->GetConstVs();
    unsigned long seed          = seismic_parameters.modelSettings()->GetSeed();
    bool ps_seis                = seismic_parameters.modelSettings()->GetPSSeismic();

    NRLib::StormContGrid              &zgrid          = seismic_parameters.zGrid();
    NRLib::StormContGrid              &twtgrid        = seismic_parameters.twtGrid();
    NRLib::RegularSurface<double>     &bottom_eclipse = seismic_parameters.bottomEclipse();
    std::vector<NRLib::StormContGrid> &rgridvec       = seismic_parameters.rGrids();

    //find max twt for seismic grid - should handle the largest offset.
    size_t time_samples_stretch;
    std::vector<double> twt_0        = seismic_parameters.GenerateTwt0ForNMO(time_samples_stretch);
    std::vector<double> z_0          = seismic_parameters.GenerateZ0ForNMO();
    std::vector<double> & offset_vec = seismic_parameters.GetOffsetVec();
    //seismic_parameters.seismicOutput()->printVector(twt_0, "twt_0.txt");
    //seismic_parameters.seismicOutput()->printVector(z_0, "z_0.txt");
    
    //setup vectors for twt, vrms, theta, refl and seismic
    std::vector<double> twt_vec     (nzrefl);   

    NRLib::Grid2D<double> twtx_reg              (twt_0.size(),         offset_vec.size());
    NRLib::Grid2D<double> timegrid_pos          (twt_0.size(),         offset_vec.size());
    NRLib::Grid2D<double> nmo_timegrid_pos      (time_samples_stretch, offset_vec.size());
    NRLib::Grid2D<double> nmo_timegrid_stack_pos;
    NRLib::Grid2D<double> nmo_timeshiftgrid_pos;
    NRLib::Grid2D<double> nmo_timeshiftgrid_stack_pos;
    NRLib::Grid2D<double> nmo_depthgrid_pos;
    NRLib::Grid2D<double> nmo_depthgrid_stack_pos;

    std::vector<double> twts_0;
    if (model_settings->GetTwtFileName() != "") {
      twts_0 = seismic_parameters.GenerateTWT0Shift(twt_0[0], time_samples_stretch);
    }

    if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
      nmo_timegrid_stack_pos.Resize(time_samples_stretch, 1);
    }
    if (seismic_parameters.GetTimeshiftOutput()) {
      nmo_timeshiftgrid_pos.      Resize(twts_0.size(), offset_vec.size());
      nmo_timeshiftgrid_stack_pos.Resize(twts_0.size(), 1);
    }
    if (seismic_parameters.GetDepthOutput()){
      nmo_depthgrid_pos.      Resize(z_0.size(), offset_vec.size());
      nmo_depthgrid_stack_pos.Resize(z_0.size(), 1);
    }

    //prepare segy and storm files
    NMOOutput nmo_output(seismic_parameters, twt_0, z_0, twts_0, time_samples_stretch);

    printTime();
    std::cout << "Generating synthetic NMO seismic for offsets: ";
    for (size_t i = 0; i < offset_vec.size(); ++i){
      std::cout << offset_vec[i] << " ";
    }
    std::cout << "\n";
    float monitor_size, next_monitor;
    monitorInitialize(nx, ny, monitor_size, next_monitor);
    int n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step;
    bool ilxl_loop = false;
    //find index min and max and whether loop over i,j or il,xl:
    seismic_parameters.findLoopIndeces(n_xl, il_min, il_max, il_step, xl_min, xl_max, xl_step, ilxl_loop);
    NRLib::SegyGeometry *geometry = seismic_parameters.segyGeometry();

    int il_steps = 0;
    int xl_steps = 0;
    //----------------------LOOP OVER I,J OR IL,XL---------------------------------      
    for (int il = il_min; il <= il_max; il += il_step) {
    //for (int il = 1350; il < 1352; il += il_step) {
      ++il_steps;
      xl_steps = 0;
      for (int xl = xl_min; xl <= xl_max; xl +=xl_step) { 
      //for (int xl = 1280; xl < 1282; xl += xl_step) {        
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
        //std::cout << "il, xl, i, j, x, y = " << il << " " << xl << " " << i << " " << j << " " << x << " " << y << "\n";
        //----------------------BEGIN GEN SEIS WITH NMO FOR I,J---------------------------------
        if (generateTraceOk(seismic_parameters, i, j)) {

          //get twt at all layers from twtgrid
          for (size_t k = 0; k < nzrefl; ++k) {
            twt_vec[k]   = twtgrid(i,j,k);
          }
          size_t max_sample;
          generateNMOSeismicTrace(seismic_parameters,
                                  twt_vec,
                                  twt_0,
                                  offset_vec,
                                  timegrid_pos,
                                  nmo_timegrid_pos,
                                  twtx_reg,
                                  i,j,
                                  seed,
                                  zgrid,
                                  rgridvec,
                                  max_sample);


          //stacking of offsets:
          if (seismic_parameters.GetStackOutput() || seismic_parameters.GetStormOutput()) {
            float noffset_inv = static_cast<float>(1.0 / offset_vec.size());
            for (size_t k = 0; k < nmo_timegrid_stack_pos.GetNI(); ++k) {
              nmo_timegrid_stack_pos(k,0) = 0.0;
              for (size_t off = 0; off < offset_vec.size(); ++off) {
                nmo_timegrid_stack_pos(k,0) += noffset_inv * nmo_timegrid_pos(k,off);
              }
            }
          }

          //depth conversion:
          if (nmo_output.GetNMODepthSegyOk() || nmo_output.GetNMODepthStackSegyOk() || seismic_parameters.GetDepthStormOutput()) {
            std::vector<double> zgrid_vec_extrapol(nzrefl+2);
            std::vector<double> twt_vec_extrapol(nzrefl+2);
            extrapolZandTwtVec(zgrid_vec_extrapol, twt_vec_extrapol, twt_vec, zgrid, bottom_eclipse.GetZ(x,y), constvp[2], constvs[2], i, j, ps_seis);
            if (nmo_output.GetNMODepthSegyOk()) {
              convertSeis(twt_vec_extrapol, twt_0, zgrid_vec_extrapol, z_0, nmo_timegrid_pos, nmo_depthgrid_pos, max_sample);
            }
            if (nmo_output.GetNMODepthStackSegyOk() || seismic_parameters.GetDepthStormOutput()){
              convertSeis(twt_vec_extrapol, twt_0, zgrid_vec_extrapol, z_0, nmo_timegrid_stack_pos, nmo_depthgrid_stack_pos, max_sample);
            }
          }

          //timeshift:
          if (nmo_output.GetNMOTimeshiftSegyOk() || nmo_output.GetNMOTimeshiftStackSegyOk() || seismic_parameters.GetTimeshiftStormOutput()) {
            std::vector<double> timeshiftgrid_vec_extrapol(nzrefl+1);
            std::vector<double> twt_vec_extrapol(nzrefl+1);
            timeshiftgrid_vec_extrapol[0] = 0;
            twt_vec_extrapol[0]           = 0;
            NRLib::StormContGrid &twt_timeshift = seismic_parameters.twtShiftGrid();
            for (size_t k = 0; k < nzrefl; ++k) {
              timeshiftgrid_vec_extrapol[k+1] = twt_timeshift(i,j,k);
              twt_vec_extrapol[k+1]           = twt_vec[k];
            }
            if (nmo_output.GetNMOTimeshiftSegyOk()) {
              convertSeis(twt_vec_extrapol, twt_0, timeshiftgrid_vec_extrapol, twts_0, nmo_timegrid_pos, nmo_timeshiftgrid_pos, max_sample);
            }
            if (nmo_output.GetNMOTimeshiftStackSegyOk() || seismic_parameters.GetTimeshiftStormOutput()){
              convertSeis(twt_vec_extrapol, twt_0, timeshiftgrid_vec_extrapol, twts_0, nmo_timegrid_stack_pos, nmo_timeshiftgrid_stack_pos, max_sample);
            } 
          }

          //print trace to segy and store in storm grid
          nmo_output.AddTrace(seismic_parameters,
                              timegrid_pos,
                              nmo_timegrid_pos,
                              nmo_timegrid_stack_pos,
                              nmo_depthgrid_pos,
                              nmo_depthgrid_stack_pos,
                              nmo_timeshiftgrid_pos,
                              nmo_timeshiftgrid_stack_pos,
                              twtx_reg,
                              x, y, i, j);
          }
        //-----------------OUTSIDE ECLIPSE GRID, WHERE NO SEISMIC GENERATED, ZERO TRACE-----------------
        else {
          //print zero trace to segy and store zero trace in storm grid
          nmo_output.AddZeroTrace(seismic_parameters, x, y, i, j);
        }
        monitor(n_xl, il_steps, xl_steps, monitor_size, next_monitor);
      }
    }

    printTime();

    //write storm grid if requested
    nmo_output.WriteSeismicStorm(seismic_parameters);

    seismic_parameters.deleteZandRandTWTGrids();
    seismic_parameters.deleteElasticParameterGrids();
    seismic_parameters.deleteWavelet();
    seismic_parameters.deleteGeometryAndOutput();
  }
}

void SeismicForward::generateNMOSeismicTrace(SeismicParameters             &seismic_parameters,
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
                                             size_t                        &max_sample)
{
  size_t nx     = seismic_parameters.seismicGeometry()->nx();
  double dt     = seismic_parameters.seismicGeometry()->dt();
  double tmin   = twt_0[0] - 0.5*dt;
  size_t nzrefl = seismic_parameters.seismicGeometry()->zreflectorcount();

  std::vector<size_t> n_min(offset_vec.size());
  std::vector<size_t> n_max(offset_vec.size());
  std::vector<double> vrms_vec    (nzrefl);
  std::vector<double> vrms_vec_reg(twt_0.size());

  NRLib::Grid2D<double> theta_pos(nzrefl, offset_vec.size()); 
  NRLib::Grid2D<double> refl_pos (nzrefl, offset_vec.size());
  NRLib::Grid2D<double> twtx     (nzrefl, offset_vec.size());

  double wavelet_scale        = seismic_parameters.waveletScale();
  Wavelet * wavelet           = seismic_parameters.wavelet();
  double        deviation     = seismic_parameters.modelSettings()->GetStandardDeviation();

  NRLib::RegularSurface<double>     &toptime        = seismic_parameters.topTime();
  

  //calculate vrms per reflection and regularly sampled  
  seismic_parameters.findVrmsPos(vrms_vec, vrms_vec_reg, twt_0, i, j);

  //find min and max sample for seismic - for each offset.
  seismic_parameters.getSeisLimits(twt_0, vrms_vec, offset_vec, n_min, n_max);

  //find theta - for each layer for each offset:
  findNMOTheta(theta_pos, twt_vec, vrms_vec, offset_vec);

  //find reflection coeff - for each layer for each offset:
  seismic_parameters.findNMOReflections(refl_pos, theta_pos, offset_vec, i, j);

  //keep reflections for zero offset if output on storm
  if (seismic_parameters.modelSettings()->GetOutputReflections()){
    for (size_t k = 0; k < nzrefl; ++k) {
      rgridvec[0](i,j,k) = float(refl_pos(k,0));
    }
  }
  //add noise to reflections
  if (seismic_parameters.modelSettings()->GetWhiteNoise()) {
    SeismicRegridding::addNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos); //nb, make unique seed when i and j loop is made
    //keep reflections for zero offset if output on storm and white noise
    if (seismic_parameters.modelSettings()->GetOutputReflections()) {
      for (size_t k = 0; k < nzrefl; ++k) {
        rgridvec[1](i,j,k) = float(refl_pos(k,0));
      }
    }
  }

  //find twtx - for each layer for each offset:
  findTWTx(twtx, twt_vec, vrms_vec, offset_vec);

  //generate seismic
  seisConvolutionNMO(timegrid_pos,
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

  //sample twtx regularly - for each layer for each offset:
  findTWTx(twtx_reg, twt_0, vrms_vec_reg, offset_vec);

  //NMO correction:
  NMOCorrect(twt_0, timegrid_pos, twtx_reg, nmo_timegrid_pos, n_min, n_max, max_sample);

  //////debug print
  //seismic_parameters.seismicOutput()->printVector(twt_0, "twt_0.txt");
  //seismic_parameters.seismicOutput()->printVector(twt_vec, "twt_vec.txt");
  ////seismic_parameters.seismicOutput()->printVector(vrms_vec, "vrms_vec.txt");
  //seismic_parameters.seismicOutput()->printVector(vrms_vec_reg, "vrms_vec_reg.txt");
  ////seismic_parameters.seismicOutput()->printMatrix(theta_pos, "theta_pos.txt");
  ////seismic_parameters.seismicOutput()->printMatrix(refl_pos, "refl_pos.txt");
  //seismic_parameters.seismicOutput()->printMatrix(twtx, "twtx.txt");
  //seismic_parameters.seismicOutput()->printMatrix(timegrid_pos, "timegrid_pos.txt");
  //seismic_parameters.seismicOutput()->printMatrix(twtx_reg, "twtx_reg.txt");
  //seismic_parameters.seismicOutput()->printMatrix(nmo_timegrid_pos, "nmo_timegrid_pos.txt");
}

void SeismicForward::generateSeismicTrace(SeismicParameters             &seismic_parameters,
                                          const std::vector<double>     &twt_vec,
                                          const std::vector<double>     &twt_0,
                                          const std::vector<double>     &theta_vec,
                                          NRLib::Grid2D<double>         &timegrid_pos,
                                          size_t                         i,
                                          size_t                         j,
                                          unsigned long                  seed,
                                          const NRLib::StormContGrid    &zgrid,
                                          std::vector<NRLib::StormContGrid> &rgridvec)
{
  size_t nx     = seismic_parameters.seismicGeometry()->nx();
  size_t nt     = seismic_parameters.seismicGeometry()->nt();
  double dt     = seismic_parameters.seismicGeometry()->dt();
  double tmin   = twt_0[0] - 0.5*dt;
  size_t nzrefl = seismic_parameters.seismicGeometry()->zreflectorcount();

  size_t n_min = 0;
  size_t n_max = nt;
  
  NRLib::Grid2D<double> refl_pos (nzrefl, theta_vec.size());

  double wavelet_scale        = seismic_parameters.waveletScale();
  Wavelet * wavelet           = seismic_parameters.wavelet();

  NRLib::RegularSurface<double>     &toptime        = seismic_parameters.topTime();

  seismic_parameters.findReflections(refl_pos, theta_vec, i, j);

  //keep reflections for zero offset if output on storm
  if (seismic_parameters.modelSettings()->GetOutputReflections()){
    for (size_t k = 0; k < nzrefl; ++k) {
      rgridvec[0](i,j,k) = float(refl_pos(k,0));
    }
  }
  //add noise to reflections
  if (seismic_parameters.modelSettings()->GetWhiteNoise()) {
    double deviation = seismic_parameters.modelSettings()->GetStandardDeviation();
    SeismicRegridding::addNoiseToReflectionsPos(seed+long(i+nx*j), deviation, refl_pos); //nb, make unique seed when i and j loop is made
    //keep reflections for zero offset if output on storm and white noise
    if (seismic_parameters.modelSettings()->GetOutputReflections()) {
      for (size_t k = 0; k < nzrefl; ++k) {
        rgridvec[1](i,j,k) = float(refl_pos(k,0));
      }
    }
  }

  //generate seismic
  seisConvolution(timegrid_pos,
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

}

void SeismicForward::monitorInitialize(size_t nx,
                                       size_t ny,
                                       float &monitor_size,
                                       float &next_monitor)
{
  monitor_size = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
}

void SeismicForward::monitor(size_t n_xl,
                             size_t il_steps,
                             size_t xl_steps,
                             float monitor_size,
                             float &next_monitor)
{
  //int size_loop = n_xl * il_steps + xl_steps + 1;
  if (n_xl * il_steps + xl_steps + 1 >= static_cast<size_t>(next_monitor)) {
    next_monitor += monitor_size;
    std::cout << "^";
    fflush(stdout);
    if (next_monitor > monitor_size*51){
      std::cout << "\n";
    }
  }  
}


void SeismicForward::printTime() 
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  std::cout << "Time: "
    << (now->tm_hour) << ':' 
    << (now->tm_min) << ':'
    <<  now->tm_sec
    << "\n";
}

bool SeismicForward::generateTraceOk(SeismicParameters &seismic_parameters,
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



std::vector<double> SeismicForward::linInterp1D(const std::vector<double> &x_in,
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

std::vector<double> SeismicForward::splineInterp1D(const std::vector<double> &x_in,
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

void SeismicForward::extrapolZandTwtVec(std::vector<double>        &zgrid_vec_extrapol,
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


void SeismicForward::convertSeis(const std::vector<double>   &twt_vec,
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

  std::vector<double> zt_reg = linInterp1D(twt_vec, zgrid_vec, twt_0);
  zt_reg.resize(max_sample);

  for (size_t off = 0; off < seismic.GetNJ(); off++) {
    for (size_t k = 0; k < max_sample; k++) {
      seismic_vec[k] = seismic(k, off);
    }
    conv_seismic_vec = splineInterp1D(zt_reg, seismic_vec, z_0, 0);
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
    data_vec_out = splineInterp1D(t_vec_in, data_vec_in, t_vec_out, 0);
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


void SeismicForward::findNMOTheta(NRLib::Grid2D<double>     &thetagrid,
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


void SeismicForward::findTWTx(NRLib::Grid2D<double>     &twtx_grid,
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



void SeismicForward::seisConvolutionNMO(NRLib::Grid2D<double>               &timegrid_pos,
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

void SeismicForward::seisConvolution(NRLib::Grid2D<double>               &timegrid_pos,
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




void SeismicForward::generateSeismic(std::vector<NRLib::StormContGrid> &rgridvec,
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
                                     bool timeshift_output)
{
  size_t nt = timegridvec[0].GetNK();
  size_t nc = rgridvec[0].GetNK();
  size_t nz = depthgridvec[0].GetNK();
  std::vector<double> seis(rgridvec.size());
  double x, y, z;
  size_t nx, ny;
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();

  if (time_output == true) {
    nx = timegridvec[0].GetNI();
    ny = timegridvec[0].GetNJ();
  } else {
    nx = depthgridvec[0].GetNI();
    ny = depthgridvec[0].GetNJ();
  }

  printf("\nComputing synthetic seismic:");
  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
  double topt = 0.0;
  //  double start, stop;
  //  start=omp_get_wtime();
  //  #ifdef WITH_OMP
  // int n_threads = omp_get_num_procs();
  // int n_threads = 1;
  //  printf("Antall tr�der: %d \n", n_threads);
  //#pragma omp parallel for private(x,y,z,topt)  num_threads(n_threads)
  //#endif
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      zgrid.FindCenterOfCell(i, j, 0, x, y, z);
      topt = toptime.GetZ(x, y);
      if (time_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double t = t0 + 0.5 * dt;

          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timegridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            t = t + dt;
          }
        }
        else {
          for (size_t k = 0; k < nt; k++){
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timegridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (depth_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double z = z0 + 0.5 * dz;

          std::vector<double> zvec(zgrid.GetNK() + 2);
          std::vector<double> twt(zgrid.GetNK() + 2);

          zvec[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= zgrid.GetNK(); k++) {
            zvec[k] = zgrid(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }
          zvec[zvec.size() - 1] = bot.GetZ(x, y);
          twt[twt.size() - 1] = twt[twt.size() - 2] + 2000.0 * (zvec[zvec.size() - 1] - zvec[zvec.size() - 2]) / constvp[2];


          for (size_t k = 0; k < nz; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(z, zvec, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              depthgridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            z = z + dz;

          }
        } else {
          for (size_t k = 0; k < nz; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              depthgridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (timeshift_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double tt = t0 + 0.5 * dt;

          std::vector<double> twt_shift(twt_timeshift.GetNK() + 1);
          std::vector<double> twt(twtgrid.GetNK() + 1);

          twt_shift[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= twt_timeshift.GetNK(); k++) {
            twt_shift[k] = twt_timeshift(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }

          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(tt, twt_shift, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timeshiftgridvec[l](i, j, k) = static_cast<float>(seis[l]);
            }
            tt = tt + dt;

          }
        }
        else {
          for (size_t k = 0; k < nt; k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              timeshiftgridvec[l](i, j, k) = 0.0;
            }
          }
        }
      }

      if (ny * i + j + 1 >= static_cast<size_t>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }
  std::cout << "\n";
  // stop=omp_get_wtime();
  // double par_time=stop-start;

  // printf("Time: %f \n", par_time);
}


void SeismicForward::generateSeismicOnFile(std::vector<NRLib::StormContGrid> &rgridvec,
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
                                           bool timeshift_output)
{

  //size_t nt = timegridvec[0].GetNK();
  size_t nc = rgridvec[0].GetNK();
  // size_t nz = depthgridvec[0].GetNK();
  std::vector<double> seis(rgridvec.size());
  //double dz;
  double x, y, z;
  //size_t nx, ny;
  //double rickerLimit = 1150.0/wavelet->GetPeakFrequency();
  double rickerLimit = wavelet->GetDepthAdjustmentFactor();
  //  if(time_output == true) {
  //    nx = timegridvec[0].GetNI();
  //    ny = timegridvec[0].GetNJ();
  //  }
  //  else {
  //    nx = depthgridvec[0].GetNI();
  //    ny = depthgridvec[0].GetNJ();
  //  }
  int n_angles = int(rgridvec.size());
  std::ofstream *time_file = new std::ofstream[n_angles];
  std::ofstream *depth_file = new std::ofstream[n_angles];
  std::ofstream *timeshift_file = new std::ofstream[n_angles];
  if (time_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "time_" + NRLib::ToString(i);
      NRLib::OpenWrite(time_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  if (depth_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "depth_" + NRLib::ToString(i);
      // printf("%s\n", filename);
      NRLib::OpenWrite(depth_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  if (timeshift_output == true) {

    for (int i = 0; i < n_angles; i++) {
      std::string filename = "timeshift_" + NRLib::ToString(i);
      NRLib::OpenWrite(timeshift_file[i], filename, std::ios::out | std::ios::binary);
    }
  }
  printf("\nComputing synthetic seismic:");
  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";
  double topt = 0.0;

  for (size_t i = 0; i < static_cast<size_t>(nx); i++) {
    for (size_t j = 0; j < static_cast<size_t>(ny); j++) {
      zgrid.FindCenterOfCell(i, j, 0, x, y, z);
      topt = toptime.GetZ(x, y);
      if (time_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double t = t0 + 0.5 * dt;

          for (size_t k = 0; k < static_cast<size_t>(nt); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(time_file[l], static_cast<float>(seis[l]));
            }
            t = t + dt;
          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nt); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(time_file[l], 0.0);
            }

      }

      if (depth_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double z = z0 + 0.5 * dz;

          std::vector<double> zvec(zgrid.GetNK() + 2);
          std::vector<double> twt(zgrid.GetNK() + 2);


          zvec[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= zgrid.GetNK(); k++) {
            zvec[k] = zgrid(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }
          zvec[zvec.size() - 1] = bot.GetZ(x, y);
          twt[twt.size() - 1] = twt[twt.size() - 2] + 2000.0 * (zvec[zvec.size() - 1] - zvec[zvec.size() - 2]) / constvp[2];


          for (size_t k = 0; k < static_cast<size_t>(nz); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(z, zvec, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(depth_file[l], static_cast<float>(seis[l]));
            }
            z = z + dz;

          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nz); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(depth_file[l], 0.0);
            }
      }

      if (timeshift_output == true) {
        if (toptime.IsMissing(topt) == false) {
          double tt = t0 + 0.5 * dt;

          std::vector<double> twt_shift(twt_timeshift.GetNK() + 1);
          std::vector<double> twt(twtgrid.GetNK() + 1);

          twt_shift[0] = 0.0;
          twt[0] = 0.0;
          for (size_t k = 1; k <= twt_timeshift.GetNK(); k++) {
            twt_shift[k] = twt_timeshift(i, j, k - 1);
            twt[k] = twtgrid(i, j, k - 1);
          }

          for (size_t k = 0; k < static_cast<size_t>(nt); k++) {
            for (size_t l = 0; l < rgridvec.size(); l++) {
              seis[l] = 0.0;
            }
            double t = findTFromZ(tt, twt_shift, twt);
            for (size_t kk = 0; kk < nc; kk++) {
              if (fabs(twtgrid(i, j, kk) - t) < rickerLimit) {
                double ricker = waveletScale * wavelet->FindWaveletPoint(twtgrid(i, j, kk) - t);
                for (size_t l = 0; l < rgridvec.size(); l++) {
                  seis[l] += rgridvec[l](i, j, kk) * ricker;
                }
              }
            }
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(timeshift_file[l], static_cast<float>(seis[l]));
            }
            tt = tt + dt;

          }
        } else
          for (size_t k = 0; k < static_cast<size_t>(nt); k++)
            for (size_t l = 0; l < rgridvec.size(); l++) {
              NRLib::WriteBinaryFloat(timeshift_file[l], 0.0);
            }

      }


      if (ny * i + j + 1 >= static_cast<size_t>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }

    }
  }
  std::cout << "\n";
  // printf("seismic done\n");
  delete[] time_file;
  delete[] depth_file;
  delete[] timeshift_file;
  // stop=omp_get_wtime();
  // double par_time=stop-start;

  // printf("Time: %f \n", par_time);
}

double SeismicForward::findTFromZ(double               z,
                                  std::vector<double> &zvec,
                                  std::vector<double> &tvec)
{
  double t;
  size_t i = 0;
  while (i < zvec.size() - 1 && z > zvec[i]) {
    i++;
  }

  if (i > 0) {
    double a = (zvec[i] - z) / (zvec[i] - zvec[i - 1]);
    t = a * tvec[i - 1] + (1 - a) * tvec[i];
  } else {
    t = tvec[0];
  }

  return t;
}