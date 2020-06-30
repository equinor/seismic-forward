#include "nrlib/eclipsegrid/eclipsegrid.hpp"
#include "nrlib/random/randomgenerator.hpp"
#include "nrlib/random/normal.hpp"
//#include "nrlib/iotools/fileio.hpp"

#include "physics/wavelet.hpp"

#include "tbb/compat/thread"

#include "seismic_regridding.hpp"
#include "tasklist.hpp"

//#include <fstream>

#include <ctime>

#ifdef WITH_OMP
#include <omp.h>
#endif

#include "nrlib/surface/regularsurfacerotated.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/surface/surfaceio.hpp"

//-----------------------------------------------------------------------------------
void SeismicRegridding::MakeSeismicRegridding(SeismicParameters & seismic_parameters,
                                              ModelSettings     * model_settings,
                                              size_t              n_threads)
//-----------------------------------------------------------------------------------
{
  //time_t t1 = time(0);   // get time now
  NRLib::LogKit::WriteHeader("Find depth values");
  FindZValues(seismic_parameters,
              model_settings,
              n_threads);
  //seismic_parameters.PrintElapsedTime(t1, "finding Zvalues");

  //t1 = time(0);
  NRLib::LogKit::WriteHeader("Fill parameter grids");
  FindParameters(seismic_parameters,
                 model_settings,
                 n_threads);

  PostProcess(seismic_parameters,
              model_settings);
  //seismic_parameters.PrintElapsedTime(t1, "finding elastic parameters");

  seismic_parameters.DeleteEclipseGrid();
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nDeleting Eclipse grid to free memory.\n");
  NRLib::LogKit::WriteHeader("Make TWT grid");

  NRLib::RegularSurface<double> & toptime = seismic_parameters.GetTopTime();
  NRLib::RegularSurface<double> & bottime = seismic_parameters.GetBottomTime();
  Wavelet                       * wavelet = seismic_parameters.GetWavelet();

  FindTWT(seismic_parameters, model_settings, toptime, bottime, n_threads);

  NRLib::LogKit::WriteHeader("Export grids");

  //---generate, write and delete vrms grid if writing is requested---------
  if (model_settings->GetNMOCorr() && model_settings->GetOutputVrms()){
    if (model_settings->GetPSSeismic()) {
      NRLib::StormContGrid & twtssgrid = seismic_parameters.GetTwtSSGrid();
      NRLib::StormContGrid & twtppgrid = seismic_parameters.GetTwtPPGrid();
      NRLib::StormContGrid & vpgrid    = seismic_parameters.GetVpGrid();
      NRLib::StormContGrid & vsgrid    = seismic_parameters.GetVsGrid();

      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFind Vrms using Vp.");
      FindVrms(seismic_parameters,
               model_settings,
               vpgrid,
               twtppgrid);
      seismic_parameters.GetSeismicOutput()->WriteVrms(seismic_parameters, "PP");

      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFind Vrms using Vs.");
      FindVrms(seismic_parameters,
               model_settings,
               vsgrid,
               twtssgrid);
      seismic_parameters.GetSeismicOutput()->WriteVrms(seismic_parameters, "SS");
      seismic_parameters.DeleteVrmsGrid();
    }
    else {
      NRLib::StormContGrid & twtgrid = seismic_parameters.GetTwtGrid();
      NRLib::StormContGrid & vpgrid  = seismic_parameters.GetVpGrid();

      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFind Vrms using Vp.");
      FindVrms(seismic_parameters,
               model_settings,
               vpgrid,
               twtgrid);
      seismic_parameters.GetSeismicOutput()->WriteVrms(seismic_parameters);
      seismic_parameters.DeleteVrmsGrid();
    }
  }
  //---add wavelet above and below toptime and bottime-------------
  toptime.Add(-1 * wavelet->GetTwtWavelet()); // add one wavelet length to bot and subtract from top
  bottime.Add(     wavelet->GetTwtWavelet());

  double tmin = toptime.Min();
  double dt   = seismic_parameters.GetSeismicGeometry()->dt();
  size_t ns   = static_cast<size_t>(floor(tmin / dt + 0.5));
  tmin        = ns * dt;
  double tmax = bottime.Max();
  size_t nt   = static_cast<size_t>(floor((tmax - tmin) / dt + 0.5)) + 1;
  seismic_parameters.GetSeismicGeometry()->setNt(nt);
  seismic_parameters.GetSeismicGeometry()->setTRange(tmin, tmax);

  //---write toptime and bottime---------------------------------
  if (model_settings->GetOutputTimeSurfaces()) {
    seismic_parameters.GetSeismicOutput()->WriteTimeSurfaces(seismic_parameters);
  }

  bool interpolate = model_settings->GetResamplParamToSegyInterpol();
  time_t t1 = time(0);   // get time now

  //---resample and write ELASTIC parameters in segy------------
  if (model_settings->GetOutputElasticParametersTimeSegy() || model_settings->GetOutputElasticParametersDepthSegy()) {
    if (interpolate)
      printf("\nStart resampling elastic parameters and write to SegY - with interpolation.");
    else
      printf("\nStart resampling elastic parameters and write to SegY - index based.");
  }
  if (model_settings->GetOutputElasticParametersTimeSegy()) {
    WriteElasticParametersSegy(seismic_parameters, n_threads, true);
  }
  if (model_settings->GetOutputElasticParametersDepthSegy()) {
    WriteElasticParametersSegy(seismic_parameters, n_threads, false);
  }

  //---resample, write and delete EXTRA parameter grids---------
  if (model_settings->GetOutputExtraParametersTimeSegy() || model_settings->GetOutputExtraParametersDepthSegy()) {
    if (interpolate)
      printf("\nStart resampling extra parameters and write to SegY - with interpolation.");
    else
      printf("\nStart resampling extra parameters and write to SegY - index based.");
  }
  if (model_settings->GetOutputExtraParametersTimeSegy()) {
    WriteExtraParametersSegy(seismic_parameters, n_threads, true);
  }
  if (model_settings->GetOutputExtraParametersDepthSegy()) {
    WriteExtraParametersSegy(seismic_parameters, n_threads, false);
  }
  seismic_parameters.DeleteExtraParameterGrids();

  if (model_settings->GetOutputElasticParametersTimeSegy()  ||
      model_settings->GetOutputElasticParametersDepthSegy() ||
      model_settings->GetOutputExtraParametersTimeSegy()    ||
      model_settings->GetOutputExtraParametersDepthSegy()) {
    seismic_parameters.PrintElapsedTime(t1, "resampling parameters and write to SegY.");
  }

  //---write elastic parameters, z values and twt on storm format---
  if (model_settings->GetOutputVp()) {
    seismic_parameters.GetSeismicOutput()->WriteVpVsRho(seismic_parameters);
  }
  if (model_settings->GetOutputZvalues()) {
    seismic_parameters.GetSeismicOutput()->WriteZValues(seismic_parameters);
  }
  if (model_settings->GetOutputTwt()) {
    seismic_parameters.GetSeismicOutput()->WriteTwt(seismic_parameters);
  }
}

//-------------------------------------------------------------------------
void SeismicRegridding::FindZValues(SeismicParameters & seismic_parameters,
                                    ModelSettings     * model_settings,
                                    size_t              n_threads)
//-------------------------------------------------------------------------
{
  NRLib::StormContGrid                & zgrid            = seismic_parameters.GetZGrid();
  const NRLib::EclipseGeometry        & geometry         = seismic_parameters.GetEclipseGrid().GetGeometry();
  const size_t                          top_k            = seismic_parameters.GetTopK();
  const NRLib::RegularSurface<double> & topeclipse       = seismic_parameters.GetTopEclipse();
  const NRLib::RegularSurface<double> & boteclipse       = seismic_parameters.GetBottomEclipse();

  const bool                            rem_neg_delta    = model_settings->GetRemoveNegativeDeltaZ();
  const bool                            use_corner_point = model_settings->GetUseCornerpointInterpol();

  const bool                            bilinear         = false;
  const double                          missing          = -999.0;



  const double                          z_top_shift      = model_settings->GetZWaveletTop();
  const double                          z_bot_shift      = model_settings->GetZWaveletBot();

  NRLib::RegularSurface<double>         orig_top(topeclipse);
  NRLib::RegularSurface<double>         orig_bot(boteclipse);

  //orig_top.Add(     z_top_shift);
  //orig_bot.Add(-1 * z_bot_shift);

  geometry.FindRegularGridOfZValues(zgrid,
                                    orig_top,
                                    orig_bot,
                                    top_k,
                                    n_threads,
                                    use_corner_point,
                                    bilinear,
                                    missing);

  //
  // There should be no more negative dz values anymore, but we check anyway ...
  //
  RemoveNegativeDz(zgrid,
                   n_threads,
                   rem_neg_delta,
                   missing);
}

//----------------------------------------------------------------------------
void SeismicRegridding::RemoveNegativeDz(NRLib::StormContGrid & zgrid,
                                         const size_t           n_threads,
                                         const bool             rem_neg_delta,
                                         const double           missing)
//----------------------------------------------------------------------------
{
  const size_t ni = zgrid.GetNI();
  const size_t nj = zgrid.GetNJ();
  const size_t nk = zgrid.GetNK();

  std::vector<std::vector<double> > negative_dz_pts;

  double dz_min = 0.0;
  int    min_i  = -1;
  int    min_j  = -1;
  int    min_k  = -1;

#ifdef WITH_OMP
  int chunk_size = 1;
  chunk_size = std::max(1, int(ni/(2*static_cast<int>(n_threads))));
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
  for (size_t i = 0; i < ni ; i++) {
    std::vector<std::vector<double> > neg_dz_pts;
    for (size_t j = 0; j < nj; j++) {
      for (int k = static_cast<int>(nk - 2); k >= 0; --k) {
        size_t kk = static_cast<size_t>(k);
        float  z1 = zgrid(i, j, kk    );
        float  z2 = zgrid(i, j, kk + 1);

        if (z1 != missing && z2 != missing && z1 > z2) {
          if (rem_neg_delta) {
            zgrid(i, j, kk) = z2;
          }
          double x, y, z;
          zgrid.FindCenterOfCell(i, j, kk, x, y, z);

          // Logging negative values
          std::vector<double> neg(5);
          neg[0] = static_cast<double>(k);
          neg[1] = x;
          neg[2] = y;
          neg[3] = zgrid(i, j, k);
          neg[4] = z2 - z1;
          neg_dz_pts.push_back(neg);
          if (z2 - z1 < dz_min) {
            dz_min = z2 - z1;
            min_i  = i;
            min_j  = j;
            min_k  = k;
          }
        }
      }
    }
#ifdef WITH_OMP
#pragma omp critical
#endif
    {
      negative_dz_pts.reserve(negative_dz_pts.size() + neg_dz_pts.size());
      negative_dz_pts.insert(negative_dz_pts.end(), neg_dz_pts.begin(), neg_dz_pts.end());
    }
  }

  size_t n = negative_dz_pts.size();

  if (n > 0) {
    if (rem_neg_delta)
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nNumber of negative dz found and removed   : %5d", n);
    else {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nNumber of negative dz found               : %5d", n);
      TaskList::AddTask("Check section 3: Negative dz values found when building z-grid");
    }
    double x,y,z;
    zgrid.FindCenterOfCell(min_i, min_j, min_k, x, y, z);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nLargest negative value                    : %5.2f\n", dz_min);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nLargest negative value location (x,y,z)   : (%.2f, %.2f, %.2f)\n", x, y, z);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nLargest negative value location (i,j,k)   : (%d, %d, %d)\n", min_i, min_j, min_k);
  }
  else {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nNo crossing depth values found!\n");
  }

  bool debug_neg_dz = true;
  if (n > 0 && debug_neg_dz) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nExporting points where there are negative dz values (z-value is the conflicting value)!\n");
    std::fstream fout;
    NRLib::OpenWrite(fout, "negative_dz_points.rxat");
    fout << "Float Negative dz\n"
         << "Discrete Layer"
         << std::endl;
    for (size_t i = 0; i < negative_dz_pts.size(); i++) {
      fout << std::fixed
        << std::setprecision(2)
        << std::setw(12) << negative_dz_pts[i][1] << " "
        << std::setw(12) << negative_dz_pts[i][2] << " "
        << std::setw(8)  << negative_dz_pts[i][3] << " "
        << std::setw(8)  << negative_dz_pts[i][4] << " "
        << std::setw(8)  << static_cast<int>(negative_dz_pts[i][0])
           << std::endl;
    }
    fout.close();
  }
}


//----------------------------------------------------------------------
void SeismicRegridding::FindParameters(SeismicParameters & seismic_parameters,
                                       ModelSettings     * model_settings,
                                       size_t              n_threads)
//----------------------------------------------------------------------
{
  NRLib::StormContGrid               & vpgrid                   = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid               & vsgrid                   = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid               & rhogrid                  = seismic_parameters.GetRhoGrid();
  std::vector<NRLib::StormContGrid*>   extra_parameter_grid     = seismic_parameters.GetExtraParametersGrids();

  const NRLib::EclipseGrid           & egrid                    = seismic_parameters.GetEclipseGrid();

  size_t                               topk                     = seismic_parameters.GetTopK();
  size_t                               botk                     = seismic_parameters.GetBottomK();

  double                               zlimit                   = model_settings->GetZeroThicknessLimit();
  std::vector<double>                  constvp                  = model_settings->GetConstVp();
  std::vector<double>                  constvs                  = model_settings->GetConstVs();
  std::vector<double>                  constrho                 = model_settings->GetConstRho();
  std::vector<std::string>             names                    = model_settings->GetParameterNames();
  std::vector<double>                  extra_parameter_defaults = model_settings->GetExtraParameterDefaultValues();
  std::vector<std::string>             extra_parameter_names;

  // Only resample extra parameters if requested for output segy.
  if (model_settings->GetOutputExtraParametersTimeSegy() || model_settings->GetOutputExtraParametersDepthSegy()) {
    extra_parameter_names = model_settings->GetExtraParameterNames();
  }
  size_t n_extra_params = extra_parameter_names.size();

  const NRLib::EclipseGeometry       & eclipse_geometry         = egrid.GetGeometry();
  //---for parallelisation
  //use copy-constructor, need copy as values are filled in.
  NRLib::Grid<double>                  eclipse_vp               = egrid.GetParameter(names[0]);
  NRLib::Grid<double>                  eclipse_vs               = egrid.GetParameter(names[1]);
  NRLib::Grid<double>                  eclipse_rho              = egrid.GetParameter(names[2]);

  std::vector<NRLib::Grid<double> >    eclipse_extra_params;

  for (size_t i = 0; i < n_extra_params; ++i) {
    NRLib::Grid<double> one_parameter_grid = egrid.GetParameter(extra_parameter_names[i]);
    eclipse_extra_params.push_back(one_parameter_grid);
  }

  size_t nijk = egrid.GetNI()*egrid.GetNJ()*egrid.GetNK();
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFilling inactive cells in Eclipse grid above and in reservoir.\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nTotal number of grid cells: %d x %d x %d = %d\n",egrid.GetNI(),egrid.GetNJ(),egrid.GetNK(), nijk);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nOB1 = Using default value for overburden above first layer.");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nOB2 = Using default value for overburden in and below first layer.");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nRES = Using default value for reservoir in reservoir.");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nZRO = Using cell value from cell above for zero thickness cells.\n");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nGridName            OB1     OB2     RES     ZRO");
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n-----------------------------------------------\n");

  //-----prepare eclipsegrid - include default values and value above where delta < zlimit
  FillInGridValues("Vp" , eclipse_geometry, eclipse_vp , constvp[0] , constvp[1] , zlimit, topk, botk);
  FillInGridValues("Vs" , eclipse_geometry, eclipse_vs , constvs[0] , constvs[1] , zlimit, topk, botk);
  FillInGridValues("Rho", eclipse_geometry, eclipse_rho, constrho[0], constrho[1], zlimit, topk, botk);
  for (size_t ii = 0; ii < n_extra_params; ++ii) {
    FillInGridValues(extra_parameter_names[ii], eclipse_geometry, eclipse_extra_params[ii], extra_parameter_defaults[ii], extra_parameter_defaults[ii], zlimit, topk, botk);
  }

  //default value in top
  for (size_t i = 0; i < vpgrid.GetNI(); i++) {
    for (size_t j = 0; j < vpgrid.GetNJ(); j++) {
      vpgrid (i, j, 0) = static_cast<float>(constvp[0]);
      vsgrid (i, j, 0) = static_cast<float>(constvs[0]);
      rhogrid(i, j, 0) = static_cast<float>(constrho[0]);
      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        (*extra_parameter_grid[ii])(i, j, 0) = 0.0;
      }
    }
  }

  //blocking - for parallelisation - if n_threads > 1
  size_t nbx = 1;
  size_t nby = 1;
  size_t nb  = 1;
  size_t nx  = egrid.GetNI() - 1;
  size_t ny  = egrid.GetNJ() - 1;
  size_t nxb = nx;
  size_t nyb = ny;
  if (n_threads > 1){
    nbx = 10;
    nby = 10;
    nb  = nbx*nby;
    nxb = static_cast<size_t>(floor(static_cast<double>(nx) / static_cast<double>(nbx) + 0.5));
    nyb = static_cast<size_t>(floor(static_cast<double>(ny) / static_cast<double>(nby) + 0.5));
  }

  double angle     = vpgrid.GetAngle();
  double cosA      = cos(angle);
  double sinA      = sin(angle);
  double x_min_rot = vpgrid.GetXMin()*cosA + vpgrid.GetYMin()*sinA;
  double y_min_rot = vpgrid.GetYMin()*cosA - vpgrid.GetXMin()*sinA;

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nResampling parameters in Eclipse grid into regular grids.\n");

#ifdef WITH_OMP
  int  chunk_size = 1;
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
  for (size_t b = 0 ; b < nb ; ++b) {

    size_t ibx  = b % nbx;
    size_t iby  = static_cast<size_t>(std::floor(static_cast<double>(b) / static_cast<double>(nbx)));
    size_t imin = std::max(ibx*nxb, static_cast<size_t>(0)); // find min and max of block
    size_t jmin = std::max(iby*nyb, static_cast<size_t>(0));
    size_t imax = nx;
    size_t jmax = ny;

    if (ibx != (nbx - 1))
      imax = std::min((ibx + 1)*nxb, nx);
    if (iby != (nby - 1))
      jmax = std::min((iby + 1)*nyb, ny);

    std::vector<double>       x_rot(4), y_rot(4);
    std::vector<bool>         inside(4);
    std::vector<NRLib::Point> pt_vp(4), pt_vs(4), pt_rho(4);

    std::vector<std::vector<NRLib::Point> > pt_extra_param(n_extra_params);
    for (size_t i = 0; i < n_extra_params; ++i)
      pt_extra_param[i] = pt_vp;

    for (size_t k = topk; k <= botk + 1; k++) {
      for (size_t i = imin; i < imax; ++i) {
        for (size_t j = jmin; j < jmax; ++j) {

          if (eclipse_geometry.IsPillarActive(i    , j    ) &&
              eclipse_geometry.IsPillarActive(i + 1, j    ) &&
              eclipse_geometry.IsPillarActive(i    , j + 1) &&
              eclipse_geometry.IsPillarActive(i + 1, j + 1) &&
              eclipse_geometry.IsPillarActive(i + 2, j    ) &&
              eclipse_geometry.IsPillarActive(i + 2, j + 1) &&
              eclipse_geometry.IsPillarActive(i    , j + 2) &&
              eclipse_geometry.IsPillarActive(i + 1, j + 2) &&
              eclipse_geometry.IsPillarActive(i + 2, j + 2)) {

            if (k <= botk) {
              for (size_t pt = 0; pt < 4; ++pt)
                pt_vp[pt] = eclipse_geometry.FindCellCenterPoint(i + int(pt % 2), j + int(floor(double(pt) / 2)), k);
            }
            else {
              for (size_t pt = 0; pt < 4; ++pt)
                pt_vp[pt] = eclipse_geometry.FindCellCenterPoint(i + int(pt % 2), j + int(floor(double(pt) / 2)), k - 1);
            }

            for (size_t pt = 0; pt < 4; ++pt)
              inside[pt] = vpgrid.IsInside(pt_vp[pt].x, pt_vp[pt].y);

            if (inside[0] || inside[1] || inside[2] || inside[3]) {

              for (size_t pt = 0; pt < 4; ++pt) {
                pt_vs[pt]  = pt_vp[pt];
                pt_rho[pt] = pt_vp[pt];
                for (size_t ii = 0; ii < n_extra_params; ++ii) {
                  pt_extra_param[ii][pt] = pt_vp[pt];
                }
              }

              if (k == botk + 1) {
                for (size_t pt = 0; pt < 4; ++pt) {
                  pt_vp[pt].z  = constvp[2];
                  pt_vs[pt].z  = constvs[2];
                  pt_rho[pt].z = constrho[2];
                  for (size_t ii = 0; ii < n_extra_params; ++ii) {
                    pt_extra_param[ii][pt].z = 0.0;
                  }
                }
              }
              else {
                for (size_t pt = 0; pt < 4; ++pt) {
                  pt_vp[pt].z  = eclipse_vp (i + int(pt % 2), j + int(floor(double(pt / 2))), k);
                  pt_vs[pt].z  = eclipse_vs (i + int(pt % 2), j + int(floor(double(pt / 2))), k);
                  pt_rho[pt].z = eclipse_rho(i + int(pt % 2), j + int(floor(double(pt / 2))), k);
                  for (size_t ii = 0; ii < n_extra_params; ++ii) {
                    pt_extra_param[ii][pt].z = eclipse_extra_params[ii](i + int(pt % 2), j + int(floor(double(pt / 2))), k);
                  }
                }
              }

              bool triangulate_124 = Is124Triangulate(pt_vp);

              std::vector<NRLib::Triangle> triangles_elastic(6);
              std::vector<NRLib::Triangle> triangles_extra_param(n_extra_params * 2);

              SetElasticTriangles(pt_vp, pt_vs, pt_rho, pt_extra_param, triangulate_124, triangles_elastic, triangles_extra_param);

              for (size_t pt = 0; pt < 4; ++pt) {
                x_rot[pt] = pt_vp[pt].x * cosA + pt_vp[pt].y *sinA;
                y_rot[pt] = pt_vp[pt].y * cosA - pt_vp[pt].x *sinA;
              }

              double cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
              double cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
              double cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
              double cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

              size_t start_ii = static_cast<size_t>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 0.5));
              size_t start_jj = static_cast<size_t>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 0.5));
              size_t end_ii   = static_cast<size_t>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 1.0));
              size_t end_jj   = static_cast<size_t>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 1.0));

              if (end_ii > vpgrid.GetNI()) {
                end_ii = vpgrid.GetNI();
              }
              if (end_jj > vpgrid.GetNJ()) {
                end_jj = vpgrid.GetNJ();
              }

              for (size_t ii = start_ii; ii < end_ii; ii++) {
                for (size_t jj = start_jj; jj < end_jj; jj++) {

                  double x, y, z;
                  vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);

                  NRLib::Point p1, p2;
                  p1.x  = x;
                  p1.y  = y;
                  p1.z  = pt_vp[0].z;
                  p2    = p1;
                  p2.z += 1000;

                  NRLib::Line  line(p1, p2, false, false);
                  NRLib::Point intersec_pt;

                  if (triangles_elastic[0].FindNearestPoint(line, intersec_pt) < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    triangles_elastic[2].FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    triangles_elastic[4].FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    for (size_t iii = 0; iii < n_extra_params; ++iii) {
                      triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                      NRLib::StormContGrid &param_grid = *(extra_parameter_grid[iii]);
                      param_grid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                  else if (triangles_elastic[1].FindNearestPoint(line, intersec_pt) < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    triangles_elastic[3].FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    triangles_elastic[5].FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    for (size_t iii = 0; iii < n_extra_params; ++iii) {
                      triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                      NRLib::StormContGrid &param_grid = *(extra_parameter_grid[iii]);
                      param_grid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  ////-------------find edges---------------------
  for (size_t k = topk; k <= botk + 1; k++) {
    for (size_t i = 0; i < egrid.GetNI() - 1; i++) {
      //bot edge
      size_t j = 0;
      if (FindBotCell(eclipse_geometry, egrid.GetNJ(), i, j)){
        FindEdges(seismic_parameters,
                  eclipse_geometry,
                  eclipse_vp,
                  eclipse_vs,
                  eclipse_rho,
                  eclipse_extra_params,
                  i, j, k,
                  false, true, false, false);
      }
      //top edge
      j = egrid.GetNJ() - 1;
      if (FindTopCell(eclipse_geometry, i, j)) {
        FindEdges(seismic_parameters,
                  eclipse_geometry,
                  eclipse_vp,
                  eclipse_vs,
                  eclipse_rho,
                  eclipse_extra_params,
                  i, j, k,
                  true, false, false, false);
      }
    }
    for (size_t j = 0; j < egrid.GetNJ() - 1; ++j) {
      //left edge
      size_t i = 0;
      if (FindLeftCell(eclipse_geometry, egrid.GetNI(), i, j)) {
        FindEdges(seismic_parameters,
                  eclipse_geometry,
                  eclipse_vp,
                  eclipse_vs,
                  eclipse_rho,
                  eclipse_extra_params,
                  i, j, k,
                  false, false, false, true);
      }
      //right edge
      i = egrid.GetNI() - 1;
      if (FindRightCell(eclipse_geometry, i, j)) {
        FindEdges(seismic_parameters,
                  eclipse_geometry,
                  eclipse_vp,
                  eclipse_vs,
                  eclipse_rho,
                  eclipse_extra_params,
                  i, j, k,
                  false, false, true, false);
      }
    }
    //-------------find corners---------------------
    //bot left
    size_t i = 0;
    size_t j = 0;
    std::vector<NRLib::Point> pt_vp(4);
    FindCornerCellPoints(eclipse_geometry,
                         pt_vp,
                         i, j, k,
                         botk);
    FindCorners(seismic_parameters,
                eclipse_geometry,
                eclipse_vp,
                eclipse_vs,
                eclipse_rho,
                eclipse_extra_params,
                i, j, k, pt_vp);
    //top left
    j = egrid.GetNJ() - 1;
    FindCornerCellPoints(eclipse_geometry,
                         pt_vp,
                         i, j, k,
                         botk);
    FindCorners(seismic_parameters,
                eclipse_geometry,
                eclipse_vp,
                eclipse_vs,
                eclipse_rho,
                eclipse_extra_params,
                i, j, k, pt_vp);
    //top right
    i = egrid.GetNI() - 1;
    FindCornerCellPoints(eclipse_geometry,
                         pt_vp,
                         i, j, k,
                         botk);
    FindCorners(seismic_parameters,
                eclipse_geometry,
                eclipse_vp,
                eclipse_vs,
                eclipse_rho,
                eclipse_extra_params,
                i, j, k, pt_vp);
    //bot right
    j = 0;
    FindCornerCellPoints(eclipse_geometry,
                         pt_vp,
                         i, j, k,
                         botk);
    FindCorners(seismic_parameters,
                eclipse_geometry,
                eclipse_vp,
                eclipse_vs,
                eclipse_rho,
                eclipse_extra_params,
                i, j, k, pt_vp);
  }
}

//-------------------------------------------------------------------------------------
void SeismicRegridding::FillInGridValues(const std::string            & text,
                                         const NRLib::EclipseGeometry & geometry,
                                         NRLib::Grid<double>          & grid_copy,
                                         double                         default_top,    // default value above
                                         double                         default_value,  // default value inside
                                         double                         zlimit,         // zero thickness limit
                                         size_t                         topk,
                                         size_t                         botk)
//-------------------------------------------------------------------------------------
{
  int nzlimit = 0;
  int ndeftop = 0;
  int ndefins = 0;
  int ndef    = 0;

  for (size_t k = topk ; k <= botk ; k++) {
    for (size_t i = 0; i < geometry.GetNI(); i++) {
      for (size_t j = 0; j < geometry.GetNJ(); j++) {
        if (!geometry.IsActive(i, j, k)) {
          if (k > 0 && k > topk) {
            if (geometry.GetDZ(i, j, k) < zlimit) {
              grid_copy(i, j, k) = grid_copy(i, j, k - 1);
              nzlimit++;
            }
            else if (grid_copy(i, j, k - 1) == default_top) {
              grid_copy(i, j, k) = default_top;
              ndefins++;
            }
            else {
              grid_copy(i, j, k) = default_value;
              ndef++;
            }
          }
          else {
            grid_copy(i, j, k) = default_top;
            ndeftop++;
          }
        }
      }
    }
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "%-15s %7d %7d %7d %7d\n",text.c_str(), ndeftop, ndefins, ndef, nzlimit);
}

//-----------------------------------------------------------------------
bool SeismicRegridding::Is124Triangulate(std::vector<NRLib::Point> pt_vp)
//-----------------------------------------------------------------------
{
  bool triangulate_124 = true;
  NRLib::Point vec1, vec2;
  vec1   = pt_vp[0] - pt_vp[1];
  vec1.z = 0;
  vec2   = pt_vp[3] - pt_vp[1];
  vec2.z = 0;
  double delaunay_angle = vec1.GetAngle(vec2);
  vec1   = pt_vp[0] - pt_vp[2];
  vec1.z = 0;
  vec2   = pt_vp[3] - pt_vp[2];
  vec2.z = 0;
  delaunay_angle += vec1.GetAngle(vec2);
  if (delaunay_angle > NRLib::Pi) {
    triangulate_124 = false;
  }
  return triangulate_124;
}

//----------------------------------------------------------------------------------------------------------
void SeismicRegridding::SetElasticTriangles(std::vector<NRLib::Point>               & pt_vp,
                                            std::vector<NRLib::Point>               & pt_vs,
                                            std::vector<NRLib::Point>               & pt_rho,
                                            std::vector<std::vector<NRLib::Point> > & pt_extra_param,
                                            bool                                      triangulate_124,
                                            std::vector<NRLib::Triangle>            & triangles_elastic,
                                            std::vector<NRLib::Triangle>            & triangles_extra_param)
//----------------------------------------------------------------------------------------------------------
{
  if (triangulate_124) {
    triangles_elastic[0].SetCornerPoints(pt_vp [0], pt_vp [1], pt_vp [3]);
    triangles_elastic[1].SetCornerPoints(pt_vp [0], pt_vp [2], pt_vp [3]);
    triangles_elastic[2].SetCornerPoints(pt_vs [0], pt_vs [1], pt_vs [3]);
    triangles_elastic[3].SetCornerPoints(pt_vs [0], pt_vs [2], pt_vs [3]);
    triangles_elastic[4].SetCornerPoints(pt_rho[0], pt_rho[1], pt_rho[3]);
    triangles_elastic[5].SetCornerPoints(pt_rho[0], pt_rho[2], pt_rho[3]);
    for (size_t ii = 0; ii < triangles_extra_param.size()/2; ++ii){
      triangles_extra_param[ii*2    ].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][1], pt_extra_param[ii][3]);
      triangles_extra_param[ii*2 + 1].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][2], pt_extra_param[ii][3]);
    }
  }
  else {
    triangles_elastic[0].SetCornerPoints(pt_vp [0], pt_vp [1], pt_vp [2]);
    triangles_elastic[1].SetCornerPoints(pt_vp [1], pt_vp [2], pt_vp [3]);
    triangles_elastic[2].SetCornerPoints(pt_vs [0], pt_vs [1], pt_vs [2]);
    triangles_elastic[3].SetCornerPoints(pt_vs [1], pt_vs [2], pt_vs [3]);
    triangles_elastic[4].SetCornerPoints(pt_rho[0], pt_rho[1], pt_rho[2]);
    triangles_elastic[5].SetCornerPoints(pt_rho[1], pt_rho[2], pt_rho[3]);
    for (size_t ii = 0; ii < triangles_extra_param.size()/2; ++ii){
      triangles_extra_param[ii*2    ].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][1], pt_extra_param[ii][2]);
      triangles_extra_param[ii*2 + 1].SetCornerPoints(pt_extra_param[ii][1], pt_extra_param[ii][2], pt_extra_param[ii][3]);
    }
  }
}

//---------------------------------------------------------------------------
bool SeismicRegridding::FindTopCell(const NRLib::EclipseGeometry & geometry,
                                    size_t                         i,
                                    size_t                       & jj)
//---------------------------------------------------------------------------
{
  int j = static_cast<int>(jj);
  while (j >= 0 && !(geometry.IsPillarActive(i    ,     j) &&
                     geometry.IsPillarActive(i + 1,     j) &&
                     geometry.IsPillarActive(i    , j + 1) &&
                     geometry.IsPillarActive(i + 1, j + 1) &&
                     geometry.IsPillarActive(i + 2,     j) &&
                     geometry.IsPillarActive(i + 2, j + 1))) {
    j--;
  }
  if (j >= 0){
    jj = static_cast<size_t>(j);
    return true;
  }
  else
    return false;
}

//---------------------------------------------------------------------------
bool SeismicRegridding::FindBotCell(const NRLib::EclipseGeometry & geometry,
                                    size_t                         nj,
                                    size_t                         i,
                                    size_t                       & j)
//---------------------------------------------------------------------------
{
  while (j < nj && !(geometry.IsPillarActive(i    ,     j) &&
                     geometry.IsPillarActive(i + 1,     j) &&
                     geometry.IsPillarActive(i    , j + 1) &&
                     geometry.IsPillarActive(i + 1, j + 1) &&
                     geometry.IsPillarActive(i + 2,     j) &&
                     geometry.IsPillarActive(i + 2, j + 1))) {
    j++;
  }
  if (j < nj)
    return true;
  else
    return false;
}

//---------------------------------------------------------------------------
bool SeismicRegridding::FindLeftCell(const NRLib::EclipseGeometry & geometry,
                                     size_t                         ni,
                                     size_t                       & i,
                                     size_t                         j)
//---------------------------------------------------------------------------
{
  while (i < ni && !(geometry.IsPillarActive(i    ,     j) &&
                     geometry.IsPillarActive(i    , j + 1) &&
                     geometry.IsPillarActive(i + 1,     j) &&
                     geometry.IsPillarActive(i + 1, j + 1) &&
                     geometry.IsPillarActive(i    , j + 2) &&
                     geometry.IsPillarActive(i + 1, j + 2))) {
    i++;
  }
  if (i < ni)
    return true;
  else
    return false;
}

//---------------------------------------------------------------------------
bool SeismicRegridding::FindRightCell(const NRLib::EclipseGeometry & geometry,
                                      size_t                       & ii,
                                      size_t                         j)
//---------------------------------------------------------------------------
{
  int i = static_cast<int>(ii);
  while (i >= 0 && !(geometry.IsPillarActive(i    ,     j) &&
                     geometry.IsPillarActive(i    , j + 1) &&
                     geometry.IsPillarActive(i + 1,     j) &&
                     geometry.IsPillarActive(i + 1, j + 1) &&
                     geometry.IsPillarActive(i    , j + 2) &&
                     geometry.IsPillarActive(i + 1, j + 2))) {
    i--;
  }
  if (i >= 0){
    ii = static_cast<size_t>(i);
    return true;
  }
  else
    return false;
}

//-------------------------------------------------------------------------------------------
void SeismicRegridding::FindEdges(SeismicParameters                   & seismic_parameters,
                                  const NRLib::EclipseGeometry        & eclipse_geometry,
                                  const NRLib::Grid<double>           & eclipse_vp,
                                  const NRLib::Grid<double>           & eclipse_vs,
                                  const NRLib::Grid<double>           & eclipse_rho,
                                  std::vector<NRLib::Grid<double> >   & eclipse_extra_params,
                                  size_t                                i,
                                  size_t                                j,
                                  size_t                                k,
                                  bool                                  top,
                                  bool                                  bot,
                                  bool                                  right,
                                  bool                                  left)
//-------------------------------------------------------------------------------------------
{
  NRLib::StormContGrid                  & vpgrid                         = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid                  & vsgrid                         = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid                  & rhogrid                        = seismic_parameters.GetRhoGrid();
  std::vector<NRLib::StormContGrid*>      extra_parameter_grid           = seismic_parameters.GetExtraParametersGrids();

  std::vector<double>                     constvp                        = seismic_parameters.GetModelSettings()->GetConstVp();
  std::vector<double>                     constvs                        = seismic_parameters.GetModelSettings()->GetConstVs();
  std::vector<double>                     constrho                       = seismic_parameters.GetModelSettings()->GetConstRho();
  std::vector<double>                     extra_parameter_default_values = seismic_parameters.GetModelSettings()->GetExtraParameterDefaultValues();

  size_t                                  n_extra_params                 = eclipse_extra_params.size();

  size_t                                  topk                           = seismic_parameters.GetTopK();
  size_t                                  botk                           = seismic_parameters.GetBottomK();

  double                                  angle                          = vpgrid.GetAngle();
  double                                  cosA                           = cos(angle);
  double                                  sinA                           = sin(angle);
  double                                  x_min_rot                      = vpgrid.GetXMin()*cosA + vpgrid.GetYMin()*sinA;
  double                                  y_min_rot                      = vpgrid.GetYMin()*cosA - vpgrid.GetXMin()*sinA;

  NRLib::Point                            mid_edge1;
  NRLib::Point                            mid_edge2;

  std::vector<double>                     x_rot(4);
  std::vector<double>                     y_rot(4);
  std::vector<bool>                       inside(4);

  std::vector<NRLib::Point>               pt_vp(4);
  std::vector<NRLib::Point>               pt_vs(4);
  std::vector<NRLib::Point>               pt_rho(4);

  std::vector<std::vector<NRLib::Point> > pt_extra_param(n_extra_params);
  for (size_t ii = 0; ii < n_extra_params; ++ii)
    pt_extra_param[ii] = pt_vp;

  std::vector<size_t> a_corn(4), b_corn(4), c_corn(4);
  GetCornerPointDir(a_corn, b_corn, c_corn, left, right, bot, top);

  size_t ic = i;
  size_t jc = j;
  if (bot || top)
    ic = i+1;
  else if (left || right)
    jc = j+1;

  if (k <= botk) {
    pt_vp[0] = eclipse_geometry.FindCellCenterPoint(i,  j,  k);
    pt_vp[1] = eclipse_geometry.FindCellCenterPoint(ic, jc, k);

    mid_edge1 = 0.5*(eclipse_geometry.FindCornerPoint(i,  j,  k, a_corn[0], b_corn[0], c_corn[0])
                   + eclipse_geometry.FindCornerPoint(i,  j,  k, a_corn[1], b_corn[1], c_corn[1]))
              + 0.5*(eclipse_geometry.FindCornerPoint(i,  j,  k, a_corn[2], b_corn[2], c_corn[2])
                   + eclipse_geometry.FindCornerPoint(i,  j,  k, a_corn[3], b_corn[3], c_corn[3]));

    mid_edge2 = 0.5*(eclipse_geometry.FindCornerPoint(ic, jc, k, a_corn[0], b_corn[0], c_corn[0])
                   + eclipse_geometry.FindCornerPoint(ic, jc, k, a_corn[1], b_corn[1], c_corn[1]))
              + 0.5*(eclipse_geometry.FindCornerPoint(ic, jc, k, a_corn[2], b_corn[2], c_corn[2])
                   + eclipse_geometry.FindCornerPoint(ic, jc, k, a_corn[3], b_corn[3], c_corn[3]));
  }
  else {
    pt_vp[0] = eclipse_geometry.FindCellCenterPoint(i,  j,  k - 1);
    pt_vp[1] = eclipse_geometry.FindCellCenterPoint(ic, jc, k - 1);

    mid_edge1 = 0.5*(eclipse_geometry.FindCornerPoint(i,  j,  k - 1, a_corn[0], b_corn[0], c_corn[0])
                   + eclipse_geometry.FindCornerPoint(i,  j,  k - 1, a_corn[1], b_corn[1], c_corn[1]))
              + 0.5*(eclipse_geometry.FindCornerPoint(i,  j,  k - 1, a_corn[2], b_corn[2], c_corn[2])
                   + eclipse_geometry.FindCornerPoint(i,  j,  k - 1, a_corn[3], b_corn[3], c_corn[3]));

    mid_edge2 = 0.5*(eclipse_geometry.FindCornerPoint(ic, jc, k - 1, a_corn[0], b_corn[0], c_corn[0])
                   + eclipse_geometry.FindCornerPoint(ic, jc, k - 1, a_corn[1], b_corn[1], c_corn[1]))
              + 0.5*(eclipse_geometry.FindCornerPoint(ic, jc, k - 1, a_corn[2], b_corn[2], c_corn[2])
                   + eclipse_geometry.FindCornerPoint(ic, jc, k - 1, a_corn[3], b_corn[3], c_corn[3]));
  }

  pt_vp[2]  = mid_edge1 - pt_vp[0];
  mid_edge1 = 0.5 * mid_edge1;
  pt_vp[3]  = mid_edge2 - pt_vp[1];
  for (size_t pt = 0; pt < 4;++pt)
    inside[pt] = vpgrid.IsInside(pt_vp[pt].x, pt_vp[pt].y);
  mid_edge2 = 0.5 * mid_edge2;

  if (inside[0] || inside[1] || inside[2] || inside[3]) {
    for (size_t pt = 0; pt < 4; ++pt){
      pt_vs[pt]  = pt_vp[pt];
      pt_rho[pt] = pt_vp[pt];
      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        pt_extra_param[ii][pt] = pt_vp[pt];
      }
    }
    if (k == botk + 1) {
      for (size_t pt = 0; pt < 2; ++pt) { //nb, only loop two first points here
        pt_vp[pt].z  = constvp[2];
        pt_vs[pt].z  = constvs[2];
        pt_rho[pt].z = constrho[2];
        for (size_t ii = 0; ii < n_extra_params; ++ii) {
          pt_extra_param[ii][pt].z = 0.0;
        }
      }
    }
    else {
      pt_vp[0].z  = eclipse_vp (i , j , k);
      pt_vs[0].z  = eclipse_vs (i , j , k);
      pt_rho[0].z = eclipse_rho(i , j , k);
      pt_vp[1].z  = eclipse_vp (ic, jc, k);
      pt_vs[1].z  = eclipse_vs (ic, jc, k);
      pt_rho[1].z = eclipse_rho(ic, jc, k);

      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        pt_extra_param[ii][0].z = eclipse_extra_params[ii](i,  j,  k);
        pt_extra_param[ii][1].z = eclipse_extra_params[ii](ic, jc, k);
      }
    }
    for (size_t pt = 2; pt < 4; ++pt) { //nb, only loop two last points here
      pt_vp[pt].z  = pt_vp[pt-2].z;
      pt_vs[pt].z  = pt_vs[pt-2].z;
      pt_rho[pt].z = pt_rho[pt-2].z;
      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        pt_extra_param[ii][pt].z = pt_extra_param[ii][pt-2].z;
      }
    }

    bool triangulate_124 = Is124Triangulate(pt_vp);

    std::vector<NRLib::Triangle> triangles_elastic(6);
    std::vector<NRLib::Triangle> triangles_extra_param(n_extra_params*2);
    SetElasticTriangles(pt_vp, pt_vs, pt_rho, pt_extra_param, triangulate_124, triangles_elastic, triangles_extra_param);

    for (size_t pt = 0; pt < 4; ++pt) {
      x_rot[pt] = pt_vp[pt].x * cosA + pt_vp[pt].y *sinA;
      y_rot[pt] = pt_vp[pt].y * cosA - pt_vp[pt].x *sinA;
    }

    double cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
    double cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
    double cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
    double cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

    size_t start_ii = static_cast<size_t>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
    size_t start_jj = static_cast<size_t>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
    size_t end_ii   = static_cast<size_t>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
    size_t end_jj   = static_cast<size_t>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));

    if (end_ii > vpgrid.GetNI()) {
      end_ii = vpgrid.GetNI();
    }
    if (end_jj > vpgrid.GetNJ()) {
      end_jj = vpgrid.GetNJ();
    }

    NRLib::Polygon inside_e_cells;
    inside_e_cells.AddPoint(pt_vp[0]);
    inside_e_cells.AddPoint(pt_vp[1]);
    inside_e_cells.AddPoint(mid_edge2);
    if (k <= botk) {
      inside_e_cells.AddPoint(0.5*(eclipse_geometry.FindCornerPoint(i, j, k,     a_corn[2], b_corn[2], c_corn[2])
                                 + eclipse_geometry.FindCornerPoint(i, j, k,     a_corn[3], b_corn[3], c_corn[3])));
    }
    else {
      inside_e_cells.AddPoint(0.5*(eclipse_geometry.FindCornerPoint(i, j, k - 1, a_corn[2], b_corn[2], c_corn[2])
                                 + eclipse_geometry.FindCornerPoint(i, j, k - 1, a_corn[3], b_corn[3], c_corn[3])));
    }

    inside_e_cells.AddPoint(mid_edge1);
    for (size_t ii = start_ii; ii < end_ii; ii++) {
      for (size_t jj = start_jj; jj < end_jj; jj++) {
        double x, y, z;
        vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
        NRLib::Point p1(x, y,    0.0);
        NRLib::Point p2(x, y, 1000.0);
        if (inside_e_cells.IsInsidePolygonXY(p1)) {
          NRLib::Line  line(p1, p2, false, false);
          NRLib::Point intersec_pt;
          if (triangles_elastic[0].FindNearestPoint(line, intersec_pt) < 0.00000000001) {
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            triangles_elastic[2].FindIntersection(line, intersec_pt, true);
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            triangles_elastic[4].FindIntersection(line, intersec_pt, true);
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            for (size_t iii = 0; iii < n_extra_params; ++iii) {
              triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
              (*extra_parameter_grid[iii])(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            }
          }
          else if (triangles_elastic[1].FindNearestPoint(line, intersec_pt) < 0.00000000001) {
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            triangles_elastic[3].FindIntersection(line, intersec_pt, true);
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            triangles_elastic[5].FindIntersection(line, intersec_pt, true);
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            for (size_t iii = 0; iii < n_extra_params; ++iii) {
              triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
              (*extra_parameter_grid[iii])(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            }
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------
void SeismicRegridding::GetCornerPointDir(std::vector<size_t> & a,
                                          std::vector<size_t> & b,
                                          std::vector<size_t> & c,
                                          bool                  left,
                                          bool                  right,
                                          bool                  bot,
                                          bool                  top)
//---------------------------------------------------------------------
{
  if (top || bot || left) {
    a[0] = 0;
    a[1] = 0;
  }
  else {
    a[0] = 1;
    a[1] = 1;
  }
  if (bot || left || right) {
    b[0] = 0;
    b[1] = 0;
  }
  else {
    b[0] = 1;
    b[1] = 1;
  }
  if (top || bot || right) {
    a[2] = 1;
    a[3] = 1;
  }
  else {
    a[2] = 0;
    a[3] = 0;
  }
  if (top || left || right) {
    b[2] = 1;
    b[3] = 1;
  }
  else {
    b[2] = 0;
    b[3] = 0;
  }
  c[0] = 0;
  c[1] = 1;
  c[2] = 0;
  c[3] = 1;
}

//-----------------------------------------------------------------------------------
void SeismicRegridding::FindCornerCellPoints(const NRLib::EclipseGeometry & geometry,
                                             std::vector<NRLib::Point>    & pt_vp,
                                             size_t                         i,
                                             size_t                         j,
                                             size_t                         k,
                                             size_t                         botk)
//-----------------------------------------------------------------------------------
{
  if (k > botk)
    k = k - 1;
  if (i == 0 && j == 0) {//bot left 1243
    pt_vp[0] =        0.5 * (geometry.FindCornerPoint(i, j, k, 0, 0, 0) + geometry.FindCornerPoint(i, j, k, 0, 0, 1));
    pt_vp[1] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 1, 0, 0) + geometry.FindCornerPoint(i, j, k, 1, 0, 1)) + pt_vp[0]);
    pt_vp[3] =               geometry.FindCellCenterPoint(i, j, k);
    pt_vp[2] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 0, 1, 0) + geometry.FindCornerPoint(i, j, k, 0, 1, 1)) + pt_vp[0]);
  }
  else if (i == 0 && j > 0) {//top left 3124
    pt_vp[2] =        0.5 * (geometry.FindCornerPoint(i, j, k, 0, 1, 0) + geometry.FindCornerPoint(i, j, k, 0, 1, 1));
    pt_vp[0] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 0, 0, 0) + geometry.FindCornerPoint(i, j, k, 0, 0, 1)) + pt_vp[2]);
    pt_vp[1] =               geometry.FindCellCenterPoint(i, j, k);
    pt_vp[3] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 1, 1, 0) + geometry.FindCornerPoint(i, j, k, 1, 1, 1)) + pt_vp[2]);
  }
  else if (i > 0 && j == 0) {//bot right 2134
    pt_vp[1] =        0.5 * (geometry.FindCornerPoint(i, j, k, 1, 0, 0) + geometry.FindCornerPoint(i, j, k, 1, 0, 1));
    pt_vp[0] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 0, 0, 0) + geometry.FindCornerPoint(i, j, k, 0, 0, 1)) + pt_vp[1]);
    pt_vp[2] =               geometry.FindCellCenterPoint(i, j, k);
    pt_vp[3] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 1, 1, 0) + geometry.FindCornerPoint(i, j, k, 1, 1, 1)) + pt_vp[1]);
  }
  else {//top right 4213
    pt_vp[3] =        0.5 * (geometry.FindCornerPoint(i, j, k, 1, 1, 0) + geometry.FindCornerPoint(i, j, k, 1, 1, 1));
    pt_vp[1] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 1, 0, 0) + geometry.FindCornerPoint(i, j, k, 1, 0, 1)) + pt_vp[3]);
    pt_vp[0] =               geometry.FindCellCenterPoint(i, j, k);
    pt_vp[2] = 0.5 * (0.5 * (geometry.FindCornerPoint(i, j, k, 0, 1, 0) + geometry.FindCornerPoint(i, j, k, 0, 1, 1)) + pt_vp[3]);
  }
}


//----------------------------------------------------------------------------------------------
void SeismicRegridding::FindCorners(SeismicParameters                   & seismic_parameters,
                                    const NRLib::EclipseGeometry        & geometry,
                                    const NRLib::Grid<double>           & eclipse_vp,
                                    const NRLib::Grid<double>           & eclipse_vs,
                                    const NRLib::Grid<double>           & eclipse_rho,
                                    std::vector<NRLib::Grid<double> >   & eclipse_extra_params,
                                    size_t                                i,
                                    size_t                                j,
                                    size_t                                k,
                                    std::vector<NRLib::Point>           & pt_vp)
//---------------------------------------------------------------------------------------------
{
  NRLib::StormContGrid               & vpgrid                         = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid               & vsgrid                         = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid               & rhogrid                        = seismic_parameters.GetRhoGrid();
  std::vector<NRLib::StormContGrid*>   extra_parameter_grid           = seismic_parameters.GetExtraParametersGrids();

  std::vector<double>                  constvp                        = seismic_parameters.GetModelSettings()->GetConstVp();
  std::vector<double>                  constvs                        = seismic_parameters.GetModelSettings()->GetConstVs();
  std::vector<double>                  constrho                       = seismic_parameters.GetModelSettings()->GetConstRho();
  std::vector<double>                  extra_parameter_default_values = seismic_parameters.GetModelSettings()->GetExtraParameterDefaultValues();

  size_t                               n_extra_params                 = eclipse_extra_params.size();

  size_t                               topk                           = seismic_parameters.GetTopK();
  size_t                               botk                           = seismic_parameters.GetBottomK();

  double                               angle                          = vpgrid.GetAngle();
  double                               cosA                           = cos(angle);
  double                               sinA                           = sin(angle);
  double                               x_min_rot                      = vpgrid.GetXMin()*cosA + vpgrid.GetYMin()*sinA;
  double                               y_min_rot                      = vpgrid.GetYMin()*cosA - vpgrid.GetXMin()*sinA;

  std::vector<double>                  x_rot(4);
  std::vector<double>                  y_rot(4);
  std::vector<bool>                    inside(4);
  std::vector<NRLib::Point>            pt_vs(4);
  std::vector<NRLib::Point>            pt_rho(4);

  std::vector<std::vector<NRLib::Point> > pt_extra_param(n_extra_params);
  for (size_t ii = 0; ii < n_extra_params; ++ii)
    pt_extra_param[ii] = pt_vs;

  for (size_t pt = 0; pt < 4;++pt)
    inside[pt] = vpgrid.IsInside(pt_vp[pt].x, pt_vp[pt].y);

  if (inside[0] || inside[1] || inside[2] || inside[3]) {
    pt_vs[0]  = pt_vp[0];
    pt_rho[0] = pt_vp[0];
    for (size_t ii = 0; ii < n_extra_params; ++ii) {
      pt_extra_param[ii][0] = pt_vp[0];
    }
    if (k == botk + 1) {
      pt_vp[3].z  = constvp[2];
      pt_vs[3].z  = constvs[2];
      pt_rho[3].z = constrho[2];
      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        pt_extra_param[ii][3].z = 0.0;
      }
    }
    else {
      pt_vp[3].z  = eclipse_vp(i, j, k);
      pt_vs[3].z  = eclipse_vs(i, j, k);
      pt_rho[3].z = eclipse_rho(i, j, k);

      for (size_t ii = 0; ii < n_extra_params; ++ii) {
        pt_extra_param[ii][3].z = eclipse_extra_params[ii](i, j, k);
      }
    }

    for (size_t pt = 0; pt < 4; ++pt) {
      x_rot[pt] = pt_vp[pt].x*cosA + pt_vp[pt].y*sinA;
      y_rot[pt] = pt_vp[pt].y*cosA - pt_vp[pt].x*sinA;
    }

    double cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
    double cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
    double cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
    double cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

    size_t start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
    size_t start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
    size_t end_ii   = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
    size_t end_jj   = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));

    if (end_ii > vpgrid.GetNI()) {
      end_ii = vpgrid.GetNI();
    }
    if (end_jj > vpgrid.GetNJ()) {
      end_jj = vpgrid.GetNJ();
    }
    NRLib::Polygon inside_e_cells;
    inside_e_cells.AddPoint(pt_vp[0]);
    inside_e_cells.AddPoint(pt_vp[1]);
    inside_e_cells.AddPoint(pt_vp[3]);
    inside_e_cells.AddPoint(pt_vp[2]);
    for (size_t ii = start_ii; ii < end_ii; ii++) {
      for (size_t jj = start_jj; jj < end_jj; jj++) {
        double x, y, z;
        vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
        NRLib::Point p1(x, y, 0.0);
        if (inside_e_cells.IsInsidePolygonXY(p1)) {
          vpgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(pt_vp[3].z);
          vsgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(pt_vs[3].z);
          rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt_rho[3].z);
          for (size_t iii = 0; iii < n_extra_params; ++iii) {
            NRLib::StormContGrid &param_grid = *(extra_parameter_grid[iii]);
            param_grid(ii, jj, (k - topk) + 1) = static_cast<float>(pt_extra_param[iii][3].z);
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
void SeismicRegridding::PostProcess(SeismicParameters & seismic_parameters,
                                    ModelSettings     *  model_settings)
//--------------------------------------------------------------------------
{
  NRLib::StormContGrid & vpgrid              = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid & vsgrid              = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid & rhogrid             = seismic_parameters.GetRhoGrid();
  float                  missing             = seismic_parameters.GetMissingVal();

  std::vector<double>    constvp             = model_settings->GetConstVp();
  std::vector<double>    constvs             = model_settings->GetConstVs();
  std::vector<double>    constrho            = model_settings->GetConstRho();
  bool                   default_underburden = model_settings->GetDefaultUnderburden();

  size_t                 ni                  = vpgrid.GetNI();
  size_t                 nj                  = vpgrid.GetNJ();
  size_t                 nk                  = vpgrid.GetNK();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nFilling remaining cells in regular grids.\n");

  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  for (size_t i = 0 ; i < ni ; ++i) {
    for (size_t j = 0 ; j < nj ; ++j) {

      bool found_bot = false;
      for (size_t k = nk - 1 ; k > 0 ; --k) {
        if (found_bot && vpgrid(i, j, k) == missing) {
          vpgrid (i, j, k) = static_cast<float>(constvp [1]);
          vsgrid (i, j, k) = static_cast<float>(constvs [1]);
          rhogrid(i, j, k) = static_cast<float>(constrho[1]);
          count1++;
        }
        else if (!found_bot && vpgrid(i, j, k) != missing) {
          found_bot = true;
          if (default_underburden) {
            for (size_t kk = nk - 1; kk > k; --kk) {
              vpgrid (i, j, kk) = static_cast<float>(constvp [2]);
              vsgrid (i, j, kk) = static_cast<float>(constvs [2]);
              rhogrid(i, j, kk) = static_cast<float>(constrho[2]);
              count2++;
            }
          }
          else {
            for (size_t kk = nk - 1; kk > k; --kk) {
              vpgrid (i, j, kk) = vpgrid (i, j, k);
              vsgrid (i, j, kk) = vsgrid (i, j, k);
              rhogrid(i, j, kk) = rhogrid(i, j, k);
              count2++;
            }
          }
        }
        if (!found_bot) {
          for (size_t k = 0 ; k < nk ; ++k) {
            vpgrid (i, j, k) = static_cast<float>(constvp [1]);
            vsgrid (i, j, k) = static_cast<float>(constvs [1]);
            rhogrid(i, j, k) = static_cast<float>(constrho[1]);
            count3++;
          }
        }
      }
    }
  }
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSetting undefined cells in trace equal to default reservoir value    : %10d", count3);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSetting cells in bottom layer equal to default reservoir value       : %10d", count1);
  if (default_underburden)
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSetting cells in trace equal to default underburden                  : %10d\n", count2);
  else
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSetting cells in trace equal to layer below                          : %10d\n", count2);
}

//---------------------------------------------------------------------------------
void SeismicRegridding::FindTWT(SeismicParameters             & seismic_parameters,
                                ModelSettings                 * model_settings,
                                NRLib::RegularSurface<double> & toptime,
                                NRLib::RegularSurface<double> & bottime,
                                size_t                          n_threads)
//---------------------------------------------------------------------------------
{
  NRLib::StormContGrid & vpgrid      = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid & vsgrid      = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid & twtgrid     = seismic_parameters.GetTwtGrid();
  NRLib::StormContGrid & twtssgrid   = seismic_parameters.GetTwtSSGrid();
  NRLib::StormContGrid & twtppgrid   = seismic_parameters.GetTwtPPGrid();
  NRLib::StormContGrid & zgrid       = seismic_parameters.GetZGrid();

  bool                   ps_seismic  = model_settings->GetPSSeismic();
  bool                   nmo_seismic = model_settings->GetNMOCorr();
  double                 v_w         = model_settings->GetVw();
  double                 z_w         = model_settings->GetZw();

  size_t                 ni          = vpgrid.GetNI();
  size_t                 nj          = vpgrid.GetNJ();
  size_t                 nk          = twtgrid.GetNK();
  double                 dx1         = vpgrid.GetDX();
  double                 dy1         = vpgrid.GetDY();
  double                 dx2         = bottime.GetDX();
  double                 dy2         = bottime.GetDY();

  if (ps_seismic) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nAssuming PS seismic.\n");
  }
  if (nmo_seismic) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nAssuming NMO.\n");
  }

#ifdef WITH_OMP
  int  chunk_size = 1;
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif
  for (size_t i = 0; i < ni ; i++) {
    for (size_t j = 0; j < nj ; j++) {
      double x, y, z;
      vpgrid.FindCenterOfCell(i, j, 0, x, y, z);
      twtgrid(i, j, 0) = static_cast<float>(toptime.GetZ(x, y));
      if (ps_seismic && nmo_seismic) {
        float a = 2.0;
        twtppgrid(i, j, 0) = 2.0f/(a + 1.0f)*(twtgrid(i, j, 0) + 1000.0f*(a - 1.0f)*static_cast<float>(z_w/v_w));
        twtssgrid(i, j, 0) = 2.0f*twtgrid(i, j, 0) - twtppgrid(i, j, 0);
      }
      if (!toptime.IsMissing(twtgrid(i, j, 0))) {
        for (size_t k = 1 ; k < nk ; k++) {
          float dz   = zgrid(i, j, k) - zgrid(i, j, k - 1);
          float tfac = 1000.0f;
          if(ps_seismic) {
            twtgrid(i, j, k) = twtgrid(i, j, k - 1) + dz*tfac/vpgrid(i, j, k + 1) + dz*tfac/vsgrid(i, j, k + 1);
          }
          else {
            tfac *= 2.0f;
            twtgrid(i, j, k) = twtgrid(i, j, k - 1) + dz*tfac/vpgrid(i, j, k + 1);
          }
          if (ps_seismic && nmo_seismic){
            tfac *= 2.0f;
            twtppgrid(i, j, k) = twtppgrid(i, j, k - 1) + dz*tfac/vpgrid(i, j, k + 1);
            twtssgrid(i, j, k) = twtssgrid(i, j, k - 1) + dz*tfac/vsgrid(i, j, k + 1);
          }
        }

        double xstart = x - dx1;
        double ystart = y - dy1;
        double xend   = x + dx1;
        double yend   = y + dy1;
        x = xstart;
        y = ystart;
        while (x < xend) {
          y = ystart;
          while (y < yend) {
            size_t ii, jj;
            bottime.FindIndex(x, y, ii, jj);
            bottime(ii, jj) = twtgrid(i, j, nk - 1);
            y = y + dy2;
          }
          x = x + dx2;
        }
      }
      else {
        for (size_t k = 0; k < nk; k++) {
          twtgrid(i, j, k) = -999.0f;
        }
        if (ps_seismic && nmo_seismic) {
          for (size_t k = 0; k < nk; k++) {
            twtppgrid(i, j, k) = -999.0f;
            twtssgrid(i, j, k) = -999.0f;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------
void SeismicRegridding::FindVrms(SeismicParameters          & seismic_parameters,
                                 ModelSettings              * model_settings,
                                 const NRLib::StormContGrid & vgrid,
                                 const NRLib::StormContGrid & twtgrid)
//-------------------------------------------------------------------------------
{
  NRLib::StormContGrid & zgrid    = seismic_parameters.GetZGrid();
  NRLib::StormContGrid & vrmsgrid = seismic_parameters.GetVrmsGrid();

  double                 v_w      = seismic_parameters.GetModelSettings()->GetVw();
  double                 z_w      = seismic_parameters.GetModelSettings()->GetZw();
  double                 twt_w    = 2000*z_w/v_w;

  double                 v_over;
  double                 tmp;
  double                 tmp0;

  for (size_t i = 0; i < vrmsgrid.GetNI(); ++i) {
    for (size_t j = 0; j < vrmsgrid.GetNJ(); ++j) {
      if (twtgrid(i,j,0) == -999.0f) {
        for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
          vrmsgrid(i,j,k) = -999.0f;
        }
      }
      else {
        v_over = 2000*(zgrid(i,j,0) - z_w)/(twtgrid(i,j,0) - 2000*z_w/v_w);
        tmp0 = v_w*v_w*twt_w + v_over*v_over*(twtgrid(i,j,0) - twt_w);
        for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
          tmp = tmp0;
          for (size_t l = 1; l <= k; ++l) {
            tmp += vgrid(i,j,l)* vgrid(i,j,l)*(twtgrid(i,j,l) - twtgrid(i,j,l-1));
          }
          tmp = tmp / twtgrid(i,j,k);
          vrmsgrid(i,j,k) = static_cast<float>(std::sqrt(tmp));
        }
      }
    }
  }
}


//============================ Not looked through =====================================

void SeismicRegridding::WriteElasticParametersSegy(SeismicParameters & seismic_parameters,
                                                   size_t              n_threads,
                                                   bool                time)
{
  std::vector<NRLib::StormContGrid*> input_grid(0);
  std::vector<std::string>           filenames(0);

  NRLib::StormContGrid &vpgrid  = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid &vsgrid  = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid &rhogrid = seismic_parameters.GetRhoGrid();

  input_grid.push_back(&vpgrid);
  input_grid.push_back(&vsgrid);
  input_grid.push_back(&rhogrid);
  if (time) {
    filenames.push_back("vp_time");
    filenames.push_back("vs_time");
    filenames.push_back("rho_time");

    WriteParametersTimeSegy(seismic_parameters, n_threads, input_grid, filenames);
  }
  else {
    filenames.push_back("vp_depth");
    filenames.push_back("vs_depth");
    filenames.push_back("rho_depth");

    WriteParametersDepthSegy(seismic_parameters, n_threads, input_grid, filenames);
  }
  std::cout << "\n";
}

void SeismicRegridding::WriteExtraParametersSegy(SeismicParameters &seismic_parameters,
                                                 size_t             n_threads,
                                                 bool               time)
{
  std::vector<NRLib::StormContGrid*> extra_parameter_grid = seismic_parameters.GetExtraParametersGrids();
  std::vector<std::string> filenames                      = seismic_parameters.GetModelSettings()->GetExtraParameterNames();

  if (time) {
    for (size_t i = 0; i < filenames.size(); ++i) {
      filenames[i] = filenames[i] + "_time";
    }
    WriteParametersTimeSegy(seismic_parameters, n_threads, extra_parameter_grid, filenames);
  }
  else {
    for (size_t i = 0; i < filenames.size(); ++i) {
      filenames[i] = filenames[i] + "_depth";
    }
    WriteParametersDepthSegy(seismic_parameters, n_threads, extra_parameter_grid, filenames);
  }
  std::cout << "\n";
}

void SeismicRegridding::WriteParametersDepthSegy(SeismicParameters                  &seismic_parameters,
                                                 size_t                              n_threads,
                                                 std::vector<NRLib::StormContGrid*>  input_grid,
                                                 std::vector<std::string>            filenames)
{
  bool time = false;
  size_t nz = seismic_parameters.GetSeismicGeometry()->nz();
  double dz = seismic_parameters.GetSeismicGeometry()->dz();
  double z0 = seismic_parameters.GetSeismicGeometry()->z0();

  NRLib::StormContGrid &zgrid = seismic_parameters.GetZGrid();

  //find min z in a sample
  size_t nzmin = static_cast<size_t>(floor(z0) / dz + 0.5);
  double zmin_sampl = nzmin *dz;
  std::vector<double> z_0(nz);
  for (size_t k = 0; k < nz; ++k) {
    z_0[k] = zmin_sampl + (k)* dz;
  }

  WriteParametersSegyInParallel(seismic_parameters, n_threads, input_grid, filenames, z_0, zgrid, time);
}


void SeismicRegridding::WriteParametersTimeSegy(SeismicParameters                  &seismic_parameters,
                                                size_t                              n_threads,
                                                std::vector<NRLib::StormContGrid*>  input_grid,
                                                std::vector<std::string>            filenames)
{
  bool time = true;
  size_t nt = seismic_parameters.GetSeismicGeometry()->nt();
  double dt = seismic_parameters.GetSeismicGeometry()->dt();

  NRLib::RegularSurface<double> &toptime = seismic_parameters.GetTopTime();
  NRLib::StormContGrid          &twtgrid = seismic_parameters.GetTwtGrid();
  double t_min = toptime.Min();

  std::vector<double> twt_0(nt);
  for (size_t k = 0; k < nt; ++k) {
    twt_0[k] = t_min + (k) * dt;
  }

  WriteParametersSegyInParallel(seismic_parameters, n_threads, input_grid, filenames, twt_0, twtgrid, time);
}

void SeismicRegridding::WriteParametersSegyInParallel(SeismicParameters                  &seismic_parameters,
                                                      size_t                              n_threads,
                                                      std::vector<NRLib::StormContGrid*>  input_grid,
                                                      std::vector<std::string>            filenames,
                                                      std::vector<double>                &time_or_depth_vec_reg,
                                                      NRLib::StormContGrid               &time_or_depth_grid,
                                                      bool                                time)
{
  NRLib::RegularSurface<double> &toptime = seismic_parameters.GetTopTime();

  ResamplOutput resampl_output(seismic_parameters, time, time_or_depth_vec_reg.size());
  for (size_t i = 0; i < filenames.size(); ++i) {
    resampl_output.AddResampleCase(filenames[i], *(input_grid[i]), time, time_or_depth_vec_reg, seismic_parameters);
  }
  size_t n_traces;
  tbb::concurrent_queue<Trace*> traces = seismic_parameters.FindTracesInForward(n_traces);
  size_t queue_capacity = seismic_parameters.GetModelSettings()->GetTracesInMemory();

  tbb::concurrent_queue<ResamplTrace*> empty_queue;
  tbb::concurrent_bounded_queue<ResamplTrace*> result_queue;
  result_queue.set_capacity(queue_capacity);
  std::vector<std::thread*> worker_thread;

  GenResamplParam parameters(seismic_parameters, time_or_depth_vec_reg, time_or_depth_grid, toptime, time_or_depth_vec_reg.size(), n_traces, time, empty_queue, result_queue, traces);

  if (n_threads > 1) {
    for (size_t i = 0; i < n_threads - 1; ++i) {
      worker_thread.push_back(new std::thread(GenerateParameterGridForOutputQueue, &parameters, &resampl_output));
    }
    std::thread write_thread(WriteResampledParameter, &parameters, &resampl_output);

    for (size_t i = 0; i < n_threads - 1; ++i) {
      worker_thread[i]->join();
      delete worker_thread[i];
    }
    write_thread.join();
  }
  else {
    float monitor_size, next_monitor;
    seismic_parameters.MonitorInitialize(n_traces, monitor_size, next_monitor);
    for (size_t i = 0; i < n_traces; ++i) {
      Trace *trace;
      parameters.traces.try_pop(trace);
      GenerateParameterGridForOutput(&parameters, trace, &resampl_output);
      ResamplTrace *resampl_trace;
      parameters.result_queue.try_pop(resampl_trace);
      resampl_output.AddTrace(seismic_parameters, parameters.time_or_depth_vec_reg, resampl_trace->GetTraces(), resampl_trace->GetX(), resampl_trace->GetY());
      seismic_parameters.Monitor(i, monitor_size, next_monitor);
      delete trace;
      delete resampl_trace;
    }
  }
  ResamplTrace *resampl_trace;
  while (empty_queue.try_pop(resampl_trace)) {
    delete resampl_trace;
  }
}

void SeismicRegridding::GenerateParameterGridForOutputQueue(GenResamplParam *params,
                                                            ResamplOutput   *resampl_output)
{
  Trace *trace;
  while (params->traces.try_pop(trace)) {
    GenerateParameterGridForOutput(params, trace, resampl_output);
  }
}

void SeismicRegridding::GenerateParameterGridForOutput(GenResamplParam *params,
                                                       Trace           *trace,
                                                       ResamplOutput   *resampl_output)
{
  ResamplTrace *resampl_trace;
  if (!params->empty_queue.try_pop(resampl_trace)) {
    resampl_trace = new ResamplTrace(resampl_output->GetTraces());
  }
  resampl_trace->SetJobID(trace);
  size_t i = trace->GetI();
  size_t j = trace->GetJ();

  double x, y, z;
  std::vector<NRLib::StormContGrid*> input_grid = resampl_output->GetInputGrid();
  input_grid[0]->FindCenterOfCell(i, j, 0, x, y, z);
  double topt = params->toptime.GetZ(x, y);
  bool toptime_missing = params->toptime.IsMissing(topt);
  bool interpolate = params->seismic_parameters.GetModelSettings()->GetResamplParamToSegyInterpol();

  std::vector<double> linear_interp, input_vec(input_grid[0]->GetNK() - 1), input_t(params->time_or_depth_grid->GetNK());

  NRLib::StormContGrid &time_or_depth_grid_ref = *(params->time_or_depth_grid);

  std::vector<NRLib::Grid2D<double> > &output_vec = resampl_trace->GetTraces();
  if (!toptime_missing) { //check whether there are values in input_grid in this pillar - if not, cells in output_grid will be zero
    if (interpolate) {
      for (size_t k = 0; k < params->time_or_depth_grid->GetNK(); ++k) {
        input_t[k] = time_or_depth_grid_ref(i, j, k);
        //input_t[k] = params->time_or_depth_grid(i, j, k);
      }
      for (size_t l = 0; l < output_vec.size(); ++l) {
        NRLib::StormContGrid &input_grid_ref = *(input_grid[l]);
        for (size_t k = 0; k < params->time_or_depth_grid->GetNK(); ++k) {
          input_vec[k] = input_grid_ref(i, j, k);
        }
        linear_interp = params->seismic_parameters.LinInterp1D(input_t, input_vec, params->time_or_depth_vec_reg);
        for (size_t k = 0; k < linear_interp.size(); ++k) {
          if (params->time_or_depth_vec_reg[k] < input_t[0])
            output_vec[l](k, 0) = input_grid_ref(i, j, 0);
          else if (params->time_or_depth_vec_reg[k] > input_t[input_t.size() - 1]) {
            output_vec[l](k, 0) = input_grid_ref(i, j, input_grid_ref.GetNK() - 1);
          }
          else
            output_vec[l](k, 0) = linear_interp[k];
        }
      }
    }
    else {
      for (size_t k = 0; k < output_vec[0].GetNI(); k++) {
        //find cell index in time or depth grid
        double location = params->time_or_depth_vec_reg[k];
        size_t location_index = FindCellIndex(i, j, location, *(params->time_or_depth_grid));
        if (location_index == 999999) {                 // if location is above all values in pillar of time_or_depth_grid,
          location_index = input_grid[0]->GetNK() - 1;  // output_grid is given the value of the bottom cell of input_gridndex) << " " << input_grid[0](i, j, location_index - 1) << "\n";
        }
        for (size_t l = 0; l < output_vec.size(); ++l) {
          output_vec[l](k, 0) = (*input_grid[l])(i, j, location_index);
        }
      }
    }
  }
  else {
    for (size_t k = 0; k < output_vec[0].GetNI(); k++) {
      for (size_t l = 0; l < output_vec.size(); ++l) {
        output_vec[l](k, 0) = 0.0;
      }
    }
  }
  params->result_queue.push(resampl_trace);
}

size_t SeismicRegridding::FindCellIndex(size_t                 i,
                                        size_t                 j,
                                        double                 target_k,
                                        NRLib::StormContGrid & grid)
{
  size_t found_k = 999999;
  size_t nz = grid.GetNK();
  for (size_t k = 0; k < nz; k++) {
    if (grid(i, j, k) > target_k) {
      found_k = k;
      break;
    }
  }
  return found_k;
}

void SeismicRegridding::WriteResampledParameter(GenResamplParam * params,
                                                ResamplOutput   * resampl_output)
{
  float monitor_size, next_monitor;
  params->seismic_parameters.MonitorInitialize(params->n_traces, monitor_size, next_monitor);
  size_t trace = 0;
  std::map<size_t, ResamplTrace*> finished_jobs;
  std::map<size_t, ResamplTrace*>::iterator it;
  while (trace < params->n_traces) {
    ResamplTrace *resampl_trace;
    if (params->result_queue.try_pop(resampl_trace)) {
      finished_jobs.insert(std::pair<size_t, ResamplTrace*>(resampl_trace->GetJobNumber(), resampl_trace));
    }
    it = finished_jobs.find(trace);
    while (it != finished_jobs.end()) {
      resampl_output->AddTrace(params->seismic_parameters, params->time_or_depth_vec_reg, finished_jobs[trace]->GetTraces(), finished_jobs[trace]->GetX(), finished_jobs[trace]->GetY());
      params->empty_queue.push(finished_jobs[trace]);
      finished_jobs.erase(it);
      params->seismic_parameters.Monitor(trace, monitor_size, next_monitor);
      ++trace;
      it = finished_jobs.find(trace);
    }
  }
}
