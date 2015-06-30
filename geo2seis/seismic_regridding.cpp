#include <seismic_regridding.hpp>

#include <nrlib/geometry/triangle.hpp>
#include <nrlib/eclipsegrid/eclipsegrid.hpp>

#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>
#include <physics/wavelet.hpp>
#include <nrlib/random/random.hpp>
#include <nrlib/random/normal.hpp>

#include <seismic_geometry.hpp>

void SeismicRegridding::seismicRegridding(SeismicParameters &seismic_parameters) {
  printf("Start finding Zvalues.\n");
  findZValues(seismic_parameters);
  printf("Zvalues found.\n");

  printf("Start finding elastic parameters.\n");
  findVpAndR(seismic_parameters);
  printf("Elastic parameters found.\n");


  seismic_parameters.deleteEclipseGrid();

  if (seismic_parameters.modelSettings()->GetNMOCorr() == false) {
    // if we are adding noise to reflection write noise-free reflections to file before adding noise
    if (seismic_parameters.modelSettings()->GetOutputReflections()) {
      seismic_parameters.seismicOutput()->writeReflections(seismic_parameters, false);
    }

    //-----------------Add white noise to reflections----------------------------
    if (seismic_parameters.modelSettings()->GetWhiteNoise()) {
      unsigned long seed = seismic_parameters.modelSettings()->GetSeed();
      double deviation = seismic_parameters.modelSettings()->GetStandardDeviation();
      addNoiseToReflections(seed, deviation, seismic_parameters.rGrids());
    }

    if (seismic_parameters.modelSettings()->GetOutputReflections() && seismic_parameters.modelSettings()->GetWhiteNoise()) {
      seismic_parameters.seismicOutput()->writeReflections(seismic_parameters, true);
    }
  }

  NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
  NRLib::RegularSurface<double> &bottime = seismic_parameters.bottomTime();
  NRLib::StormContGrid &vpgrid           = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid           = seismic_parameters.vsGrid();
  NRLib::StormContGrid &twtgrid          = seismic_parameters.twtGrid();
  NRLib::StormContGrid &zgrid            = seismic_parameters.zGrid();
  bool find_for_ps = seismic_parameters.modelSettings()->GetPSSeismic();
  findTWT(vpgrid, vsgrid, twtgrid, zgrid, toptime, bottime, find_for_ps);

  if (seismic_parameters.modelSettings()->GetNMOCorr()){  
    if (seismic_parameters.modelSettings()->GetOutputVrms()){
      findVrms(seismic_parameters);
      printf("Write rms velocity.\n");
      seismic_parameters.seismicOutput()->writeVrms(seismic_parameters);
      seismic_parameters.deleteVrmsGrid();
    }
    seismic_parameters.deleteVrmsGrid();
  }

  std::vector<double> constvp = seismic_parameters.modelSettings()->GetConstVp();
  Wavelet* wavelet = seismic_parameters.wavelet();

  toptime.Add(-2000 / constvp[0] * wavelet->GetDepthAdjustmentFactor()); // add one wavelet length to bot and subtract from top
  bottime.Add(2000 / constvp[2] * wavelet->GetDepthAdjustmentFactor());

  double tmin = toptime.Min();
  double tmax = bottime.Max();
  size_t nt = static_cast<size_t>(floor((tmax - tmin) / seismic_parameters.seismicGeometry()->dt()+0.5));
  if (seismic_parameters.modelSettings()->GetNLayersFileName() != "") {
    NRLib::StormContGrid tmpgrid(seismic_parameters.modelSettings()->GetNLayersFileName());
    nt = tmpgrid.GetNK();
  }
  seismic_parameters.seismicGeometry()->setNt(nt);
  seismic_parameters.seismicGeometry()->setTRange(tmin, tmax);
  

  //---------------Print toptime and bottime------------------------
  if (seismic_parameters.modelSettings()->GetOutputTimeSurfaces()) {
    seismic_parameters.seismicOutput()->writeTimeSurfaces(seismic_parameters);
  }

  //----------generating grid for output in segy format-----------------

  //time grid:
  if (seismic_parameters.modelSettings()->GetOutputElasticParametersTimeSegy()) {
    seismic_parameters.seismicOutput()->writeElasticParametersTimeSegy(seismic_parameters);
  }

  if (seismic_parameters.modelSettings()->GetOutputExtraParametersTimeSegy()) {
    seismic_parameters.seismicOutput()->writeExtraParametersTimeSegy(seismic_parameters);
  }

  //depth grid:
  if (seismic_parameters.modelSettings()->GetOutputElasticParametersDepthSegy()) {
    seismic_parameters.seismicOutput()->writeElasticParametersDepthSegy(seismic_parameters);
  }

  if (seismic_parameters.modelSettings()->GetOutputExtraParametersDepthSegy()) {
    seismic_parameters.seismicOutput()->writeExtraParametersDepthSegy(seismic_parameters);
  }

  if (seismic_parameters.modelSettings()->GetOutputVp()) {
    seismic_parameters.seismicOutput()->writeVpVsRho(seismic_parameters);
  }

  if (seismic_parameters.modelSettings()->GetOutputZvalues()) {
    seismic_parameters.seismicOutput()->writeZValues(seismic_parameters);
  }

  if (seismic_parameters.modelSettings()->GetOutputTwt()) {
    seismic_parameters.seismicOutput()->writeTwt(seismic_parameters);
  }


  //printf("remove grids\n");
  seismic_parameters.deleteExtraParameterGrids();

}


void SeismicRegridding::findZValues(SeismicParameters &seismic_parameters) {
  NRLib::StormContGrid &zgrid = seismic_parameters.zGrid();
  const NRLib::EclipseGeometry &geometry = seismic_parameters.eclipseGrid().GetGeometry();
  size_t top_k = seismic_parameters.topK();
  bool use_corner_point = seismic_parameters.modelSettings()->GetUseCornerpointInterpol();

  double xmin  = zgrid.GetXMin();
  double ymin  = zgrid.GetYMin();
  double dx    = zgrid.GetDX();
  double dy    = zgrid.GetDY();
  double angle = zgrid.GetAngle();

  for (size_t k = 0; k < zgrid.GetNK() - 1; k++) {
    NRLib::Grid2D<double> values(zgrid.GetNI(), zgrid.GetNJ(), 0);
    if (use_corner_point) {
      geometry.FindLayerSurfaceCornerpoint(values, k + top_k, 0, dx, dy, xmin, ymin, angle, 0);
    } else {
      geometry.FindLayerSurface(values, k + top_k, 0, dx, dy, xmin, ymin, angle, 0);
    }
    for (size_t i = 0; i < zgrid.GetNI(); i++) {
      for (size_t j = 0; j < zgrid.GetNJ(); j++) {
        zgrid(i, j, k) = static_cast<float>(values(i, j));
      }
    }
  }

  size_t k = zgrid.GetNK() - 2;
  NRLib::Grid2D<double> values(zgrid.GetNI(), zgrid.GetNJ(), 0);
  if (use_corner_point) {
    geometry.FindLayerSurfaceCornerpoint(values, k + top_k, 1, dx, dy, xmin, ymin, angle, 0);
  } else {
    geometry.FindLayerSurface(values, k + top_k, 1, dx, dy, xmin, ymin, angle, 0);
  }

  for (size_t i = 0; i < zgrid.GetNI(); i++) {
    for (size_t j = 0; j < zgrid.GetNJ(); j++) {
      zgrid(i, j, k + 1) = static_cast<float>(values(i, j));
    }
  }

}

void SeismicRegridding::findVrms(SeismicParameters &seismic_parameters){
  double v_w = seismic_parameters.modelSettings()->GetVw();
  double z_w = seismic_parameters.modelSettings()->GetZw();
  NRLib::StormContGrid &zgrid             = seismic_parameters.zGrid();
  NRLib::StormContGrid &twtgrid           = seismic_parameters.twtGrid();
  NRLib::StormContGrid &vrmsgrid          = seismic_parameters.vrmsGrid();
  NRLib::StormContGrid &vpgrid            = seismic_parameters.vpGrid();

  double v_over;
  double twt_w = 2000*z_w/v_w;
  double tmp, tmp0;
  for (size_t i = 0; i < vrmsgrid.GetNI(); ++i) {
    for (size_t j = 0; j < vrmsgrid.GetNJ(); ++j) {
      if (twtgrid(i,j,0) == -999.0) {
        for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
          vrmsgrid(i,j,k) = -999.0;
        }
      }
      else {
        v_over = 2000*(zgrid(i,j,0) - z_w)/(twtgrid(i,j,0) - 2000*z_w/v_w);
        tmp0 = v_w*v_w*twt_w + v_over*v_over*(twtgrid(i,j,0) - twt_w);
        for (size_t k = 0; k < vrmsgrid.GetNK(); ++k) {
          tmp = tmp0;
          for (size_t l = 1; l <= k; ++l) {
            tmp += vpgrid(i,j,l)* vpgrid(i,j,l)*(twtgrid(i,j,l) - twtgrid(i,j,l-1));
          }
          tmp = tmp / twtgrid(i,j,k);
          vrmsgrid(i,j,k) = float(std::sqrt(tmp));
        }
      }
    }
  }
}


void SeismicRegridding::findVpAndR(SeismicParameters &seismic_parameters) {
  NRLib::StormContGrid &vpgrid  = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid  = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();

  std::vector<NRLib::StormContGrid> &rgrid                = seismic_parameters.rGrids();
  std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();
  const NRLib::EclipseGrid          &egrid                = seismic_parameters.eclipseGrid();

  size_t topk    = seismic_parameters.topK();
  size_t botk    = seismic_parameters.bottomK();
  double theta_0 = seismic_parameters.theta0();
  double dtheta  = seismic_parameters.dTheta();
  std::vector<double> constvp   = seismic_parameters.modelSettings()->GetConstVp();
  std::vector<double> constvs   = seismic_parameters.modelSettings()->GetConstVs();
  std::vector<double> constrho  = seismic_parameters.modelSettings()->GetConstRho();

  std::vector<std::string> names                     = seismic_parameters.modelSettings()->GetParameterNames();
  std::vector<double> extra_parameter_default_values = seismic_parameters.modelSettings()->GetExtraParameterDefaultValues();
  std::vector<std::string> extra_parameter_names     = seismic_parameters.modelSettings()->GetExtraParameterNames();
  NRLib::StormContGrid &zgrid                        = seismic_parameters.zGrid();
  //NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();

  double zlimit = seismic_parameters.modelSettings()->GetZeroThicknessLimit();
  bool ps_seis  = seismic_parameters.modelSettings()->GetPSSeismic();

  const NRLib::EclipseGeometry &geometry = egrid.GetGeometry();
  NRLib::Point pt1vp, pt2vp, pt3vp, pt4vp;
  NRLib::Point pt1vs, pt2vs, pt3vs, pt4vs;
  NRLib::Point pt1rho, pt2rho, pt3rho, pt4rho;
  std::vector<NRLib::Point> pt1_extra_param(extra_parameter_names.size()), pt2_extra_param(extra_parameter_names.size());
  std::vector<NRLib::Point> pt3_extra_param(extra_parameter_names.size()), pt4_extra_param(extra_parameter_names.size());

  double diffvp, meanvp, diffvs, meanvs, diffrho, meanrho;
  std::vector<Zoeppritz *> zoeppritz(rgrid.size());
  size_t n_angle = rgrid.size();
  if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
    double theta;
    for (size_t i = 0; i < n_angle; i++) {
      theta = theta_0 + i * dtheta;
      if (ps_seis) {
        zoeppritz[i] = new ZoeppritzPS();
      } else {
        zoeppritz[i] = new ZoeppritzPP();
      }
      zoeppritz[i]->ComputeConstants(theta);
    }
  }
  NRLib::Grid2D<double> value_above_vp(egrid.GetNI(), egrid.GetNJ(), constvp[0]);
  NRLib::Grid2D<double> value_above_vs(egrid.GetNI(), egrid.GetNJ(), constvs[0]);
  NRLib::Grid2D<double> value_above_rho(egrid.GetNI(), egrid.GetNJ(), constrho[0]);
  std::vector<NRLib::Grid2D<double> > value_above_extra_param;
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    NRLib::Grid2D<double> value_above_grid(egrid.GetNI(), egrid.GetNJ(), 0.0);
    value_above_extra_param.push_back(value_above_grid);
  }

  double vp_angle = vpgrid.GetAngle();
  double cosvpangle = cos(vp_angle);
  double sinvpangle = sin(vp_angle);
  double x_min_rot = vpgrid.GetXMin() * cos(vp_angle) + vpgrid.GetYMin() * sin(vp_angle);
  double y_min_rot = vpgrid.GetYMin() * cos(vp_angle) - vpgrid.GetXMin() * sin(vp_angle);
  double cell_min_x, cell_max_x, cell_min_y, cell_max_y, x1_rot, x2_rot, x3_rot, x4_rot, y1_rot, y2_rot, y3_rot, y4_rot;
  size_t start_ii, start_jj, end_ii, end_jj;

  const NRLib::Grid<double> &vp_grid = egrid.GetParameter(names[0]);
  const NRLib::Grid<double> &vs_grid = egrid.GetParameter(names[1]);
  const NRLib::Grid<double> &rho_grid = egrid.GetParameter(names[2]);
  std::vector<NRLib::Grid<double> > parameter_grid_from_eclipse;
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    const NRLib::Grid<double> &one_parameter_grid = egrid.GetParameter(extra_parameter_names[i]);
    parameter_grid_from_eclipse.push_back(one_parameter_grid);
  }

  for (size_t i = 0; i < vpgrid.GetNI(); i++) {
    for (size_t j = 0; j < vpgrid.GetNJ(); j++) {
      vpgrid(i, j, 0) = static_cast<float>(constvp[0]);
      vsgrid(i, j, 0) = static_cast<float>(constvs[0]);
      rhogrid(i, j, 0) = static_cast<float>(constrho[0]);
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        extra_parameter_grid[ii](i, j, 0) = 0.0;
      }
    }
  }
  for (size_t k = topk; k <= botk + 1; k++) {
    //   printf("k = %d \n", k);
    for (size_t i = 0; i < egrid.GetNI() - 1; i++) {
      for (size_t j = 0; j < egrid.GetNJ() - 1; j++) {
        if (geometry.IsPillarActive(i, j) && geometry.IsPillarActive(i + 1, j) && geometry.IsPillarActive(i, j + 1) && geometry.IsPillarActive(i + 1, j + 1) &&
          geometry.IsPillarActive(i + 2, j) && geometry.IsPillarActive(i + 2, j + 1) && geometry.IsPillarActive(i, j + 2) && geometry.IsPillarActive(i + 1, j + 2) &&
          geometry.IsPillarActive(i + 2, j + 2)) {
            if (k <= botk) {
              pt1vp = geometry.FindCellCenterPoint(i, j, k);
              pt2vp = geometry.FindCellCenterPoint(i + 1, j, k);
              pt3vp = geometry.FindCellCenterPoint(i, j + 1, k);
              pt4vp = geometry.FindCellCenterPoint(i + 1, j + 1, k);
            } else {
              pt1vp = geometry.FindCellCenterPoint(i, j, k - 1);
              pt2vp = geometry.FindCellCenterPoint(i + 1, j, k - 1);
              pt3vp = geometry.FindCellCenterPoint(i, j + 1, k - 1);
              pt4vp = geometry.FindCellCenterPoint(i + 1, j + 1, k - 1);
            }
            int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
            int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
            int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
            int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
            if (inside1 || inside2 || inside3 || inside4) {
              pt1vs = pt1vp;
              pt1rho = pt1vp;
              pt2vs = pt2vp;
              pt2rho = pt2vp;
              pt3vs = pt3vp;
              pt3rho = pt3vp;
              pt4vs = pt4vp;
              pt4rho = pt4vp;
              for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                pt1_extra_param[ii] = pt1vp;
                pt2_extra_param[ii] = pt2vp;
                pt3_extra_param[ii] = pt3vp;
                pt4_extra_param[ii] = pt4vp;
              }
              if (k == botk + 1) {
                pt1vp.z = constvp[2];
                pt2vp.z = constvp[2];
                pt3vp.z = constvp[2];
                pt4vp.z = constvp[2];
                pt1vs.z = constvs[2];
                pt2vs.z = constvs[2];
                pt3vs.z = constvs[2];
                pt4vs.z = constvs[2];
                pt1rho.z = constrho[2];
                pt2rho.z = constrho[2];
                pt3rho.z = constrho[2];
                pt4rho.z = constrho[2];
                for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                  pt1_extra_param[ii].z = 0.0;
                  pt2_extra_param[ii].z = 0.0;
                  pt3_extra_param[ii].z = 0.0;
                  pt4_extra_param[ii].z = 0.0;
                }
              } else {
                findPointZValue(i, j, k, pt1vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
                findPointZValue(i, j, k, pt1vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
                findPointZValue(i, j, k, pt1rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

                findPointZValue(i + 1, j, k, pt2vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
                findPointZValue(i + 1, j, k, pt2vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
                findPointZValue(i + 1, j, k, pt2rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

                findPointZValue(i, j + 1, k, pt3vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
                findPointZValue(i, j + 1, k, pt3vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
                findPointZValue(i, j + 1, k, pt3rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

                findPointZValue(i + 1, j + 1, k, pt4vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
                findPointZValue(i + 1, j + 1, k, pt4vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
                findPointZValue(i + 1, j + 1, k, pt4rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

                for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                  findPointZValue(i, j, k, pt1_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
                  findPointZValue(i + 1, j, k, pt2_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
                  findPointZValue(i, j + 1, k, pt3_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
                  findPointZValue(i + 1, j + 1, k, pt4_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
                }
              }
              value_above_vp(i, j) = pt1vp.z;
              value_above_vs(i, j) = pt1vs.z;
              value_above_rho(i, j) = pt1rho.z;
              value_above_vp(i + 1, j) = pt2vp.z;
              value_above_vs(i + 1, j) = pt2vs.z;
              value_above_rho(i + 1, j) = pt2rho.z;
              value_above_vp(i, j + 1) = pt3vp.z;
              value_above_vs(i, j + 1) = pt3vs.z;
              value_above_rho(i, j + 1) = pt3rho.z;
              value_above_vp(i + 1, j + 1) = pt4vp.z;
              value_above_vs(i + 1, j + 1) = pt4vs.z;
              value_above_rho(i + 1, j + 1) = pt4rho.z;
              for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                value_above_extra_param[ii](i, j) = pt1_extra_param[ii].z;
                value_above_extra_param[ii](i + 1, j) = pt2_extra_param[ii].z;
                value_above_extra_param[ii](i, j + 1) = pt3_extra_param[ii].z;
                value_above_extra_param[ii](i + 1, j + 1) = pt4_extra_param[ii].z;
              }
              bool triangulate_124 = true;
              NRLib::Point vec1, vec2;
              vec1 = pt1vp - pt2vp;
              vec1.z = 0;
              vec2 = pt4vp - pt2vp;
              vec2.z = 0;
              double delaunay_angle = vec1.GetAngle(vec2);
              vec1 = pt1vp - pt3vp;
              vec1.z = 0;
              vec2 = pt4vp - pt3vp;
              vec2.z = 0;
              delaunay_angle += vec1.GetAngle(vec2);
              if (delaunay_angle > NRLib::Pi) {
                triangulate_124 = false;
              }

              NRLib::Triangle triangle3(pt1vp, pt2vp, pt3vp);
              NRLib::Triangle triangle4(pt2vp, pt3vp, pt4vp);
              NRLib::Triangle triangle7(pt1vs, pt2vs, pt3vs);
              NRLib::Triangle triangle8(pt2vs, pt3vs, pt4vs);
              NRLib::Triangle triangle11(pt1rho, pt2rho, pt3rho);
              NRLib::Triangle triangle12(pt2rho, pt3rho, pt4rho);
              std::vector<NRLib::Triangle> triangles_extra_param;
              for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                NRLib::Triangle triangle_temp(pt1_extra_param[ii], pt2_extra_param[ii], pt3_extra_param[ii]);
                triangles_extra_param.push_back(triangle_temp);
                triangle_temp.SetCornerPoints(pt2_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
                triangles_extra_param.push_back(triangle_temp);
              }
              if (triangulate_124) {
                triangle3.SetCornerPoints(pt1vp, pt2vp, pt4vp);
                triangle4.SetCornerPoints(pt1vp, pt3vp, pt4vp);
                triangle7.SetCornerPoints(pt1vs, pt2vs, pt4vs);
                triangle8.SetCornerPoints(pt1vs, pt3vs, pt4vs);
                triangle11.SetCornerPoints(pt1rho, pt2rho, pt4rho);
                triangle12.SetCornerPoints(pt1rho, pt3rho, pt4rho);
                for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                  triangles_extra_param[ii * 2].SetCornerPoints(pt1_extra_param[ii], pt2_extra_param[ii], pt4_extra_param[ii]);
                  triangles_extra_param[ii * 2 + 1].SetCornerPoints(pt1_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
                }
              }

              x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
              y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
              x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
              y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
              x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
              y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
              x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
              y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

              cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
              cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
              cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
              cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

              start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 0.5));
              start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 0.5));
              end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 1.0));
              end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 1.0));
              // start_ii=static_cast<unsigned int>(max(0.0,(cell_min_x-x_min_rot)/vpgrid.GetDX()-1.0));
              //  start_jj=static_cast<unsigned int>(max(0.0,(cell_min_y-y_min_rot)/vpgrid.GetDY()-1.0));
              //  end_ii=static_cast<unsigned int>(max(0.0,(cell_max_x-x_min_rot)/vpgrid.GetDX()+1.5));
              //  end_jj=static_cast<unsigned int>(max(0.0,(cell_max_y-y_min_rot)/vpgrid.GetDY()+1.5));
              if (end_ii > vpgrid.GetNI()) {
                end_ii = vpgrid.GetNI();
              }
              if (end_jj > vpgrid.GetNJ()) {
                end_jj = vpgrid.GetNJ();
              }
              for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
                for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
                  double x, y, z;
                  vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
                  //         if(toptime.IsMissing(toptime.GetZ(x,y))==false){
                  if (true) { // earlier: x >= minx && x <= maxx && y >= miny && y <=maxy
                    NRLib::Point p1, p2;
                    p1.x = x;
                    p1.y = y;
                    p1.z = pt1vp.z;
                    p2 = p1;
                    p2.z += 1000;
                    NRLib::Line line(p1, p2, false, false);
                    NRLib::Point intersec_pt;
                    // bool intersect = triangle3.FindIntersection(line, intersec_pt, true);
                    double dist = triangle3.FindNearestPoint(line, intersec_pt); // To avoid numerical instabilities when point is on edge of triangle
                    //if(intersect==true){
                    if (dist < 0.00000000001) {
                      vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                      meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                      triangle7.FindIntersection(line, intersec_pt, true);              //<--m� ha denne!
                      vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z); //<--m� ha denne!
                      diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                      meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                      triangle11.FindIntersection(line, intersec_pt, true);
                      rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                      meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                      for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                        triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                        extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      }
                    } else {
                      dist = triangle4.FindNearestPoint(line, intersec_pt);
                      // intersect = triangle4.FindIntersection(line, intersec_pt, true);
                      // if(intersect == true){
                      if (dist < 0.00000000001) {
                        vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                        if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                          diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                          meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                        }
                        triangle8.FindIntersection(line, intersec_pt, true);
                        vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                        if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                          diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                          meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                        }
                        triangle12.FindIntersection(line, intersec_pt, true);
                        rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                        if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                          diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                          meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                        }
                        for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                          triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                          extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                        }
                      }
                    }
                    // if(intersect){
                    if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                      if (dist < 0.00000000001) {
                        for (size_t l = 0; l < n_angle; l++) {
                          rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
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


    NRLib::Point mid_edge1, mid_edge2;
    int j_out = 0;
    int i_out;
    //------------------------Filling in j_out=0 bot half edge------------------------------------

    for (size_t i = 0; i < egrid.GetNI() - 1; i++) { //Filling in j_out=0 bot half edge
      j_out = 0;
      //if(geometry.IsPillarActive(i,j_out) && geometry.IsPillarActive(i+1,j_out) && geometry.IsPillarActive(i, j_out+1) && geometry.IsPillarActive(i+1, j_out+1) &&
      //         geometry.IsPillarActive(i+2,j_out) && geometry.IsPillarActive(i+2, j_out+1) ){
      while (j_out < static_cast<int>(egrid.GetNJ()) && !(geometry.IsPillarActive(i, j_out) && geometry.IsPillarActive(i + 1, j_out) && geometry.IsPillarActive(i, j_out + 1) && geometry.IsPillarActive(i + 1, j_out + 1) &&
        geometry.IsPillarActive(i + 2, j_out) && geometry.IsPillarActive(i + 2, j_out + 1))) {
          j_out++;
      }
      if (j_out < static_cast<int>(egrid.GetNJ())) {
        if (k <= botk) {
          pt1vp = geometry.FindCellCenterPoint(i, j_out, k);
          pt2vp = geometry.FindCellCenterPoint(i + 1, j_out, k);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i, j_out, k, 0, 0, 0) + geometry.FindCornerPoint(i, j_out, k, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i, j_out, k, 1, 0, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k, 0, 0, 0) + geometry.FindCornerPoint(i + 1, j_out, k, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i + 1, j_out, k, 1, 0, 1));
        } else {
          pt1vp = geometry.FindCellCenterPoint(i, j_out, k - 1);
          pt2vp = geometry.FindCellCenterPoint(i + 1, j_out, k - 1);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 1, 0, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i + 1, j_out, k - 1, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i + 1, j_out, k - 1, 1, 0, 1));
        }

        //mid_edge2 = geometry.FindCornerPoint(i+1,j_out,k,0,0,0)+geometry.FindCornerPoint(i+1,j_out,k,1,0,0);
        int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
        int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
        pt3vp = mid_edge1 - pt1vp;
        int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
        mid_edge1 = 0.5 * mid_edge1;
        pt4vp = mid_edge2 - pt2vp;
        int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
        mid_edge2 = 0.5 * mid_edge2;
        if (inside1 || inside2 || inside3 || inside4) {
          pt1vs = pt1vp;
          pt1rho = pt1vp;
          pt2vs = pt2vp;
          pt2rho = pt2vp;
          pt3vs = pt3vp;
          pt3rho = pt3vp;
          pt4vs = pt4vp;
          pt4rho = pt4vp;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt1_extra_param[ii] = pt1vp;
            pt2_extra_param[ii] = pt2vp;
            pt3_extra_param[ii] = pt3vp;
            pt4_extra_param[ii] = pt4vp;
          }
          if (k == botk + 1) {
            pt1vp.z = constvp[2];
            pt2vp.z = constvp[2];
            pt1vs.z = constvs[2];
            pt2vs.z = constvs[2];
            pt1rho.z = constrho[2];
            pt2rho.z = constrho[2];
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              pt1_extra_param[ii].z = 0.0;
              pt2_extra_param[ii].z = 0.0;
            }
          } else {
            findPointZValue(i, j_out, k, pt1vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i, j_out, k, pt1vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i, j_out, k, pt1rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            findPointZValue(i + 1, j_out, k, pt2vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i + 1, j_out, k, pt2vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i + 1, j_out, k, pt2rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              findPointZValue(i, j_out, k, pt1_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
              findPointZValue(i + 1, j_out, k, pt2_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
            }
          }
          pt3vp.z = pt1vp.z;
          pt4vp.z = pt2vp.z;
          pt3vs.z = pt1vs.z;
          pt4vs.z = pt2vs.z;
          pt3rho.z = pt1rho.z;
          pt4rho.z = pt2rho.z;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt3_extra_param[ii].z = pt1_extra_param[ii].z;
            pt4_extra_param[ii].z = pt2_extra_param[ii].z;
          }

          bool triangulate_124 = true;
          NRLib::Point vec1, vec2;
          vec1 = pt1vp - pt2vp;
          vec1.z = 0;
          vec2 = pt4vp - pt2vp;
          vec2.z = 0;
          double delaunay_angle = vec1.GetAngle(vec2);
          vec1 = pt1vp - pt3vp;
          vec1.z = 0;
          vec2 = pt4vp - pt3vp;
          vec2.z = 0;
          delaunay_angle += vec1.GetAngle(vec2);
          if (delaunay_angle > NRLib::Pi) {
            triangulate_124 = false;
          }

          NRLib::Triangle triangle3(pt1vp, pt2vp, pt3vp);
          NRLib::Triangle triangle4(pt2vp, pt3vp, pt4vp);
          NRLib::Triangle triangle7(pt1vs, pt2vs, pt3vs);
          NRLib::Triangle triangle8(pt2vs, pt3vs, pt4vs);
          NRLib::Triangle triangle11(pt1rho, pt2rho, pt3rho);
          NRLib::Triangle triangle12(pt2rho, pt3rho, pt4rho);
          std::vector<NRLib::Triangle> triangles_extra_param;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            NRLib::Triangle triangle_temp(pt1_extra_param[ii], pt2_extra_param[ii], pt3_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
            triangle_temp.SetCornerPoints(pt2_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
          }
          if (triangulate_124) {
            triangle3.SetCornerPoints(pt1vp, pt2vp, pt4vp);
            triangle4.SetCornerPoints(pt1vp, pt3vp, pt4vp);
            triangle7.SetCornerPoints(pt1vs, pt2vs, pt4vs);
            triangle8.SetCornerPoints(pt1vs, pt3vs, pt4vs);
            triangle11.SetCornerPoints(pt1rho, pt2rho, pt4rho);
            triangle12.SetCornerPoints(pt1rho, pt3rho, pt4rho);
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              triangles_extra_param[ii * 2].SetCornerPoints(pt1_extra_param[ii], pt2_extra_param[ii], pt4_extra_param[ii]);
              triangles_extra_param[ii * 2 + 1].SetCornerPoints(pt1_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            }
          }
          x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
          y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
          x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
          y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
          x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
          y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
          x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
          y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

          cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
          cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
          cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
          cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

          start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
          start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
          end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
          end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
          if (end_ii > vpgrid.GetNI()) {
            end_ii = vpgrid.GetNI();
          }
          if (end_jj > vpgrid.GetNJ()) {
            end_jj = vpgrid.GetNJ();
          }
          NRLib::Polygon inside_e_cells;
          inside_e_cells.AddPoint(pt1vp);
          inside_e_cells.AddPoint(pt2vp);
          inside_e_cells.AddPoint(mid_edge2);
          // if(k > 0)
          //   inside_e_cells.AddPoint(geometry.FindCornerPoint(i,j_out,k-1,1,0,1));
          // else
          if (k <= botk) {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i, j_out, k, 1, 0, 1)));
          } else {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 1, 0, 1)));
          }
          inside_e_cells.AddPoint(mid_edge1);
          for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
            for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
              double x, y, z;
              vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
              //   if(toptime.IsMissing(toptime.GetZ(x,y))==false){
              NRLib::Point p1(x, y, 0.0);
              NRLib::Point p2(x, y, 1000.0);
              if (inside_e_cells.IsInsidePolygonXY(p1)) {
                NRLib::Line line(p1, p2, false, false);
                NRLib::Point intersec_pt;
                // bool intersect = triangle3.FindIntersection(line, intersec_pt, true);
                double dist = triangle3.FindNearestPoint(line, intersec_pt);
                // if(intersect==true){
                if (dist < 0.00000000001) {
                  vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                  meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                  triangle7.FindIntersection(line, intersec_pt, true);
                  vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                  meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                  triangle11.FindIntersection(line, intersec_pt, true);
                  rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                  meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                  for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                    triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                    extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  }
                } else {
                  dist = triangle4.FindNearestPoint(line, intersec_pt);
                  //intersect = triangle4.FindIntersection(line, intersec_pt, true);
                  // if(intersect == true){
                  if (dist < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                    meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                    triangle8.FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                    meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                    triangle12.FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                    meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                    for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                      triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                      extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                }
                // if(intersect){
                if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                  if (dist < 0.00000000001) {
                    for (size_t l = 0; l < n_angle; l++) {
                      rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
                    }
                  }
                }
              }
            }
          }
        }
      }
    } //end for loop over i along j=0 layer




    j_out = static_cast<int>(egrid.GetNJ()) - 1;
    //------------------------Filling in j_out=egrid.GetNJ()-1 top half edge------------------------------------

    for (size_t i = 0; i < egrid.GetNI() - 1; i++) { //Filling in j_out=egrid.GetNJ()-1 top half edge
      j_out = static_cast<int>(egrid.GetNJ()) - 1;
      //if(geometry.IsPillarActive(i,j_out) && geometry.IsPillarActive(i+1,j_out) && geometry.IsPillarActive(i, j_out+1) && geometry.IsPillarActive(i+1, j_out+1) &&
      //         geometry.IsPillarActive(i+2,j_out) && geometry.IsPillarActive(i+2, j_out+1) ){
      while (j_out >= 0 && !(geometry.IsPillarActive(i, j_out) && geometry.IsPillarActive(i + 1, j_out) && geometry.IsPillarActive(i, j_out + 1) && geometry.IsPillarActive(i + 1, j_out + 1) &&
        geometry.IsPillarActive(i + 2, j_out) && geometry.IsPillarActive(i + 2, j_out + 1))) {
          j_out--;
      }
      if (j_out >= 0) {
        if (k <= botk) {
          pt1vp = geometry.FindCellCenterPoint(i, j_out, k);
          pt2vp = geometry.FindCellCenterPoint(i + 1, j_out, k);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i, j_out, k, 0, 1, 0) + geometry.FindCornerPoint(i, j_out, k, 0, 1, 1)) + 0.5 * (geometry.FindCornerPoint(i, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i, j_out, k, 1, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k, 0, 1, 0) + geometry.FindCornerPoint(i + 1, j_out, k, 0, 1, 1)) + 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i + 1, j_out, k, 1, 1, 1));
        } else {
          pt1vp = geometry.FindCellCenterPoint(i, j_out, k - 1);
          pt2vp = geometry.FindCellCenterPoint(i + 1, j_out, k - 1);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 0, 1, 1)) + 0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 1, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i + 1, j_out, k - 1, 0, 1, 1)) + 0.5 * (geometry.FindCornerPoint(i + 1, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i + 1, j_out, k - 1, 1, 1, 1));
        }
        // mid_edge2 = geometry.FindCornerPoint(i,j_out+1,k,0,0,0)+geometry.FindCornerPoint(i,j_out+1,k,0,1,0);

        int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
        int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
        pt3vp = mid_edge1 - pt1vp;
        int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
        mid_edge1 = 0.5 * mid_edge1;
        pt4vp = mid_edge2 - pt2vp;
        int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
        mid_edge2 = 0.5 * mid_edge2;
        if (inside1 || inside2 || inside3 || inside4) {
          pt1vs = pt1vp;
          pt1rho = pt1vp;
          pt2vs = pt2vp;
          pt2rho = pt2vp;
          pt3vs = pt3vp;
          pt3rho = pt3vp;
          pt4vs = pt4vp;
          pt4rho = pt4vp;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt1_extra_param[ii] = pt1vp;
            pt2_extra_param[ii] = pt2vp;
            pt3_extra_param[ii] = pt3vp;
            pt4_extra_param[ii] = pt4vp;
          }
          if (k == botk + 1) {
            pt1vp.z = constvp[2];
            pt2vp.z = constvp[2];
            pt1vs.z = constvs[2];
            pt2vs.z = constvs[2];
            pt1rho.z = constrho[2];
            pt2rho.z = constrho[2];
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              pt1_extra_param[ii].z = 0.0;
              pt2_extra_param[ii].z = 0.0;
            }
          } else {
            findPointZValue(i, j_out, k, pt1vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i, j_out, k, pt1vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i, j_out, k, pt1rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            findPointZValue(i + 1, j_out, k, pt2vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i + 1, j_out, k, pt2vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i + 1, j_out, k, pt2rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              findPointZValue(i, j_out, k, pt1_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
              findPointZValue(i + 1, j_out, k, pt2_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
            }
          }
          pt3vp.z = pt1vp.z;
          pt4vp.z = pt2vp.z;
          pt3vs.z = pt1vs.z;
          pt4vs.z = pt2vs.z;
          pt3rho.z = pt1rho.z;
          pt4rho.z = pt2rho.z;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt3_extra_param[ii].z = pt1_extra_param[ii].z;
            pt4_extra_param[ii].z = pt2_extra_param[ii].z;
          }

          bool triangulate_124 = true;
          NRLib::Point vec1, vec2;
          vec1 = pt1vp - pt2vp;
          vec1.z = 0;
          vec2 = pt4vp - pt2vp;
          vec2.z = 0;
          double delaunay_angle = vec1.GetAngle(vec2);
          vec1 = pt1vp - pt3vp;
          vec1.z = 0;
          vec2 = pt4vp - pt3vp;
          vec2.z = 0;
          delaunay_angle += vec1.GetAngle(vec2);
          if (delaunay_angle > NRLib::Pi) {
            triangulate_124 = false;
          }

          NRLib::Triangle triangle3(pt1vp, pt2vp, pt3vp);
          NRLib::Triangle triangle4(pt2vp, pt3vp, pt4vp);
          NRLib::Triangle triangle7(pt1vs, pt2vs, pt3vs);
          NRLib::Triangle triangle8(pt2vs, pt3vs, pt4vs);
          NRLib::Triangle triangle11(pt1rho, pt2rho, pt3rho);
          NRLib::Triangle triangle12(pt2rho, pt3rho, pt4rho);
          std::vector<NRLib::Triangle> triangles_extra_param;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            NRLib::Triangle triangle_temp(pt1_extra_param[ii], pt2_extra_param[ii], pt3_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
            triangle_temp.SetCornerPoints(pt2_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
          }
          if (triangulate_124) {
            triangle3.SetCornerPoints(pt1vp, pt2vp, pt4vp);
            triangle4.SetCornerPoints(pt1vp, pt3vp, pt4vp);
            triangle7.SetCornerPoints(pt1vs, pt2vs, pt4vs);
            triangle8.SetCornerPoints(pt1vs, pt3vs, pt4vs);
            triangle11.SetCornerPoints(pt1rho, pt2rho, pt4rho);
            triangle12.SetCornerPoints(pt1rho, pt3rho, pt4rho);
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              triangles_extra_param[ii * 2].SetCornerPoints(pt1_extra_param[ii], pt2_extra_param[ii], pt4_extra_param[ii]);
              triangles_extra_param[ii * 2 + 1].SetCornerPoints(pt1_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            }
          }
          x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
          y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
          x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
          y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
          x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
          y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
          x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
          y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

          cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
          cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
          cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
          cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

          start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
          start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
          end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
          end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
          if (end_ii > vpgrid.GetNI()) {
            end_ii = vpgrid.GetNI();
          }
          if (end_jj > vpgrid.GetNJ()) {
            end_jj = vpgrid.GetNJ();
          }
          NRLib::Polygon inside_e_cells;
          inside_e_cells.AddPoint(pt1vp);
          inside_e_cells.AddPoint(pt2vp);
          inside_e_cells.AddPoint(mid_edge2);
          // if(k > 0)
          //  inside_e_cells.AddPoint(geometry.FindCornerPoint(i,j,k-1,1,1,1));
          //  else
          if (k <= botk) {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i, j_out, k, 1, 1, 1)));
          } else {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i, j_out, k - 1, 1, 1, 1)));
          }
          inside_e_cells.AddPoint(mid_edge1);
          for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
            for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
              double x, y, z;
              vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
              //   if(toptime.IsMissing(toptime.GetZ(x,y))==false){
              NRLib::Point p1(x, y, 0.0);
              NRLib::Point p2(x, y, 1000.0);
              if (inside_e_cells.IsInsidePolygonXY(p1)) {
                NRLib::Line line(p1, p2, false, false);
                NRLib::Point intersec_pt;
                //bool intersect = triangle3.FindIntersection(line, intersec_pt, true);
                double dist = triangle3.FindNearestPoint(line, intersec_pt);
                if (dist < 0.00000000001) {
                  //if(intersect==true){
                  vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                  meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                  triangle7.FindIntersection(line, intersec_pt, true);
                  vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                  meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                  triangle11.FindIntersection(line, intersec_pt, true);
                  rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                  meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                  for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                    triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                    extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  }
                } else {
                  //intersect = triangle4.FindIntersection(line, intersec_pt, true);
                  dist = triangle4.FindNearestPoint(line, intersec_pt);
                  //if(intersect == true){
                  if (dist < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                    meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                    triangle8.FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                    meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                    triangle12.FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                    meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                    for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                      triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                      extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                }
                // if(intersect){
                if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                  if (dist < 0.00000000001) {
                    for (size_t l = 0; l < n_angle; l++) {
                      rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
                    }
                  }
                }
              }
            }
          }
        }
      }
    } //end for loop over j along i=0 layer



    i_out = 0;
    //------------------------Filling in i_out=egrid.GetNI()-1 right half edge------------------------------------

    for (size_t j = 0; j < egrid.GetNJ() - 1; j++) { //Filling in i_out=egrid.GetNI()-1 right half edge
      i_out = 0;
      //  if(geometry.IsPillarActive(i_out,j) && geometry.IsPillarActive(i_out,j+1) && geometry.IsPillarActive(i_out+1, j) && geometry.IsPillarActive(i_out+1, j+1) &&
      //          geometry.IsPillarActive(i_out,j+2) && geometry.IsPillarActive(i_out+1, j+2) ){
      while (i_out < static_cast<int>(egrid.GetNI()) && !(geometry.IsPillarActive(i_out, j) && geometry.IsPillarActive(i_out, j + 1) && geometry.IsPillarActive(i_out + 1, j) && geometry.IsPillarActive(i_out + 1, j + 1) &&
        geometry.IsPillarActive(i_out, j + 2) && geometry.IsPillarActive(i_out + 1, j + 2))) {
          i_out++;
      }
      if (i_out < static_cast<int>(egrid.GetNI())) {
        if (k <= botk) {
          pt1vp = geometry.FindCellCenterPoint(i_out, j, k);
          pt2vp = geometry.FindCellCenterPoint(i_out, j + 1, k);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i_out, j, k, 0, 0, 0) + geometry.FindCornerPoint(i_out, j, k, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j, k, 0, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k, 0, 0, 0) + geometry.FindCornerPoint(i_out, j + 1, k, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j + 1, k, 0, 1, 1));
        } else {
          pt1vp = geometry.FindCellCenterPoint(i_out, j, k - 1);
          pt2vp = geometry.FindCellCenterPoint(i_out, j + 1, k - 1);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 0, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i_out, j + 1, k - 1, 0, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j + 1, k - 1, 0, 1, 1));
        }
        //mid_edge2 = geometry.FindCornerPoint(i_out,j+1,k,1,0,0)+geometry.FindCornerPoint(i_out,j+1,k,1,1,0);


        int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
        int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
        pt3vp = mid_edge1 - pt1vp;
        int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
        mid_edge1 = 0.5 * mid_edge1;
        pt4vp = mid_edge2 - pt2vp;
        int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
        mid_edge2 = 0.5 * mid_edge2;
        if (inside1 || inside2 || inside3 || inside4) {
          pt1vs = pt1vp;
          pt1rho = pt1vp;
          pt2vs = pt2vp;
          pt2rho = pt2vp;
          pt3vs = pt3vp;
          pt3rho = pt3vp;
          pt4vs = pt4vp;
          pt4rho = pt4vp;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt1_extra_param[ii] = pt1vp;
            pt2_extra_param[ii] = pt2vp;
            pt3_extra_param[ii] = pt3vp;
            pt4_extra_param[ii] = pt4vp;
          }
          if (k == botk + 1) {
            pt1vp.z = constvp[2];
            pt2vp.z = constvp[2];
            pt1vs.z = constvs[2];
            pt2vs.z = constvs[2];
            pt1rho.z = constrho[2];
            pt2rho.z = constrho[2];
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              pt1_extra_param[ii].z = 0.0;
              pt2_extra_param[ii].z = 0.0;
            }
          } else {
            findPointZValue(i_out, j, k, pt1vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i_out, j, k, pt1vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i_out, j, k, pt1rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            findPointZValue(i_out, j + 1, k, pt2vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i_out, j + 1, k, pt2vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i_out, j + 1, k, pt2rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              findPointZValue(i_out, j, k, pt1_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
              findPointZValue(i_out, j + 1, k, pt2_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
            }
          }
          pt3vp.z = pt1vp.z;
          pt4vp.z = pt2vp.z;
          pt3vs.z = pt1vs.z;
          pt4vs.z = pt2vs.z;
          pt3rho.z = pt1rho.z;
          pt4rho.z = pt2rho.z;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt3_extra_param[ii].z = pt1_extra_param[ii].z;
            pt4_extra_param[ii].z = pt2_extra_param[ii].z;
          }

          bool triangulate_124 = true;
          NRLib::Point vec1, vec2;
          vec1 = pt1vp - pt2vp;
          vec1.z = 0;
          vec2 = pt4vp - pt2vp;
          vec2.z = 0;
          double delaunay_angle = vec1.GetAngle(vec2);
          vec1 = pt1vp - pt3vp;
          vec1.z = 0;
          vec2 = pt4vp - pt3vp;
          vec2.z = 0;
          delaunay_angle += vec1.GetAngle(vec2);
          if (delaunay_angle > NRLib::Pi) {
            triangulate_124 = false;
          }
          NRLib::Triangle triangle3(pt1vp, pt2vp, pt3vp);
          NRLib::Triangle triangle4(pt2vp, pt3vp, pt4vp);
          NRLib::Triangle triangle7(pt1vs, pt2vs, pt3vs);
          NRLib::Triangle triangle8(pt2vs, pt3vs, pt4vs);
          NRLib::Triangle triangle11(pt1rho, pt2rho, pt3rho);
          NRLib::Triangle triangle12(pt2rho, pt3rho, pt4rho);
          std::vector<NRLib::Triangle> triangles_extra_param;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            NRLib::Triangle triangle_temp(pt1_extra_param[ii], pt2_extra_param[ii], pt3_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
            triangle_temp.SetCornerPoints(pt2_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
          }
          if (triangulate_124) {
            triangle3.SetCornerPoints(pt1vp, pt2vp, pt4vp);
            triangle4.SetCornerPoints(pt1vp, pt3vp, pt4vp);
            triangle7.SetCornerPoints(pt1vs, pt2vs, pt4vs);
            triangle8.SetCornerPoints(pt1vs, pt3vs, pt4vs);
            triangle11.SetCornerPoints(pt1rho, pt2rho, pt4rho);
            triangle12.SetCornerPoints(pt1rho, pt3rho, pt4rho);
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              triangles_extra_param[ii * 2].SetCornerPoints(pt1_extra_param[ii], pt2_extra_param[ii], pt4_extra_param[ii]);
              triangles_extra_param[ii * 2 + 1].SetCornerPoints(pt1_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            }
          }
          x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
          y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
          x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
          y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
          x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
          y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
          x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
          y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

          cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
          cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
          cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
          cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

          start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
          start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
          end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
          end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
          if (end_ii > vpgrid.GetNI()) {
            end_ii = vpgrid.GetNI();
          }
          if (end_jj > vpgrid.GetNJ()) {
            end_jj = vpgrid.GetNJ();
          }
          NRLib::Polygon inside_e_cells;
          inside_e_cells.AddPoint(pt1vp);
          inside_e_cells.AddPoint(pt2vp);
          inside_e_cells.AddPoint(mid_edge2);
          //if(k > 0)
          //  inside_e_cells.AddPoint(geometry.FindCornerPoint(i,j,k-1,0,1,1));
          //else
          if (k <= botk) {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i_out, j, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j, k, 0, 1, 1)));
          } else {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 0, 1, 1)));
          }
          inside_e_cells.AddPoint(mid_edge1);
          for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
            for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
              double x, y, z;
              vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
              //   if(toptime.IsMissing(toptime.GetZ(x,y))==false){
              NRLib::Point p1(x, y, 0.0);
              NRLib::Point p2(x, y, 1000.0);
              if (inside_e_cells.IsInsidePolygonXY(p1)) {
                NRLib::Line line(p1, p2, false, false);
                NRLib::Point intersec_pt;
                //bool intersect = triangle3.FindIntersection(line, intersec_pt, true);
                double dist = triangle3.FindNearestPoint(line, intersec_pt);
                //if(intersect==true){
                if (dist < 0.00000000001) {
                  vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                  meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                  triangle7.FindIntersection(line, intersec_pt, true);
                  vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                  meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                  triangle11.FindIntersection(line, intersec_pt, true);
                  rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                  meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                  for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                    triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                    extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  }
                } else {
                  //intersect = triangle4.FindIntersection(line, intersec_pt, true);
                  dist = triangle4.FindNearestPoint(line, intersec_pt);
                  //if(intersect == true){
                  if (dist < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                    meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                    triangle8.FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                    meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                    triangle12.FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                    meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                    for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                      triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                      extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                }
                //if(intersect){
                if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                  if (dist < 0.00000000001) {
                    for (size_t l = 0; l < n_angle; l++) {
                      rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
                    }
                  }
                }
              }
            }
          }
        }
      }
    } //end for loop over j along i=0 layer

    i_out = static_cast<int>(egrid.GetNI()) - 1;
    //------------------------Filling in i_out=egrid.GetNI()-1 right half edge------------------------------------

    for (size_t j = 0; j < egrid.GetNJ() - 1; j++) { //Filling in i_out=egrid.GetNI()-1 right half edge
      i_out = static_cast<int>(egrid.GetNI()) - 1;
      // if(geometry.IsPillarActive(i_out,j) && geometry.IsPillarActive(i_out,j+1) && geometry.IsPillarActive(i_out+1, j) && geometry.IsPillarActive(i_out+1, j+1) &&
      //         geometry.IsPillarActive(i_out,j+2) && geometry.IsPillarActive(i_out+1, j+2) ){
      while (i_out >= 0 && !(geometry.IsPillarActive(i_out, j) && geometry.IsPillarActive(i_out, j + 1) && geometry.IsPillarActive(i_out + 1, j) && geometry.IsPillarActive(i_out + 1, j + 1) &&
        geometry.IsPillarActive(i_out, j + 2) && geometry.IsPillarActive(i_out + 1, j + 2))) {
          i_out--;
      }
      if (i_out >= 0) {
        if (k <= botk) {
          pt1vp = geometry.FindCellCenterPoint(i_out, j, k);
          pt2vp = geometry.FindCellCenterPoint(i_out, j + 1, k);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i_out, j, k, 1, 0, 0) + geometry.FindCornerPoint(i_out, j, k, 1, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j, k, 1, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k, 1, 0, 0) + geometry.FindCornerPoint(i_out, j + 1, k, 1, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j + 1, k, 1, 1, 1));
        } else {
          pt1vp = geometry.FindCellCenterPoint(i_out, j, k - 1);
          pt2vp = geometry.FindCellCenterPoint(i_out, j + 1, k - 1);
          mid_edge1 = 0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 1, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 1, 1, 1));
          mid_edge2 = 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i_out, j + 1, k - 1, 1, 0, 1)) + 0.5 * (geometry.FindCornerPoint(i_out, j + 1, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j + 1, k - 1, 1, 1, 1));
        }
        //mid_edge2 = geometry.FindCornerPoint(i_out,j+1,k,1,0,0)+geometry.FindCornerPoint(i_out,j+1,k,1,1,0);

        int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
        int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
        pt3vp = mid_edge1 - pt1vp;
        int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
        mid_edge1 = 0.5 * mid_edge1;
        pt4vp = mid_edge2 - pt2vp;
        int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
        mid_edge2 = 0.5 * mid_edge2;
        if (inside1 || inside2 || inside3 || inside4) {
          pt1vs = pt1vp;
          pt1rho = pt1vp;
          pt2vs = pt2vp;
          pt2rho = pt2vp;
          pt3vs = pt3vp;
          pt3rho = pt3vp;
          pt4vs = pt4vp;
          pt4rho = pt4vp;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt1_extra_param[ii] = pt1vp;
            pt2_extra_param[ii] = pt2vp;
            pt3_extra_param[ii] = pt3vp;
            pt4_extra_param[ii] = pt4vp;
          }
          if (k == botk + 1) {
            pt1vp.z = constvp[2];
            pt2vp.z = constvp[2];
            pt1vs.z = constvs[2];
            pt2vs.z = constvs[2];
            pt1rho.z = constrho[2];
            pt2rho.z = constrho[2];
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              pt1_extra_param[ii].z = 0.0;
              pt2_extra_param[ii].z = 0.0;
            }
          } else {
            findPointZValue(i_out, j, k, pt1vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i_out, j, k, pt1vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i_out, j, k, pt1rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            findPointZValue(i_out, j + 1, k, pt2vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
            findPointZValue(i_out, j + 1, k, pt2vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
            findPointZValue(i_out, j + 1, k, pt2rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              findPointZValue(i_out, j, k, pt1_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
              findPointZValue(i_out, j + 1, k, pt2_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
            }
          }
          pt3vp.z = pt1vp.z;
          pt4vp.z = pt2vp.z;
          pt3vs.z = pt1vs.z;
          pt4vs.z = pt2vs.z;
          pt3rho.z = pt1rho.z;
          pt4rho.z = pt2rho.z;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            pt3_extra_param[ii].z = pt1_extra_param[ii].z;
            pt4_extra_param[ii].z = pt2_extra_param[ii].z;
          }

          bool triangulate_124 = true;
          NRLib::Point vec1, vec2;
          vec1 = pt1vp - pt2vp;
          vec1.z = 0;
          vec2 = pt4vp - pt2vp;
          vec2.z = 0;
          double delaunay_angle = vec1.GetAngle(vec2);
          vec1 = pt1vp - pt3vp;
          vec1.z = 0;
          vec2 = pt4vp - pt3vp;
          vec2.z = 0;
          delaunay_angle += vec1.GetAngle(vec2);
          if (delaunay_angle > NRLib::Pi) {
            triangulate_124 = false;
          }
          NRLib::Triangle triangle3(pt1vp, pt2vp, pt3vp);
          NRLib::Triangle triangle4(pt2vp, pt3vp, pt4vp);
          NRLib::Triangle triangle7(pt1vs, pt2vs, pt3vs);
          NRLib::Triangle triangle8(pt2vs, pt3vs, pt4vs);
          NRLib::Triangle triangle11(pt1rho, pt2rho, pt3rho);
          NRLib::Triangle triangle12(pt2rho, pt3rho, pt4rho);
          std::vector<NRLib::Triangle> triangles_extra_param;
          for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
            NRLib::Triangle triangle_temp(pt1_extra_param[ii], pt2_extra_param[ii], pt3_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
            triangle_temp.SetCornerPoints(pt2_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            triangles_extra_param.push_back(triangle_temp);
          }
          if (triangulate_124) {
            triangle3.SetCornerPoints(pt1vp, pt2vp, pt4vp);
            triangle4.SetCornerPoints(pt1vp, pt3vp, pt4vp);
            triangle7.SetCornerPoints(pt1vs, pt2vs, pt4vs);
            triangle8.SetCornerPoints(pt1vs, pt3vs, pt4vs);
            triangle11.SetCornerPoints(pt1rho, pt2rho, pt4rho);
            triangle12.SetCornerPoints(pt1rho, pt3rho, pt4rho);
            for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
              triangles_extra_param[ii * 2].SetCornerPoints(pt1_extra_param[ii], pt2_extra_param[ii], pt4_extra_param[ii]);
              triangles_extra_param[ii * 2 + 1].SetCornerPoints(pt1_extra_param[ii], pt3_extra_param[ii], pt4_extra_param[ii]);
            }
          }
          x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
          y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
          x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
          y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
          x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
          y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
          x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
          y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

          cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
          cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
          cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
          cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

          start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
          start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
          end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
          end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
          if (end_ii > vpgrid.GetNI()) {
            end_ii = vpgrid.GetNI();
          }
          if (end_jj > vpgrid.GetNJ()) {
            end_jj = vpgrid.GetNJ();
          }
          NRLib::Polygon inside_e_cells;
          inside_e_cells.AddPoint(pt1vp);
          inside_e_cells.AddPoint(pt2vp);
          inside_e_cells.AddPoint(mid_edge2);
          // if( k > 0)
          //   inside_e_cells.AddPoint(geometry.FindCornerPoint(i,j,k-1,1,1,1));
          //  else
          if (k <= botk) {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i_out, j, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j, k, 1, 1, 1)));
          } else {
            inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i_out, j, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j, k - 1, 1, 1, 1)));
          }
          inside_e_cells.AddPoint(mid_edge1);
          for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
            for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
              double x, y, z;
              vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
              //       if(toptime.IsMissing(toptime.GetZ(x,y))==false){
              NRLib::Point p1(x, y, 0.0);
              NRLib::Point p2(x, y, 1000.0);
              if (inside_e_cells.IsInsidePolygonXY(p1)) {
                NRLib::Line line(p1, p2, false, false);
                NRLib::Point intersec_pt;
                //bool intersect = triangle3.FindIntersection(line, intersec_pt, true);
                double dist = triangle3.FindNearestPoint(line, intersec_pt);
                if (dist < 0.00000000001) {
                  //if(intersect==true){
                  vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                  meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                  triangle7.FindIntersection(line, intersec_pt, true);
                  vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                  meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                  triangle11.FindIntersection(line, intersec_pt, true);
                  rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                  meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                  for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                    triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                    extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                  }
                } else {
                  dist = triangle4.FindNearestPoint(line, intersec_pt);
                  //intersect = triangle4.FindIntersection(line, intersec_pt, true);
                  // if(intersect == true){
                  if (dist < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
                    meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
                    triangle8.FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
                    meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
                    triangle12.FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
                    meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
                    for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                      triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                      extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                }
                //if(intersect){
                if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
                  if (dist < 0.00000000001) {
                    for (size_t l = 0; l < n_angle; l++) {
                      rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
                    }
                  }
                }
              }
            }
          }
        }
      }
    } //end for loop over j along i=egrid.GetNI()-1 layer

    //------------------------Bot left corner------------------------------------

    //Bot left corner
    i_out = 0;
    j_out = 0;
    if (k <= botk) {
      pt1vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 1));
      pt2vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 1)) + pt1vp);
      pt4vp = geometry.FindCellCenterPoint(i_out, j_out, k);
      pt3vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 1)) + pt1vp);
    } else {
      pt1vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 1));
      pt2vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 1)) + pt1vp);
      pt4vp = geometry.FindCellCenterPoint(i_out, j_out, k - 1);
      pt3vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 1)) + pt1vp);
    }
    int inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
    int inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
    int inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
    int inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
    if (inside1 || inside2 || inside3 || inside4) {
      pt1vs = pt1vp;
      pt1rho = pt1vp;
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        pt1_extra_param[ii] = pt1vp;
      }
      if (k == botk + 1) {
        pt4vp.z = constvp[2];
        pt4vs.z = constvs[2];
        pt4rho.z = constrho[2];
        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          pt4_extra_param[ii].z = 0.0;
        }
      } else {
        findPointZValue(i_out, j_out, k, pt4vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          findPointZValue(i_out, j_out, k, pt4_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
        }
      }
      x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
      y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
      x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
      y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
      x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
      y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
      x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
      y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

      cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
      cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
      cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
      cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

      start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
      start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
      end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
      end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
      if (end_ii > vpgrid.GetNI()) {
        end_ii = vpgrid.GetNI();
      }
      if (end_jj > vpgrid.GetNJ()) {
        end_jj = vpgrid.GetNJ();
      }
      NRLib::Polygon inside_e_cells;
      inside_e_cells.AddPoint(pt1vp);
      inside_e_cells.AddPoint(pt2vp);
      inside_e_cells.AddPoint(pt4vp);
      inside_e_cells.AddPoint(pt3vp);
      for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
        for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
          double x, y, z;
          vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
          //   if(toptime.IsMissing(toptime.GetZ(x,y))==false){
          NRLib::Point p1(x, y, 0.0);
          if (inside_e_cells.IsInsidePolygonXY(p1)) {
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vp.z);
            diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
            meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vs.z);
            diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
            meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4rho.z);
            diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
            meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
            for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
              extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(pt4_extra_param[ii].z);
            }
            if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
              for (size_t l = 0; l < n_angle; l++) {
                rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
              }
            }
          }
        }
      }
    }


    //------------------------Bot right corner------------------------------------

    //Bot right corner
    i_out = static_cast<int>(egrid.GetNI()) - 1;
    j_out = 0;
    if (k <= botk) {
      pt2vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 1));
      pt1vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 1)) + pt2vp);
      pt3vp = geometry.FindCellCenterPoint(i_out, j_out, k);
      pt4vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 1)) + pt2vp);
    } else {
      pt2vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 1));
      pt1vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 1)) + pt2vp);
      pt3vp = geometry.FindCellCenterPoint(i_out, j_out, k - 1);
      pt4vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 1)) + pt2vp);
    }
    inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
    inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
    inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
    inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
    if (inside1 || inside2 || inside3 || inside4) {
      pt1vs = pt1vp;
      pt1rho = pt1vp;
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        pt1_extra_param[ii] = pt1vp;
      }
      if (k == botk + 1) {
        pt4vp.z = constvp[2];
        pt4vs.z = constvs[2];
        pt4rho.z = constrho[2];
        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          pt4_extra_param[ii].z = 0.0;
        }
      } else {
        findPointZValue(i_out, j_out, k, pt4vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          findPointZValue(i_out, j_out, k, pt4_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
        }
      }
      x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
      y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
      x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
      y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
      x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
      y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
      x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
      y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

      cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
      cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
      cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
      cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

      start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
      start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
      end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
      end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
      if (end_ii > vpgrid.GetNI()) {
        end_ii = vpgrid.GetNI();
      }
      if (end_jj > vpgrid.GetNJ()) {
        end_jj = vpgrid.GetNJ();
      }
      NRLib::Polygon inside_e_cells;
      inside_e_cells.AddPoint(pt1vp);
      inside_e_cells.AddPoint(pt2vp);
      inside_e_cells.AddPoint(pt4vp);
      inside_e_cells.AddPoint(pt3vp);
      for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
        for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
          double x, y, z;
          vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
          //    if(toptime.IsMissing(toptime.GetZ(x,y))==false){
          NRLib::Point p1(x, y, 0.0);
          if (inside_e_cells.IsInsidePolygonXY(p1)) {
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vp.z);
            diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
            meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vs.z);
            diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
            meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4rho.z);
            diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
            meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
            for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
              extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(pt4_extra_param[ii].z);
            }
            if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
              for (size_t l = 0; l < n_angle; l++) {
                rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
              }
            }
          }
        }
      }
    }

    //------------------------Top right corner------------------------------------

    //Top right corner
    i_out = static_cast<int>(egrid.GetNI()) - 1;
    j_out = static_cast<int>(egrid.GetNJ()) - 1;
    if (k <= botk) {
      pt4vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 1));
      pt2vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 0, 1)) + pt4vp);
      pt1vp = geometry.FindCellCenterPoint(i_out, j_out, k);
      pt3vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 1)) + pt4vp);
    } else {
      pt4vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 1));
      pt2vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 0, 1)) + pt4vp);
      pt1vp = geometry.FindCellCenterPoint(i_out, j_out, k - 1);
      pt3vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 1)) + pt4vp);
    }
    inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
    inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
    inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
    inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
    if (inside1 || inside2 || inside3 || inside4) {
      pt1vs = pt1vp;
      pt1rho = pt1vp;
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        pt1_extra_param[ii] = pt1vp;
      }
      if (k == botk + 1) {
        pt4vp.z = constvp[2];
        pt4vs.z = constvs[2];
        pt4rho.z = constrho[2];
        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          pt4_extra_param[ii].z = 0.0;
        }
      } else {
        findPointZValue(i_out, j_out, k, pt4vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          findPointZValue(i_out, j_out, k, pt4_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
        }
      }
      x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
      y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
      x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
      y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
      x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
      y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
      x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
      y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

      cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
      cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
      cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
      cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

      start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
      start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
      end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
      end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
      if (end_ii > vpgrid.GetNI()) {
        end_ii = vpgrid.GetNI();
      }
      if (end_jj > vpgrid.GetNJ()) {
        end_jj = vpgrid.GetNJ();
      }
      NRLib::Polygon inside_e_cells;
      inside_e_cells.AddPoint(pt1vp);
      inside_e_cells.AddPoint(pt2vp);
      inside_e_cells.AddPoint(pt4vp);
      inside_e_cells.AddPoint(pt3vp);
      for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
        for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
          double x, y, z;
          vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
          //  if(toptime.IsMissing(toptime.GetZ(x,y))==false){
          NRLib::Point p1(x, y, 0.0);
          if (inside_e_cells.IsInsidePolygonXY(p1)) {
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vp.z);
            diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
            meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vs.z);
            diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
            meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4rho.z);
            diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
            meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
            for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
              extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(pt4_extra_param[ii].z);
            }
            if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
              for (size_t l = 0; l < n_angle; l++) {
                rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
              }
            }
          }
        }
      }
    }

    //------------------------Top left corner------------------------------------

    //Top left corner
    i_out = 0;
    j_out = static_cast<int>(egrid.GetNJ()) - 1;
    if (k <= botk) {
      pt3vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 1, 1));
      pt1vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k, 0, 0, 1)) + pt3vp);
      pt2vp = geometry.FindCellCenterPoint(i_out, j_out, k);
      pt4vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k, 1, 1, 1)) + pt3vp);
    } else {
      pt3vp = 0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 1, 1));
      pt1vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 0, 0, 1)) + pt3vp);
      pt2vp = geometry.FindCellCenterPoint(i_out, j_out, k - 1);
      pt4vp = 0.5 * (0.5 * (geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 0) + geometry.FindCornerPoint(i_out, j_out, k - 1, 1, 1, 1)) + pt3vp);
    }
    inside1 = vpgrid.IsInside(pt1vp.x, pt1vp.y);
    inside2 = vpgrid.IsInside(pt2vp.x, pt2vp.y);
    inside3 = vpgrid.IsInside(pt3vp.x, pt3vp.y);
    inside4 = vpgrid.IsInside(pt4vp.x, pt4vp.y);
    if (inside1 || inside2 || inside3 || inside4) {
      pt1vs = pt1vp;
      pt1rho = pt1vp;
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        pt1_extra_param[ii] = pt1vp;
      }
      if (k == botk + 1) {
        pt4vp.z = constvp[2];
        pt4vs.z = constvs[2];
        pt4rho.z = constrho[2];
        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          pt4_extra_param[ii].z = 0.0;
        }
      } else {
        findPointZValue(i_out, j_out, k, pt4vp, geometry, vp_grid, value_above_vp, constvp[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4vs, geometry, vs_grid, value_above_vs, constvs[1], zlimit);
        findPointZValue(i_out, j_out, k, pt4rho, geometry, rho_grid, value_above_rho, constrho[1], zlimit);

        for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
          findPointZValue(i_out, j_out, k, pt4_extra_param[ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
        }
      }
      x1_rot = pt1vp.x * cosvpangle + pt1vp.y * sinvpangle;
      y1_rot = pt1vp.y * cosvpangle - pt1vp.x * sinvpangle;
      x2_rot = pt2vp.x * cosvpangle + pt2vp.y * sinvpangle;
      y2_rot = pt2vp.y * cosvpangle - pt2vp.x * sinvpangle;
      x3_rot = pt3vp.x * cosvpangle + pt3vp.y * sinvpangle;
      y3_rot = pt3vp.y * cosvpangle - pt3vp.x * sinvpangle;
      x4_rot = pt4vp.x * cosvpangle + pt4vp.y * sinvpangle;
      y4_rot = pt4vp.y * cosvpangle - pt4vp.x * sinvpangle;

      cell_min_x = min(min(x1_rot, x2_rot), min(x3_rot, x4_rot));
      cell_min_y = min(min(y1_rot, y2_rot), min(y3_rot, y4_rot));
      cell_max_x = max(max(x1_rot, x2_rot), max(x3_rot, x4_rot));
      cell_max_y = max(max(y1_rot, y2_rot), max(y3_rot, y4_rot));

      start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
      start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
      end_ii = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
      end_jj = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
      if (end_ii > vpgrid.GetNI()) {
        end_ii = vpgrid.GetNI();
      }
      if (end_jj > vpgrid.GetNJ()) {
        end_jj = vpgrid.GetNJ();
      }
      NRLib::Polygon inside_e_cells;
      inside_e_cells.AddPoint(pt1vp);
      inside_e_cells.AddPoint(pt2vp);
      inside_e_cells.AddPoint(pt4vp);
      inside_e_cells.AddPoint(pt3vp);
      for (size_t ii = start_ii; ii < end_ii; ii++) { //earlier: from 0 up to vpgrid.GetNI()
        for (size_t jj = start_jj; jj < end_jj; jj++) { //earlier: from 0 up to vpgrid.GetNJ()
          double x, y, z;
          vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
          //  if(toptime.IsMissing(toptime.GetZ(x,y))==false){
          NRLib::Point p1(x, y, 0.0);
          if (inside_e_cells.IsInsidePolygonXY(p1)) {
            //what? diff and mean commented out by me and Espen!!!
            vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vp.z);
            diffvp = vpgrid(ii, jj, (k - topk) + 1) - vpgrid(ii, jj, (k - topk));
            meanvp = 0.5 * (vpgrid(ii, jj, (k - topk) + 1) + vpgrid(ii, jj, (k - topk)));
            vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4vs.z);
            diffvs = vsgrid(ii, jj, (k - topk) + 1) - vsgrid(ii, jj, (k - topk));
            meanvs = 0.5 * (vsgrid(ii, jj, (k - topk) + 1) + vsgrid(ii, jj, (k - topk)));
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(pt4rho.z);
            diffrho = rhogrid(ii, jj, (k - topk) + 1) - rhogrid(ii, jj, (k - topk));
            meanrho = 0.5 * (rhogrid(ii, jj, (k - topk) + 1) + rhogrid(ii, jj, (k - topk)));
            for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
              extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(pt4_extra_param[ii].z);
            }
            if (seismic_parameters.modelSettings()->GetNMOCorr() == false){
              for (size_t l = 0; l < n_angle; l++) {
                rgrid[l](ii, jj, k - topk) = static_cast<float>(zoeppritz[l]->GetReflection(diffvp, meanvp, diffrho, meanrho, diffvs, meanvs));
              }
            }
          }
        }
      }
    }
  }
  for (size_t i = 0; i < n_angle; i++) {
    delete zoeppritz[i];
  }
}


void SeismicRegridding::findPointZValue(size_t i, size_t j, size_t k,
  NRLib::Point &point,
  const NRLib::EclipseGeometry &geometry,
  const NRLib::Grid<double> &grid,
  const NRLib::Grid2D<double> &value_above,
  double &default_value,
  double &zlimit) {

    if (geometry.IsActive(i, j, k)) {
      point.z = grid(i, j, k);
    } else if (geometry.GetDZ(i, j, k) < zlimit) {
      point.z = value_above(i, j);
    } else {
      point.z = default_value;
    }
}

void SeismicRegridding::addNoiseToReflections(unsigned long seed, double std_dev, std::vector<NRLib::StormContGrid> &grid_vec) {
  NRLib::Random::Initialize(seed);
  NRLib::Normal normal_distibrution(0, std_dev);

  size_t nx = grid_vec[0].GetNI();
  size_t ny = grid_vec[0].GetNJ();
  size_t nz = grid_vec[0].GetNK();
  for (size_t angle = 0; angle < grid_vec.size(); ++angle) {
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ny; ++j) {
        for (size_t k = 0; k < nz; ++k) {
          grid_vec[angle](i, j, k) += static_cast<float>(normal_distibrution.Draw());
        }
      }
    }
  }
}

void SeismicRegridding::addNoiseToReflectionsPos(unsigned long         seed, 
                                                 double                std_dev, 
                                                 NRLib::Grid2D<double> &refl) {
  NRLib::Random::Initialize(seed);
  NRLib::Normal normal_distibrution(0, std_dev);

  for (size_t i = 0; i < refl.GetNI(); ++i) {
    for (size_t j = 0; j < refl.GetNJ(); ++j) {
      refl(i, j) += static_cast<float>(normal_distibrution.Draw());
    }
  }
}




void SeismicRegridding::findTWT(NRLib::StormContGrid &vpgrid, NRLib::StormContGrid &vsgrid, NRLib::StormContGrid &twtgrid, NRLib::StormContGrid &zgrid,NRLib::RegularSurface<double> &toptime, NRLib::RegularSurface<double> &bottime, bool ps_seismic) {

  size_t nk = twtgrid.GetNK();
  double dx1 = vpgrid.GetDX();
  double dy1 = vpgrid.GetDY();
  double dx2 = bottime.GetDX();
  double dy2 = bottime.GetDY();
  for (size_t i = 0; i < vpgrid.GetNI(); i++) {
    for (size_t j = 0; j < vpgrid.GetNJ(); j++) {
      double x, y, z;
      vpgrid.FindCenterOfCell(i, j, 0, x, y, z);
      twtgrid(i, j, 0) = static_cast<float>(toptime.GetZ(x, y));
      if (toptime.IsMissing(twtgrid(i, j, 0)) == false) {
        for (size_t k = 1; k < nk; k++) {
          if(ps_seismic) {
            // twtgrid(i, j, k) = twtgrid(i, j, k-1) + static_cast<float>(2000.0*(zgrid(i, j, k)-zgrid(i, j, k-1))/(0.5*(vpgrid(i, j, 2*k)+vpgrid(i, j, (2*k-1)))));
            twtgrid(i, j, k) = twtgrid(i, j, k - 1) + static_cast<float>(1000.0 * (zgrid(i, j, k) - zgrid(i, j, k - 1)) / vpgrid(i, j, k + 1)) + static_cast<float>(1000.0 * (zgrid(i, j, k) - zgrid(i, j, k - 1)) / vsgrid(i, j, k + 1));
          } else {
            // twtgrid(i, j, k) = twtgrid(i, j, k-1) + static_cast<float>(2000.0*(zgrid(i, j, k)-zgrid(i, j, k-1))/(0.5*(vpgrid(i, j, 2*k)+vpgrid(i, j, (2*k-1)))));
            twtgrid(i, j, k) = twtgrid(i, j, k - 1) + static_cast<float>(2000.0 * (zgrid(i, j, k) - zgrid(i, j, k - 1)) / vpgrid(i, j, k + 1));
          }
        }

        size_t ii, jj;
        double xstart, xend;
        double ystart, yend;
        xstart = x - dx1;
        ystart = y - dy1;
        xend = x + dx1;
        yend = y + dy1;
        x = xstart;
        y = ystart;
        while (x < xend) {
          y = ystart;
          while (y < yend) {
            bottime.FindIndex(x, y, ii, jj);
            bottime(ii, jj) = twtgrid(i, j, nk - 1);
            y = y + dy2;
          }
          x = x + dx2;
        }
      } else {
        for (size_t k = 0; k < nk; k++) {
          twtgrid(i, j, k) = -999.0;
        }
      }
    }
  }
}
