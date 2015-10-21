#include <seismic_regridding.hpp>

#include <nrlib/geometry/triangle.hpp>
#include <nrlib/eclipsegrid/eclipsegrid.hpp>

#include <physics/zoeppritz.hpp>
#include <physics/zoeppritz_ps.hpp>
#include <physics/zoeppritz_pp.hpp>
#include <physics/wavelet.hpp>
#include <nrlib/random/random.hpp>
#include <nrlib/random/normal.hpp>
#include <ctime>
#include <seismic_geometry.hpp>

void SeismicRegridding::MakeSeismicRegridding(SeismicParameters &seismic_parameters) {
  printf("Start finding Zvalues.\n");
  FindZValues(seismic_parameters);
  printf("Zvalues found.\n");

  printf("Start finding elastic parameters.\n");
  FindVp(seismic_parameters);
  printf("Elastic parameters found.\n");

  seismic_parameters.deleteEclipseGrid();

  NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
  NRLib::RegularSurface<double> &bottime = seismic_parameters.bottomTime();
  NRLib::StormContGrid &vpgrid           = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid           = seismic_parameters.vsGrid();
  NRLib::StormContGrid &twtgrid          = seismic_parameters.twtGrid();
  NRLib::StormContGrid &zgrid            = seismic_parameters.zGrid();
  bool find_for_ps                       = seismic_parameters.modelSettings()->GetPSSeismic();
  std::vector<double> constvp            = seismic_parameters.modelSettings()->GetConstVp();
  Wavelet* wavelet                       = seismic_parameters.wavelet();

  //find twt grid
  FindTWT(vpgrid, vsgrid, twtgrid, zgrid, toptime, bottime, find_for_ps);

  //generate, write and delete vrms grid if writing is requested
  if (seismic_parameters.modelSettings()->GetNMOCorr() && seismic_parameters.modelSettings()->GetOutputVrms()){
    FindVrms(seismic_parameters);
    seismic_parameters.seismicOutput()->writeVrms(seismic_parameters);
    seismic_parameters.deleteVrmsGrid();
  }

  //add wavelet above and below toptime and bottime
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

  //write toptime and bottime
  if (seismic_parameters.modelSettings()->GetOutputTimeSurfaces()) {
    seismic_parameters.seismicOutput()->writeTimeSurfaces(seismic_parameters);
  }

  //resample, write and delete extra parameter grids
  if (seismic_parameters.modelSettings()->GetOutputExtraParametersTimeSegy()) {
    seismic_parameters.seismicOutput()->writeExtraParametersTimeSegy(seismic_parameters);
  }
  if (seismic_parameters.modelSettings()->GetOutputExtraParametersDepthSegy()) {
    seismic_parameters.seismicOutput()->writeExtraParametersDepthSegy(seismic_parameters);
  }
  seismic_parameters.deleteExtraParameterGrids();
 
  //resample and write elastic parameters in segy
  if (seismic_parameters.modelSettings()->GetOutputElasticParametersTimeSegy()) {
    seismic_parameters.seismicOutput()->writeElasticParametersTimeSegy(seismic_parameters);
  }
  if (seismic_parameters.modelSettings()->GetOutputElasticParametersDepthSegy()) {
    seismic_parameters.seismicOutput()->writeElasticParametersDepthSegy(seismic_parameters);
  }

  //write elastic parameters, z values and twt on storm format
  if (seismic_parameters.modelSettings()->GetOutputVp()) {
    seismic_parameters.seismicOutput()->writeVpVsRho(seismic_parameters);
  }
  if (seismic_parameters.modelSettings()->GetOutputZvalues()) {
    seismic_parameters.seismicOutput()->writeZValues(seismic_parameters);
  }
  if (seismic_parameters.modelSettings()->GetOutputTwt()) {
    seismic_parameters.seismicOutput()->writeTwt(seismic_parameters);
  }
}


void SeismicRegridding::FindZValues(SeismicParameters &seismic_parameters) {
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

void SeismicRegridding::FindVrms(SeismicParameters &seismic_parameters){
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


void SeismicRegridding::FindPointZValue(size_t i, size_t j, size_t k,
                                        NRLib::Point &point,
                                        const NRLib::EclipseGeometry &geometry,
                                        const NRLib::Grid<double> &grid,
                                        const NRLib::Grid2D<double> &value_above,
                                        double &default_value,
                                        double &zlimit) 
{

    if (geometry.IsActive(i, j, k)) {
      point.z = grid(i, j, k);
    } else if (geometry.GetDZ(i, j, k) < zlimit) {
      point.z = value_above(i, j);
    } else {
      point.z = default_value;
    }
}


void SeismicRegridding::AddNoiseToReflectionsPos(unsigned long         seed, 
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




void SeismicRegridding::FindTWT(NRLib::StormContGrid &vpgrid, NRLib::StormContGrid &vsgrid, NRLib::StormContGrid &twtgrid, NRLib::StormContGrid &zgrid,NRLib::RegularSurface<double> &toptime, NRLib::RegularSurface<double> &bottime, bool ps_seismic) {

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

void SeismicRegridding::FindVp(SeismicParameters &seismic_parameters)
{
  NRLib::StormContGrid              &vpgrid               = seismic_parameters.vpGrid();
  NRLib::StormContGrid              &vsgrid               = seismic_parameters.vsGrid();
  NRLib::StormContGrid              &rhogrid              = seismic_parameters.rhoGrid();
  std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();
  const NRLib::EclipseGrid          &egrid                = seismic_parameters.eclipseGrid();
  const NRLib::EclipseGeometry      &geometry             = egrid.GetGeometry();

  size_t topk   = seismic_parameters.topK();
  size_t botk   = seismic_parameters.bottomK();
  double zlimit = seismic_parameters.modelSettings()->GetZeroThicknessLimit();

  std::vector<double> constvp                        = seismic_parameters.modelSettings()->GetConstVp();
  std::vector<double> constvs                        = seismic_parameters.modelSettings()->GetConstVs();
  std::vector<double> constrho                       = seismic_parameters.modelSettings()->GetConstRho();
  std::vector<std::string> names                     = seismic_parameters.modelSettings()->GetParameterNames();
  std::vector<double> extra_parameter_default_values = seismic_parameters.modelSettings()->GetExtraParameterDefaultValues();
  std::vector<std::string> extra_parameter_names     = seismic_parameters.modelSettings()->GetExtraParameterNames();

  const NRLib::Grid<double> &vp_grid  = egrid.GetParameter(names[0]);
  const NRLib::Grid<double> &vs_grid  = egrid.GetParameter(names[1]);
  const NRLib::Grid<double> &rho_grid = egrid.GetParameter(names[2]);
  std::vector<NRLib::Grid<double> > parameter_grid_from_eclipse;
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    const NRLib::Grid<double> &one_parameter_grid = egrid.GetParameter(extra_parameter_names[i]);
    parameter_grid_from_eclipse.push_back(one_parameter_grid);
  }

  NRLib::Grid2D<double> value_above_vp (egrid.GetNI(), egrid.GetNJ(), constvp[0]);
  NRLib::Grid2D<double> value_above_vs (egrid.GetNI(), egrid.GetNJ(), constvs[0]);
  NRLib::Grid2D<double> value_above_rho(egrid.GetNI(), egrid.GetNJ(), constrho[0]);
  std::vector<NRLib::Grid2D<double> > value_above_extra_param;
  for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
    NRLib::Grid2D<double> value_above_grid(egrid.GetNI(), egrid.GetNJ(), 0.0);
    value_above_extra_param.push_back(value_above_grid);
  }

  double vp_angle   = vpgrid.GetAngle();
  double cosvpangle = cos(vp_angle);
  double sinvpangle = sin(vp_angle);
  double x_min_rot  = vpgrid.GetXMin() * cos(vp_angle) + vpgrid.GetYMin() * sin(vp_angle);
  double y_min_rot  = vpgrid.GetYMin() * cos(vp_angle) - vpgrid.GetXMin() * sin(vp_angle);

  double                    cell_min_x, cell_max_x, cell_min_y, cell_max_y;
  size_t                    start_ii, start_jj, end_ii, end_jj;
  std::vector<double>       x_rot(4), y_rot(4);
  std::vector<bool>         inside(4);
  std::vector<NRLib::Point> pt_vp(4), pt_vs(4), pt_rho(4);
  std::vector<std::vector<NRLib::Point> > pt_extra_param(extra_parameter_names.size());
  for (size_t i = 0; i < extra_parameter_names.size(); ++i)
    pt_extra_param[i] = pt_vp;

  for (size_t i = 0; i < vpgrid.GetNI(); i++) {
    for (size_t j = 0; j < vpgrid.GetNJ(); j++) {
      vpgrid(i, j, 0)  = static_cast<float>(constvp[0]);
      vsgrid(i, j, 0)  = static_cast<float>(constvs[0]);
      rhogrid(i, j, 0) = static_cast<float>(constrho[0]);
      for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
        extra_parameter_grid[ii](i, j, 0) = 0.0;
      }
    }
  }
  for (size_t k = topk; k <= botk + 1; k++) {
    for (size_t i = 0; i < egrid.GetNI() - 1; i++) {
      for (size_t j = 0; j < egrid.GetNJ() - 1; j++) {
        if (geometry.IsPillarActive(i, j)   && geometry.IsPillarActive(i + 1, j)     && geometry.IsPillarActive(i, j + 1) && geometry.IsPillarActive(i + 1, j + 1) &&
          geometry.IsPillarActive(i + 2, j) && geometry.IsPillarActive(i + 2, j + 1) && geometry.IsPillarActive(i, j + 2) && geometry.IsPillarActive(i + 1, j + 2) &&
          geometry.IsPillarActive(i + 2, j + 2)) {
          if (k <= botk) {
            for (size_t pt = 0; pt < 4; ++pt)
              pt_vp[pt] = geometry.FindCellCenterPoint(i+int(pt%2), j+int(floor(double(pt)/2)), k);
          }
          else {
            for (size_t pt = 0; pt < 4; ++pt)
              pt_vp[pt] = geometry.FindCellCenterPoint(i+int(pt%2), j+int(floor(double(pt)/2)), k - 1);
          }
          for (size_t pt = 0; pt < 4;++pt)
            inside[pt] = vpgrid.IsInside(pt_vp[pt].x, pt_vp[pt].y);
          if (inside[0] || inside[1] || inside[2] || inside[3]) {
            for (size_t pt = 0; pt < 4; ++pt){
              pt_vs[pt]  = pt_vp[pt];
              pt_rho[pt] = pt_vp[pt];
              for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                pt_extra_param[ii][pt] = pt_vp[pt];
              }
            }
            if (k == botk + 1) {
              for (size_t pt = 0; pt < 4; ++pt) {
                pt_vp[pt].z  = constvp[2];
                pt_vs[pt].z  = constvs[2];
                pt_rho[pt].z = constrho[2];
                for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                  pt_extra_param[ii][pt].z = 0.0;
                }
              }
            }
            else {
              for (size_t pt = 0; pt < 4; ++pt){
                FindPointZValue(i+int(pt%2), j+int(floor(double(pt/2))), k, pt_vp[pt],  geometry, vp_grid,  value_above_vp,  constvp[1],  zlimit);
                FindPointZValue(i+int(pt%2), j+int(floor(double(pt/2))), k, pt_vs[pt],  geometry, vs_grid,  value_above_vs,  constvs[1],  zlimit);
                FindPointZValue(i+int(pt%2), j+int(floor(double(pt/2))), k, pt_rho[pt], geometry, rho_grid, value_above_rho, constrho[1], zlimit);
                for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                  FindPointZValue(i+int(pt%2), j+int(floor(double(pt/2))), k, pt_extra_param[pt][ii], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
                }
              }
            }
            for (size_t pt = 0; pt < 4; ++pt) {
              value_above_vp (i+int(pt%2), j+int(floor(double(pt/2)))) = pt_vp[pt].z;
              value_above_vs (i+int(pt%2), j+int(floor(double(pt/2)))) = pt_vs[pt].z;
              value_above_rho(i+int(pt%2), j+int(floor(double(pt/2)))) = pt_rho[pt].z;
              for (size_t ii = 0; ii < extra_parameter_names.size(); ++ii) {
                value_above_extra_param[ii](i+int(pt%2), j+int(floor(double(pt/2)))) = pt_extra_param[ii][pt].z;
              }
            }

            bool triangulate_124 = Is124Triangulate(pt_vp);

            std::vector<NRLib::Triangle> triangles_elastic(6);
            std::vector<NRLib::Triangle> triangles_extra_param(extra_parameter_names.size()*2);
            SetElasticTriangles(pt_vp, pt_vs, pt_rho, pt_extra_param, triangulate_124, triangles_elastic, triangles_extra_param);

            for (size_t pt = 0; pt < 4; ++pt) {
              x_rot[pt] = pt_vp[pt].x * cosvpangle + pt_vp[pt].y *sinvpangle;
              y_rot[pt] = pt_vp[pt].y * cosvpangle - pt_vp[pt].x *sinvpangle;
            }

            cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
            cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
            cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
            cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

            start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 0.5));
            start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 0.5));
            end_ii   = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 1.0));
            end_jj   = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 1.0));
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
                if (true) { 
                  NRLib::Point p1, p2;
                  p1.x  = x;
                  p1.y  = y;
                  p1.z  = pt_vp[0].z;
                  p2    = p1;
                  p2.z += 1000;
                  NRLib::Line line(p1, p2, false, false);
                  NRLib::Point intersec_pt;
                  double dist = triangles_elastic[0].FindNearestPoint(line, intersec_pt); // To avoid numerical instabilities when point is on edge of triangle
                  if (dist < 0.00000000001) {
                    vpgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
                    triangles_elastic[2].FindIntersection(line, intersec_pt, true);
                    vsgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
                    triangles_elastic[4].FindIntersection(line, intersec_pt, true);
                    rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                      triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
                      extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                    }
                  }
                  else {
                    dist = triangles_elastic[1].FindNearestPoint(line, intersec_pt);
                    if (dist < 0.00000000001) {
                      vpgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      triangles_elastic[3].FindIntersection(line, intersec_pt, true);
                      vsgrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      triangles_elastic[5].FindIntersection(line, intersec_pt, true);
                      rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
                      for (size_t iii = 0; iii < extra_parameter_names.size(); ++iii) {
                        triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                        extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
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
    //-------------find edges---------------------
    for (size_t i = 0; i < egrid.GetNI() - 1; i++) {
      //bot edge
      size_t j = 0;
      if (FindBotCell(geometry, egrid.GetNJ(), i, j)){
        FindVpEdges(geometry,
                    extra_parameter_names.size(),
                    seismic_parameters,
                    value_above_vp,
                    value_above_vs,
                    value_above_rho,
                    value_above_extra_param,
                    vp_grid,
                    vs_grid,
                    rho_grid,
                    parameter_grid_from_eclipse,
                    i, j, k,
                    false, true, false, false);
      }
      //top edge
      j = egrid.GetNJ() - 1;
      if (FindTopCell(geometry, i, j)) {
        FindVpEdges(geometry,
                    extra_parameter_names.size(),
                    seismic_parameters,
                    value_above_vp,
                    value_above_vs,
                    value_above_rho,
                    value_above_extra_param,
                    vp_grid,
                    vs_grid,
                    rho_grid,
                    parameter_grid_from_eclipse,
                    i, j, k,
                    true, false, false, false);
      }
    }
    for (size_t j = 0; j < egrid.GetNJ() - 1; ++j) {
      //left edge
      size_t i = 0;
      if (FindLeftCell(geometry, egrid.GetNI(), i, j)) {
        FindVpEdges(geometry,
                    extra_parameter_names.size(),
                    seismic_parameters,
                    value_above_vp,
                    value_above_vs,
                    value_above_rho,
                    value_above_extra_param,
                    vp_grid,
                    vs_grid,
                    rho_grid,
                    parameter_grid_from_eclipse,
                    i, j, k,
                    false, false, false, true);
      }
      //right edge
      i = egrid.GetNI() - 1;
      if (FindRightCell(geometry, i, j)) {
        FindVpEdges(geometry,
                    extra_parameter_names.size(),
                    seismic_parameters,
                    value_above_vp,
                    value_above_vs,
                    value_above_rho,
                    value_above_extra_param,
                    vp_grid,
                    vs_grid,
                    rho_grid,
                    parameter_grid_from_eclipse,
                    i, j, k,
                    false, false, true, false);
      }
    }
    //-------------find corners---------------------
    //bot left
    size_t i = 0;
    size_t j = 0;
    FindCornerCellPoints(geometry,
                         pt_vp,
                         i,
                         j,
                         k,
                         botk);
    FindVpCorners(geometry,
                  extra_parameter_names.size(),
                  seismic_parameters,
                  value_above_vp,
                  value_above_vs,
                  value_above_rho,
                  value_above_extra_param,
                  vp_grid,
                  vs_grid,
                  rho_grid,
                  parameter_grid_from_eclipse,
                  i, j, k, pt_vp);
    //top left
    j = egrid.GetNJ() - 1;
    FindCornerCellPoints(geometry,
                         pt_vp,
                         i,
                         j,
                         k,
                         botk);
    FindVpCorners(geometry,
                  extra_parameter_names.size(),
                  seismic_parameters,
                  value_above_vp,
                  value_above_vs,
                  value_above_rho,
                  value_above_extra_param,
                  vp_grid,
                  vs_grid,
                  rho_grid,
                  parameter_grid_from_eclipse,
                  i, j, k, pt_vp);
    //top right
    i = egrid.GetNI() - 1;
    FindCornerCellPoints(geometry,
                         pt_vp,
                         i,
                         j,
                         k,
                         botk);
    FindVpCorners(geometry,
                  extra_parameter_names.size(),
                  seismic_parameters,
                  value_above_vp,
                  value_above_vs,
                  value_above_rho,
                  value_above_extra_param,
                  vp_grid,
                  vs_grid,
                  rho_grid,
                  parameter_grid_from_eclipse,
                  i, j, k, pt_vp);
    //bot right
    j = 0;
    FindCornerCellPoints(geometry,
                         pt_vp,
                         i,
                         j,
                         k,
                         botk);
    FindVpCorners(geometry,
                  extra_parameter_names.size(),
                  seismic_parameters,
                  value_above_vp,
                  value_above_vs,
                  value_above_rho,
                  value_above_extra_param,
                  vp_grid,
                  vs_grid,
                  rho_grid,
                  parameter_grid_from_eclipse,
                  i, j, k, pt_vp);
  }
}

void SeismicRegridding::FindVpEdges(const NRLib::EclipseGeometry        &geometry,
                                    size_t                               n_extra_param,
                                    SeismicParameters                   &seismic_parameters,
                                    NRLib::Grid2D<double>               &value_above_vp,
                                    NRLib::Grid2D<double>               &value_above_vs,
                                    NRLib::Grid2D<double>               &value_above_rho,
                                    std::vector<NRLib::Grid2D<double> > &value_above_extra_param,
                                    const NRLib::Grid<double>           &vp_grid,
                                    const NRLib::Grid<double>           &vs_grid,
                                    const NRLib::Grid<double>           &rho_grid,
                                    std::vector<NRLib::Grid<double> >   & parameter_grid_from_eclipse,
                                    size_t i, size_t j, size_t k,
                                    bool top, bool bot, bool right, bool left)
{
  NRLib::StormContGrid &vpgrid  = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid  = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();
  std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();

  std::vector<double> constvp                        = seismic_parameters.modelSettings()->GetConstVp();
  std::vector<double> constvs                        = seismic_parameters.modelSettings()->GetConstVs();
  std::vector<double> constrho                       = seismic_parameters.modelSettings()->GetConstRho();
  std::vector<double> extra_parameter_default_values = seismic_parameters.modelSettings()->GetExtraParameterDefaultValues();

  size_t topk    = seismic_parameters.topK();
  size_t botk    = seismic_parameters.bottomK();
  double zlimit  = seismic_parameters.modelSettings()->GetZeroThicknessLimit();

  double vp_angle   = vpgrid.GetAngle();
  double cosvpangle = cos(vp_angle);
  double sinvpangle = sin(vp_angle);
  double x_min_rot  = vpgrid.GetXMin() * cos(vp_angle) + vpgrid.GetYMin() * sin(vp_angle);
  double y_min_rot  = vpgrid.GetYMin() * cos(vp_angle) - vpgrid.GetXMin() * sin(vp_angle);

  NRLib::Point              mid_edge1, mid_edge2;
  double                    cell_min_x, cell_max_x, cell_min_y, cell_max_y;
  size_t                    start_ii, start_jj, end_ii, end_jj;
  std::vector<double>       x_rot(4), y_rot(4);
  std::vector<bool>         inside(4);
  std::vector<NRLib::Point> pt_vp(4), pt_vs(4), pt_rho(4);
  std::vector<std::vector<NRLib::Point> > pt_extra_param(n_extra_param);
  for (size_t ii = 0; ii < n_extra_param; ++ii)
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
    pt_vp[0] = geometry.FindCellCenterPoint(i,  j,  k);
    pt_vp[1] = geometry.FindCellCenterPoint(ic, jc, k);
    mid_edge1 = 0.5 * (geometry.FindCornerPoint(i,  j,  k, a_corn[0], b_corn[0], c_corn[0]) + geometry.FindCornerPoint(i,  j,  k, a_corn[1], b_corn[1], c_corn[1])) + 0.5 * (geometry.FindCornerPoint(i,  j,  k, a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(i,  j,  k, a_corn[3], b_corn[3], c_corn[3]));
    mid_edge2 = 0.5 * (geometry.FindCornerPoint(ic, jc, k, a_corn[0], b_corn[0], c_corn[0]) + geometry.FindCornerPoint(ic, jc, k, a_corn[1], b_corn[1], c_corn[1])) + 0.5 * (geometry.FindCornerPoint(ic, jc, k, a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(ic, jc, k, a_corn[3], b_corn[3], c_corn[3]));
  }
  else {
    pt_vp[0] = geometry.FindCellCenterPoint(i,  j,  k - 1);
    pt_vp[1] = geometry.FindCellCenterPoint(ic, jc, k - 1);
    mid_edge1 = 0.5 * (geometry.FindCornerPoint(i,  j,  k - 1, a_corn[0], b_corn[0], c_corn[0]) + geometry.FindCornerPoint(i,  j,  k - 1, a_corn[1], b_corn[1], c_corn[1])) + 0.5 * (geometry.FindCornerPoint(i,  j,  k - 1, a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(i,  j,  k - 1, a_corn[3], b_corn[3], c_corn[3]));
    mid_edge2 = 0.5 * (geometry.FindCornerPoint(ic, jc, k - 1, a_corn[0], b_corn[0], c_corn[0]) + geometry.FindCornerPoint(ic, jc, k - 1, a_corn[1], b_corn[1], c_corn[1])) + 0.5 * (geometry.FindCornerPoint(ic, jc, k - 1, a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(ic, jc, k - 1, a_corn[3], b_corn[3], c_corn[3]));
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
      for (size_t ii = 0; ii < n_extra_param; ++ii) {
        pt_extra_param[ii][pt] = pt_vp[pt];
      }
    }
    if (k == botk + 1) {
      for (size_t pt = 0; pt < 2; ++pt) { //nb, only loop two first points here
        pt_vp[pt].z  = constvp[2];
        pt_vs[pt].z  = constvs[2];
        pt_rho[pt].z = constrho[2];
        for (size_t ii = 0; ii < n_extra_param; ++ii) {
          pt_extra_param[ii][pt].z = 0.0;
        }
      }
    }
    else {
      FindPointZValue(i, j, k, pt_vp[0],  geometry, vp_grid, value_above_vp,   constvp[1],  zlimit);
      FindPointZValue(i, j, k, pt_vs[0],  geometry, vs_grid, value_above_vs,   constvs[1],  zlimit);
      FindPointZValue(i, j, k, pt_rho[0], geometry, rho_grid, value_above_rho, constrho[1], zlimit);

      FindPointZValue(ic, jc, k, pt_vp[1],  geometry, vp_grid,  value_above_vp,  constvp[1],  zlimit);
      FindPointZValue(ic, jc, k, pt_vs[1],  geometry, vs_grid,  value_above_vs,  constvs[1],  zlimit);
      FindPointZValue(ic, jc, k, pt_rho[1], geometry, rho_grid, value_above_rho, constrho[1], zlimit);

      for (size_t ii = 0; ii < n_extra_param; ++ii) {
        FindPointZValue(i,  j,  k, pt_extra_param[ii][0], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
        FindPointZValue(ic, jc, k, pt_extra_param[ii][1], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
      }
    }
    for (size_t pt = 2; pt < 4; ++pt) { //nb, only loop two last points here
      pt_vp[pt].z  = pt_vp[pt-2].z;
      pt_vs[pt].z  = pt_vs[pt-2].z;
      pt_rho[pt].z = pt_rho[pt-2].z;
      for (size_t ii = 0; ii < n_extra_param; ++ii) {
        pt_extra_param[ii][pt].z = pt_extra_param[ii][pt-2].z;
      }
    }

    bool triangulate_124 = Is124Triangulate(pt_vp);

    std::vector<NRLib::Triangle> triangles_elastic(6);
    std::vector<NRLib::Triangle> triangles_extra_param(n_extra_param*2);
    SetElasticTriangles(pt_vp, pt_vs, pt_rho, pt_extra_param, triangulate_124, triangles_elastic, triangles_extra_param);

    for (size_t pt = 0; pt < 4; ++pt) {
      x_rot[pt] = pt_vp[pt].x * cosvpangle + pt_vp[pt].y *sinvpangle;
      y_rot[pt] = pt_vp[pt].y * cosvpangle - pt_vp[pt].x *sinvpangle;
    }

    cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
    cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
    cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
    cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

    start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
    start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
    end_ii   = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
    end_jj   = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));
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
      inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j, k,     a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(i, j, k,     a_corn[3], b_corn[3], c_corn[3])));
    } else {
      inside_e_cells.AddPoint(0.5 * (geometry.FindCornerPoint(i, j, k - 1, a_corn[2], b_corn[2], c_corn[2]) + geometry.FindCornerPoint(i, j, k - 1, a_corn[3], b_corn[3], c_corn[3])));
    }

    inside_e_cells.AddPoint(mid_edge1);
    for (size_t ii = start_ii; ii < end_ii; ii++) {
      for (size_t jj = start_jj; jj < end_jj; jj++) {
        double x, y, z;
        vpgrid.FindCenterOfCell(ii, jj, 0, x, y, z);
        NRLib::Point p1(x, y, 0.0);
        NRLib::Point p2(x, y, 1000.0);
        if (inside_e_cells.IsInsidePolygonXY(p1)) {
          NRLib::Line line(p1, p2, false, false);
          NRLib::Point intersec_pt;
          double dist = triangles_elastic[0].FindNearestPoint(line, intersec_pt);
          if (dist < 0.00000000001) {
            vpgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
            triangles_elastic[2].FindIntersection(line, intersec_pt, true);
            vsgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
            triangles_elastic[4].FindIntersection(line, intersec_pt, true);
            rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            for (size_t iii = 0; iii < n_extra_param; ++iii) {
              triangles_extra_param[iii * 2].FindIntersection(line, intersec_pt, true);
              extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
            }
          }
          else {
            dist = triangles_elastic[1].FindNearestPoint(line, intersec_pt);
            if (dist < 0.00000000001) {
              vpgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
              triangles_elastic[3].FindIntersection(line, intersec_pt, true);
              vsgrid(ii, jj, (k - topk) + 1)  = static_cast<float>(intersec_pt.z);
              triangles_elastic[5].FindIntersection(line, intersec_pt, true);
              rhogrid(ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
              for (size_t iii = 0; iii < n_extra_param; ++iii) {
                triangles_extra_param[iii * 2 + 1].FindIntersection(line, intersec_pt, true);
                extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(intersec_pt.z);
              }
            }
          }
        }
      }
    }
  }
}

void SeismicRegridding::FindVpCorners(const NRLib::EclipseGeometry        &geometry,
                                      size_t                              n_extra_param,
                                      SeismicParameters                   &seismic_parameters,
                                      NRLib::Grid2D<double>               &value_above_vp,
                                      NRLib::Grid2D<double>               &value_above_vs,
                                      NRLib::Grid2D<double>               &value_above_rho,
                                      std::vector<NRLib::Grid2D<double> > &value_above_extra_param,
                                      const NRLib::Grid<double>           &vp_grid,
                                      const NRLib::Grid<double>           &vs_grid,
                                      const NRLib::Grid<double>           &rho_grid,
                                      std::vector<NRLib::Grid<double> >   &parameter_grid_from_eclipse,
                                      size_t i, size_t j, size_t k,
                                      std::vector<NRLib::Point>           &pt_vp)
{
  NRLib::StormContGrid &vpgrid  = seismic_parameters.vpGrid();
  NRLib::StormContGrid &vsgrid  = seismic_parameters.vsGrid();
  NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();
  std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();

  std::vector<double> constvp                        = seismic_parameters.modelSettings()->GetConstVp();
  std::vector<double> constvs                        = seismic_parameters.modelSettings()->GetConstVs();
  std::vector<double> constrho                       = seismic_parameters.modelSettings()->GetConstRho();
  std::vector<double> extra_parameter_default_values = seismic_parameters.modelSettings()->GetExtraParameterDefaultValues();

  size_t topk   = seismic_parameters.topK();
  size_t botk   = seismic_parameters.bottomK();
  double zlimit = seismic_parameters.modelSettings()->GetZeroThicknessLimit();

  double vp_angle   = vpgrid.GetAngle();
  double cosvpangle = cos(vp_angle);
  double sinvpangle = sin(vp_angle);
  double x_min_rot  = vpgrid.GetXMin() * cos(vp_angle) + vpgrid.GetYMin() * sin(vp_angle);
  double y_min_rot  = vpgrid.GetYMin() * cos(vp_angle) - vpgrid.GetXMin() * sin(vp_angle);

  double                    cell_min_x, cell_max_x, cell_min_y, cell_max_y;
  size_t                    start_ii, start_jj, end_ii, end_jj;
  std::vector<double>       x_rot(4), y_rot(4);
  std::vector<bool>         inside(4);
  std::vector<NRLib::Point> pt_vs(4), pt_rho(4);
  std::vector<std::vector<NRLib::Point> > pt_extra_param(n_extra_param);
  for (size_t ii = 0; ii < n_extra_param; ++ii)
    pt_extra_param[ii] = pt_vs;

  for (size_t pt = 0; pt < 4;++pt)
    inside[pt] = vpgrid.IsInside(pt_vp[pt].x, pt_vp[pt].y);
  if (inside[0] || inside[1] || inside[2] || inside[3]) {
    pt_vs[0]  = pt_vp[0];
    pt_rho[0] = pt_vp[0];
    for (size_t ii = 0; ii < n_extra_param; ++ii) {
      pt_extra_param[ii][0] = pt_vp[0];
    }
    if (k == botk + 1) {
      pt_vp[3].z  = constvp[2];
      pt_vs[3].z  = constvs[2];
      pt_rho[3].z = constrho[2];
      for (size_t ii = 0; ii < n_extra_param; ++ii) {
        pt_extra_param[ii][3].z = 0.0;
      }
    } 
    else {
      FindPointZValue(i, j, k, pt_vp[3],  geometry, vp_grid,  value_above_vp,  constvp[1],  zlimit);
      FindPointZValue(i, j, k, pt_vs[3],  geometry, vs_grid,  value_above_vs,  constvs[1],  zlimit);
      FindPointZValue(i, j, k, pt_rho[3], geometry, rho_grid, value_above_rho, constrho[1], zlimit);

      for (size_t ii = 0; ii < n_extra_param; ++ii) {
        FindPointZValue(i, j, k, pt_extra_param[ii][3], geometry, parameter_grid_from_eclipse[ii], value_above_extra_param[ii], extra_parameter_default_values[ii], zlimit);
      }
    }

    for (size_t pt = 0; pt < 4; ++pt) {
      x_rot[pt] = pt_vp[pt].x * cosvpangle + pt_vp[pt].y *sinvpangle;
      y_rot[pt] = pt_vp[pt].y * cosvpangle - pt_vp[pt].x *sinvpangle;
    }

    cell_min_x = min(min(x_rot[0], x_rot[1]), min(x_rot[2], x_rot[3]));
    cell_min_y = min(min(y_rot[0], y_rot[1]), min(y_rot[2], y_rot[3]));
    cell_max_x = max(max(x_rot[0], x_rot[1]), max(x_rot[2], x_rot[3]));
    cell_max_y = max(max(y_rot[0], y_rot[1]), max(y_rot[2], y_rot[3]));

    start_ii = static_cast<unsigned int>(max(0.0, (cell_min_x - x_min_rot) / vpgrid.GetDX() - 2.0));
    start_jj = static_cast<unsigned int>(max(0.0, (cell_min_y - y_min_rot) / vpgrid.GetDY() - 2.0));
    end_ii   = static_cast<unsigned int>(max(0.0, (cell_max_x - x_min_rot) / vpgrid.GetDX() + 2.0));
    end_jj   = static_cast<unsigned int>(max(0.0, (cell_max_y - y_min_rot) / vpgrid.GetDY() + 2.0));

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
          for (size_t iii = 0; iii < n_extra_param; ++iii) {
            extra_parameter_grid[iii](ii, jj, (k - topk) + 1) = static_cast<float>(pt_extra_param[ii][3].z);
          }
        }
      }
    }
  }
}

void SeismicRegridding::SetElasticTriangles(std::vector<NRLib::Point>               & pt_vp,
                                            std::vector<NRLib::Point>               & pt_vs,
                                            std::vector<NRLib::Point>               & pt_rho,
                                            std::vector<std::vector<NRLib::Point> > & pt_extra_param,
                                            bool                                      triangulate_124,
                                            std::vector<NRLib::Triangle>            & triangles_elastic,
                                            std::vector<NRLib::Triangle>            & triangles_extra_param)
{
  if (triangulate_124) {
    triangles_elastic[0].SetCornerPoints(pt_vp[0],  pt_vp[1],  pt_vp[3]);
    triangles_elastic[1].SetCornerPoints(pt_vp[0],  pt_vp[2],  pt_vp[3]);
    triangles_elastic[2].SetCornerPoints(pt_vs[0],  pt_vs[1],  pt_vs[3]);
    triangles_elastic[3].SetCornerPoints(pt_vs[0],  pt_vs[2],  pt_vs[3]);
    triangles_elastic[4].SetCornerPoints(pt_rho[0], pt_rho[1], pt_rho[3]);
    triangles_elastic[5].SetCornerPoints(pt_rho[0], pt_rho[2], pt_rho[3]);
    for (size_t ii = 0; ii < triangles_extra_param.size()/2; ++ii){
      triangles_extra_param[ii*2].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][1], pt_extra_param[ii][3]);
      triangles_extra_param[ii*2 + 1].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][2], pt_extra_param[ii][3]);
    }
  }
  else {
    triangles_elastic[0].SetCornerPoints(pt_vp[0],  pt_vp[1],  pt_vp[2]);
    triangles_elastic[1].SetCornerPoints(pt_vp[1],  pt_vp[2],  pt_vp[3]);
    triangles_elastic[2].SetCornerPoints(pt_vs[0],  pt_vs[1],  pt_vs[2]);
    triangles_elastic[3].SetCornerPoints(pt_vs[1],  pt_vs[2],  pt_vs[3]);
    triangles_elastic[4].SetCornerPoints(pt_rho[0], pt_rho[1], pt_rho[2]);
    triangles_elastic[5].SetCornerPoints(pt_rho[1], pt_rho[2], pt_rho[3]);
    for (size_t ii = 0; ii < triangles_extra_param.size()/2; ++ii){
      triangles_extra_param[ii*2].SetCornerPoints(pt_extra_param[ii][0], pt_extra_param[ii][1], pt_extra_param[ii][2]);
      triangles_extra_param[ii*2 + 1].SetCornerPoints(pt_extra_param[ii][1], pt_extra_param[ii][2], pt_extra_param[ii][3]);
    }
  }
}


bool SeismicRegridding::Is124Triangulate(std::vector<NRLib::Point> pt_vp)
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

void SeismicRegridding::GetCornerPointDir(std::vector<size_t> &a,
                                          std::vector<size_t> &b,
                                          std::vector<size_t> &c, 
                                          bool                 left,
                                          bool                 right,
                                          bool                 bot,
                                          bool                 top)
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

bool SeismicRegridding::FindTopCell(const NRLib::EclipseGeometry &geometry,
                                    size_t  i,
                                    size_t &jj)
{
  int j = static_cast<int>(jj);
  while (j >= 0 && !(geometry.IsPillarActive(i, j) && geometry.IsPillarActive(i + 1, j) 
    && geometry.IsPillarActive(i, j + 1) && geometry.IsPillarActive(i + 1, j + 1)
    && geometry.IsPillarActive(i + 2, j) && geometry.IsPillarActive(i + 2, j + 1))) {
      j--;
  }
  if (j >= 0){
    jj = static_cast<size_t>(j);
    return true;
  }
  else
    return false;
}
bool SeismicRegridding::FindBotCell(const NRLib::EclipseGeometry &geometry,
                                    size_t  nj,
                                    size_t  i,
                                    size_t &j)
{
  while (j < nj && !(geometry.IsPillarActive(i, j) && geometry.IsPillarActive(i + 1, j) 
    && geometry.IsPillarActive(i, j + 1) && geometry.IsPillarActive(i + 1, j + 1)
    && geometry.IsPillarActive(i + 2, j) && geometry.IsPillarActive(i + 2, j + 1))) {
      j++;
  }
  if (j < nj)
    return true;
  else
    return false;
}
bool SeismicRegridding::FindLeftCell(const NRLib::EclipseGeometry &geometry,
                                     size_t  ni,
                                     size_t &i,
                                     size_t  j)
{
  while (i < ni && !(geometry.IsPillarActive(i, j) && geometry.IsPillarActive(i, j + 1) 
    && geometry.IsPillarActive(i + 1, j) && geometry.IsPillarActive(i + 1, j + 1)
    && geometry.IsPillarActive(i, j + 2) && geometry.IsPillarActive(i + 1, j + 2))) {
      i++;
  }
  if (i < ni)
    return true;
  else
    return false;
}
bool SeismicRegridding::FindRightCell(const NRLib::EclipseGeometry &geometry,
                                      size_t &ii,
                                      size_t  j)
{
  int i = static_cast<int>(ii);
  while (i >= 0 && !(geometry.IsPillarActive(i, j) && geometry.IsPillarActive(i, j + 1)
    && geometry.IsPillarActive(i + 1, j) && geometry.IsPillarActive(i + 1, j + 1)
    && geometry.IsPillarActive(i, j + 2) && geometry.IsPillarActive(i + 1, j + 2))) {
      i--;
  }
  if (i >= 0){
    ii = static_cast<size_t>(i);
    return true;
  }
  else
    return false;
}

void SeismicRegridding::FindCornerCellPoints(const NRLib::EclipseGeometry &geometry,
                                             std::vector<NRLib::Point>    &pt_vp,
                                             size_t                        i,
                                             size_t                        j,
                                             size_t                        k,
                                             size_t                        botk)
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



