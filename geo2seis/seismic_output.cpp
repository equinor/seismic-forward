#include <seismic_output.hpp>
#include <utils/storm_writer.hpp>
#include "seismic_parameters.hpp"
#include "seismic_geometry.hpp"
#include "nrlib/geometry/interpolation.hpp"

#ifdef WITH_OMP
#include <omp.h>
#endif

SeismicOutput::SeismicOutput(ModelSettings *model_settings) {
  top_time_window_  = model_settings->GetTopTimeWindow();
  bot_time_window_  = model_settings->GetBotTimeWindow();
  time_window_      = model_settings->GetTimeWindowSpecified();
  depth_window_     = model_settings->GetDepthWindowSpecified();
  top_depth_window_ = model_settings->GetTopDepthWindow();
  bot_depth_window_ = model_settings->GetBotDepthWindow();

  prefix_ = "";
  if (model_settings->GetPrefix() != "") {
    prefix_ = model_settings->GetPrefix() + "_";
  }

  suffix_ = "";
  if (model_settings->GetSuffix() != "") {
    suffix_ = "_" + model_settings->GetSuffix();
  }

  extra_parameter_names_ = model_settings->GetExtraParameterNames();

  //-----------------Get segy indexes for print----------------------------
  inline_start_     = model_settings->GetSegyInlineStart();
  xline_start_      = model_settings->GetSegyXlineStart();
  inline_direction_ = model_settings->GetSegyInlineDirection();

  //-----------------UTM precision in segy header--------------------------
  scalco_ = model_settings->GetUtmPrecision();

  xline_x_axis_ = true;
  if (NRLib::Uppercase(inline_direction_) == "X") {
    xline_x_axis_ = false;
  } else if (NRLib::Uppercase(inline_direction_) == "Y") {
    xline_x_axis_ = true;
  }
  inline_step_ = model_settings->GetSegyInlineStep();
  xline_step_  = model_settings->GetSegyXlineStep();
}

void SeismicOutput::SetSegyGeometry(SeismicParameters   &seismic_parameters,
                                    const NRLib::Volume &vol,
                                    size_t               nx,
                                    size_t               ny)
{

  double xlstepx, xlstepy, ilstepx, ilstepy;
  double rot = vol.GetAngle();
  double dx = vol.GetLX()/nx; //nb check
  double dy = vol.GetLY()/ny; //nb check
  xlstepx = cos(rot) / dx;
  xlstepy = sin(rot) / dx;
  ilstepx = -sin(rot) / dy;
  ilstepy = cos(rot) / dy;
  if (xline_x_axis_ == false) {
    double temp = xlstepx;
    xlstepx = ilstepx;
    ilstepx = temp;
    temp = xlstepy;
    xlstepy = ilstepy;
    ilstepy = temp;
  }
  ilstepx *= inline_step_;
  ilstepy *= inline_step_;
  xlstepx *= xline_step_;
  xlstepy *= xline_step_;

  if (seismic_parameters.GetSegyGeometry() == NULL){
    const NRLib::SegyGeometry *geometry =
      new NRLib::SegyGeometry(vol.GetXMin(), vol.GetYMin(), dx, dy,
                              nx, ny, inline_start_ - 0.5, xline_start_ - 0.5,
                              ilstepx, ilstepy, xlstepx, xlstepy, rot);
    seismic_parameters.SetSegyGeometry(geometry);
    delete geometry;
  }
}

bool SeismicOutput::CheckUTMPrecision(SeismicParameters   &seismic_parameters,
                                      const NRLib::Volume &vol,
                                      size_t               nx,
                                      size_t               ny)
{
  ///---Check that the precision asked for is not too high.---
  NRLib::SegyGeometry *geometry   = seismic_parameters.GetSegyGeometry();
  double dx          = vol.GetLX()/nx; //nb check
  double dy          = vol.GetLY()/ny; //nb check
  double x_max_coord = 0;
  double y_max_coord = 0;
  double cos_rot     = geometry->GetCosRot();
  double sin_rot     = geometry->GetSinRot();
  double x0_coord    = geometry->GetX0();
  double y0_coord    = geometry->GetY0();
  double xt_coord, yt_coord, x_coord, y_coord;
  for (size_t i = 0; i < nx; i = i + (nx - 1)) {
    for (size_t j = 0; j < ny; j = j + (ny - 1)) {
      xt_coord = double((i + 0.5) * dx);
      yt_coord = double((j + 0.5) * dy);
      x_coord  = double(x0_coord + xt_coord * cos_rot - yt_coord * sin_rot);
      y_coord  = double(y0_coord + yt_coord * cos_rot + xt_coord * sin_rot);
      if (x_coord > x_max_coord) {
        x_max_coord = x_coord;
      }
      if (y_coord > y_max_coord) {
        y_max_coord = y_coord;
      }
    }
  }
  double utmx_test, utmy_test;
  if (scalco_ < 0) {
    utmx_test = x_max_coord * scalco_ * -1;
    utmy_test = y_max_coord * scalco_ * -1;
  } else {
    utmx_test = x_max_coord / scalco_;
    utmy_test = y_max_coord / scalco_;
  }
  double max_int_value = 2147483647;
  if (utmx_test > max_int_value || utmy_test > max_int_value) {
    printf("Required precision of UTM coordinates not possible. No Segy file written. Try a lower precision.\n");
    return false;
  }
  return true;
}

bool SeismicOutput::PrepareSegy(NRLib::SegY               &segyout,
                                const std::vector<double> &twt_0,
                                size_t                     n_samples,
                                std::string                fileName,
                                SeismicParameters         &seismic_parameters,
                                const std::vector<double> &offset_vec,
                                size_t                     n_traces_per_ensamble,
                                bool                       time,
                                bool                       nmo)
{
  double z_min = twt_0[0];
  double z_max = twt_0[n_samples-1];
  double dz    = twt_0[1] - twt_0[0];
  double z0     = 0.0; //vurdere å bruke z_min selv uten vindu??
  int    nz     = static_cast<int>(ceil(z_max/dz));
  if ((time == true && time_window_) || (time == false && depth_window_)) {
    if (time == true ) {
      if (top_time_window_ > z_max) {
        printf("Top window is below grid. No Segy file written.\n");
        return false;
      }
      if (bot_time_window_ < z_min) {
        printf("Bottom window is above grid. No Segy file written.\n");
        return false;
      }
      if (top_time_window_ != -9999) {
        z0 = top_time_window_;
      }
      if (bot_time_window_ != -9999) {
        nz = static_cast<int>(ceil(bot_time_window_ - z0) / dz);
        z_max = bot_time_window_;
      }
    }
    else {
      if (top_depth_window_ > z_max) {
        printf("Top window is below grid. No Segy file written.\n");
        return false;
      }
      if (bot_depth_window_ < z_min) {
        printf("Bottom window is above grid. No Segy file written.\n");
        return false;
      }
      if (top_depth_window_ != -9999) {
        z0 = top_depth_window_;
      }
      if (bot_depth_window_ != -9999) {
        nz = static_cast<int>(ceil(bot_depth_window_ - z0) / dz);
        z_max = bot_depth_window_;
      }
    }
  }
  else if (nz < 0) {
    printf("Maximum depth is negative. No Segy file written.\n");
    return false;
  }

  NRLib::SegyGeometry *geometry     = seismic_parameters.GetSegyGeometry();
  geometry->FindILXLGeometry();
  std::string filename_out = prefix_ + fileName + suffix_ + ".segy";
  NRLib::TraceHeaderFormat thf(2);
  //Textual header:
  NRLib::TextualHeader header = NRLib::TextualHeader();
  header.SetLine(0, "SEGY OUTPUT FROM Seismic Forward Modeling / Geo2Seis, ver 4.0  2015");
  std::string line = "Name: " + filename_out, line2;
  header.SetLine(1, line);
  line = "  First inline: " + NRLib::ToString(geometry->GetMinIL()) + "      Last inline: " + NRLib::ToString(geometry->GetMaxIL());
  header.SetLine(2, "CDP:");
  header.SetLine(3, line);
  line = "  First xline:  " + NRLib::ToString(geometry->GetMinXL()) + "      Last xline: " + NRLib::ToString(geometry->GetMaxXL());
  header.SetLine(4, line);
  if (nmo) {
    line = "Offset (m)    min: " + NRLib::ToString(offset_vec[0]) + "    max: " + NRLib::ToString(offset_vec[offset_vec.size() - 1]);
  }
  else {
    line = "Angle  (m)    min: " + NRLib::ToString(offset_vec[0] / NRLib::Degree) + "    max: " + NRLib::ToString(offset_vec[offset_vec.size() - 1] / NRLib::Degree);
  }
  header.SetLine(5, line);
  if (time) {
    line = "Time (ms)     min: " + NRLib::ToString(z0) + "     max: " + NRLib::ToString(z_max);
    line2 = "Time increment (ms): " + NRLib::ToString(dz);
  }
  else {
    line = "Depth (m)     min: " + NRLib::ToString(z0) + "     max: " + NRLib::ToString(z_max);
    line2 = "Depth increment (m): " + NRLib::ToString(dz);
  }
  header.SetLine(6, line);
  header.SetLine(7, line2);

  segyout.Initialize(filename_out, float(z0), static_cast<size_t>(nz), float(dz), header, thf, short(n_traces_per_ensamble));
  segyout.SetGeometry(geometry);
  segyout.SetDelayRecTime(short(z0));
  return true;
}

void SeismicOutput::WriteSegyGather(NRLib::Grid2D<double>     &data_gather,
                                    NRLib::SegY               &segyout,
                                    const std::vector<double>  twt_0,
                                    const std::vector<double>  offset_vec,
                                    bool                       time,
                                    double                     x,
                                    double                     y,
                                    bool                       nmo)
{
  float z0  = segyout.GetZ0();
  float dz  = segyout.GetDz();
  size_t nz = segyout.GetNz();

  std::vector<float> datavec(nz);

  int k;
  int firstSample = static_cast<int>(floor((twt_0[0])              / dz));
  int endData   = static_cast<int>(floor((twt_0[data_gather.GetNI()-1]) / dz));
  int firstData   = firstSample;

  int windowTop, windowBot;
  if ((time == true && time_window_) || (time == false && depth_window_)) {
    if (time == true ) {
      windowTop = static_cast<int>(floor((top_time_window_) / dz));
      windowBot = windowTop + nz;
    }
    else {
      windowTop = static_cast<int>(floor((top_depth_window_) / dz));
      windowBot = windowTop + nz;
    }
  }
  for (size_t off = 0; off < offset_vec.size(); ++off) {
    if ((time == true && time_window_) || (time == false && depth_window_)) {
      if (windowTop < firstSample) {
        for (k = windowTop; k < firstSample; k++) {
          datavec[k - windowTop] = 0.0;
        }
      }
      else {
        firstData = windowTop;
      }
      if (windowBot < endData) {
        endData = windowBot;
      }
      for (k = firstData; k < endData; k++) {
       datavec[k - windowTop] = float(data_gather(k - firstSample, off));
      }
      if (windowBot > endData) {
        for (k = endData; k < windowBot; k++) {
          datavec[k - windowTop] = 0.0;
        }
      }
    }
    else {
      if (firstData < 0) {
        firstData = 0;
      }
      if (endData > nz) {
        endData = int(nz);
      }
      for (k = 0; k < firstData; k++) {
        datavec[k] = 0.0;
      }
      for (k = firstData; k < endData; k++) {
        datavec[k] = float(data_gather(k - firstSample, off));
      }
      for (k = endData; k < nz; k++) {
        datavec[k] = 0.0;
      }
    }
    if (nmo)
      segyout.WriteTrace(x,y, datavec, NULL, 0.0, 0.0, scalco_, short(offset_vec[off]));
    else
      segyout.WriteTrace(x,y, datavec, NULL, 0.0, 0.0, scalco_, short(std::floor(offset_vec[off]/NRLib::Degree + 0.5)));
  }
}

void SeismicOutput::WriteZeroSegyGather(NRLib::SegY               &segyout,
                                        const std::vector<double>  offset_vec,
                                        double                     x,
                                        double                     y,
                                        bool                       nmo)
{
  size_t nz = segyout.GetNz();

  std::vector<float> datavec(nz);
  for (size_t k = 0; k < nz; ++k) {
    datavec[k] = 0.0;
  }
  if (nmo) {
    for (size_t off = 0; off < offset_vec.size(); ++off) {
      segyout.WriteTrace(x,y, datavec, NULL, 0.0, 0.0, scalco_, short(offset_vec[off]));
    }
  }
  else {
    for (size_t off = 0; off < offset_vec.size(); ++off) {
      segyout.WriteTrace(x,y, datavec, NULL, 0.0, 0.0, scalco_, short(std::floor(offset_vec[off]/NRLib::Degree + 0.5)));
    }
  }
}

void SeismicOutput::WriteDepthSurfaces(const NRLib::RegularSurface<double> &top_eclipse, const NRLib::RegularSurface<double> &bottom_eclipse)
{
  printf("Write depth surfaces on Storm format\n");
  std::string filename = prefix_ + "topeclipse" + suffix_ + ".storm";
  top_eclipse.WriteToFile(filename);
  filename = prefix_ + "boteclipse" + suffix_ + ".storm";
  bottom_eclipse.WriteToFile(filename);
}

void SeismicOutput::WriteTimeSurfaces(SeismicParameters &seismic_parameters)
{
  NRLib::RegularSurface<double> &toptime = seismic_parameters.GetTopTime();
  NRLib::RegularSurface<double> &bottime = seismic_parameters.GetBottomTime();

  printf("Write time surfaces on Storm format.\n");
  bottime.WriteToFile(prefix_ + "bottime" + suffix_ + ".storm");
  toptime.WriteToFile(prefix_ + "toptime" + suffix_ + ".storm");
}

void SeismicOutput::WriteReflections(SeismicParameters &seismic_parameters, double angle_or_offset)
{
  std::vector<NRLib::StormContGrid> &rgridvec = seismic_parameters.GetRGrids();

  printf("Write reflections on Storm format.\n");
  std::string reflection_string = "reflections_";
  std::string filename = prefix_ + reflection_string + NRLib::ToString(angle_or_offset) + suffix_ + ".storm";
  rgridvec[0].WriteToFile(filename);
  if (rgridvec.size() == 2)
  {
    //printf("Write reflections with noise on Storm format.\n");
    reflection_string = reflection_string + "noise_";
    std::string filename = prefix_ + reflection_string + NRLib::ToString(angle_or_offset) + suffix_ + ".storm";
    rgridvec[1].WriteToFile(filename);
  }
}

void SeismicOutput::WriteVrms(SeismicParameters    &seismic_parameters,
                              std::string           name_pp_or_ps)
{
  NRLib::StormContGrid &vrmsgrid = seismic_parameters.GetVrmsGrid();
  std::string message;
  if (name_pp_or_ps != "") {
    message = "Write vrms grid for " + name_pp_or_ps + " on Storm format.\n";
    name_pp_or_ps = "_" + name_pp_or_ps;
  }
  else {
    message = "Write vrms grid on Storm format.\n";
  }
  std::cout << message;
  std::string filename = prefix_ + "vrms" + name_pp_or_ps + suffix_ + ".storm";
  vrmsgrid.WriteToFile(filename);
}

void SeismicOutput::WriteVpVsRho(SeismicParameters &seismic_parameters)
{
  NRLib::StormContGrid &vpgrid  = seismic_parameters.GetVpGrid();
  NRLib::StormContGrid &vsgrid  = seismic_parameters.GetVsGrid();
  NRLib::StormContGrid &rhogrid = seismic_parameters.GetRhoGrid();

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nWrite elastic parameters using STORM format.");
  std::string filename = prefix_ + "vp" + suffix_ + ".storm";
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Vp : %s",filename.c_str());
  vpgrid.WriteToFile(filename);
  filename = prefix_ + "vs" + suffix_ + ".storm";
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Vs : %s",filename.c_str());
  vsgrid.WriteToFile(filename);
  filename = prefix_ + "rho" + suffix_ + ".storm";
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Rho: %s",filename.c_str());
  rhogrid.WriteToFile(filename);
}

void SeismicOutput::WriteZValues(SeismicParameters &seismic_parameters)
{
  NRLib::StormContGrid &zgrid = seismic_parameters.GetZGrid();
  std::string filename = prefix_ + "zgrid" + suffix_ + ".storm";

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nWrite z-values on STORM format to file : %s",filename.c_str());
  zgrid.WriteToFile(filename);
}

void SeismicOutput::WriteTwt(SeismicParameters &seismic_parameters)
{
  NRLib::StormContGrid &twtgrid = seismic_parameters.GetTwtGrid();
  std::string filename = prefix_ + "twt" + suffix_ + ".storm";

  printf("Write two way time on Storm format.\n");
  twtgrid.WriteToFile(filename);
}

void SeismicOutput::WriteSeismicTimeStorm(SeismicParameters &seismic_parameters, NRLib::StormContGrid &timegrid, double offset, bool is_stack)
{
  printf("Write seismic in time on Storm format.\n");
  ModelSettings *model_settings = seismic_parameters.GetModelSettings();
  std::string filename;
  if (is_stack == false) {
    filename = prefix_ + "seismic_time_" + NRLib::ToString(offset) + suffix_ + ".storm";
  }
  else {
    filename = prefix_ + "seismic_time_stack" + suffix_ + ".storm";
  }
  STORM::WriteStorm(timegrid, filename, top_time_window_, bot_time_window_, time_window_);
}

void SeismicOutput::WriteSeismicDepthStorm(SeismicParameters &seismic_parameters, NRLib::StormContGrid &depthgrid, double offset, bool is_stack)
{
  printf("Write seismic in depth on Storm format.\n");
  ModelSettings *model_settings = seismic_parameters.GetModelSettings();
  std::string filename;
  if (is_stack == false) {
    filename = prefix_ + "seismic_depth_" + NRLib::ToString(offset) + suffix_ + ".storm";
  }
  else {
    filename = prefix_ + "seismic_depth_stack" + suffix_ + ".storm";
  }
  STORM::WriteStorm(depthgrid, filename, top_depth_window_, bot_depth_window_, depth_window_);
}

void SeismicOutput::WriteSeismicTimeshiftStorm(SeismicParameters &seismic_parameters, NRLib::StormContGrid &timeshiftgrid, double offset, bool is_stack)
{
  printf("Write seismic in timeshift on Storm format.\n");
  ModelSettings *model_settings = seismic_parameters.GetModelSettings();
  std::string filename;
  if (is_stack == false) {
    filename = prefix_ + "seismic_timeshift_" + NRLib::ToString(offset) + suffix_ + ".storm";
  }
  else {
    filename = prefix_ + "seismic_timeshift_stack" + suffix_ + ".storm";
  }
  STORM::WriteStorm(timeshiftgrid, filename, top_time_window_, bot_time_window_, time_window_);
}

void SeismicOutput::PrintVector(std::vector<double> vec, std::string filename)
{
  std::cout << "debug print: " << filename << "\n";
  std::ofstream fout;
  NRLib::OpenWrite(fout, filename);
  for (size_t i = 0; i < vec.size(); ++i) {
    fout << vec[i] << std::endl;
  }
  fout.close();
}

void SeismicOutput::PrintVectorSizeT(std::vector<size_t> vec, std::string filename)
{
  std::cout << "debug print: " << filename << "\n";
  std::ofstream fout;
  NRLib::OpenWrite(fout, filename);
  for (size_t i = 0; i < vec.size(); ++i) {
    fout << vec[i] << std::endl;
  }
  fout.close();
}

void SeismicOutput::PrintMatrix(NRLib::Grid2D<double> matrix, std::string filename)
{
  std::cout << "debug print: " << filename << "\n";
  std::ofstream fout;
  NRLib::OpenWrite(fout, filename);
  for (size_t i = 0; i < matrix.GetNI(); ++i) {
    for (size_t j = 0; j < matrix.GetNJ(); ++j) {
      fout << matrix(i,j) << "  ";
    }
    fout << std::endl;
  }
  fout.close();
}
