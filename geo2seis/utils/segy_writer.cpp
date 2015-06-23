#include "segy_writer.hpp"

#include "nrlib/segy/segy.hpp"

void SEGY::writeSegy(NRLib::StormContGrid      &data,
                     std::string                fileName,
                     int                        inline_start,
                     int                        xline_start,
                     bool                       xline_x_axis,
                     int                        inline_step,
                     int                        xline_step,
                     const NRLib::SegyGeometry *geometry_in,
                     short                      scalco,
                     double                     top_window,
                     double                     bot_window,
                     bool                       window_specified)
{

  size_t i, j;
  int k;
  NRLib::TextualHeader header = NRLib::TextualHeader::standardHeader();
  size_t nx = data.GetNI();
  size_t ny = data.GetNJ();
  double xlstepx, xlstepy, ilstepx, ilstepy;
  double rot = data.GetAngle();
  double dx = data.GetDX();
  double dy = data.GetDY();
  xlstepx = cos(rot) / dx;
  xlstepy = sin(rot) / dx;
  ilstepx = -sin(rot) / dy;
  ilstepy = cos(rot) / dy;
  if (xline_x_axis == false) {
    double temp = xlstepx;
    xlstepx = ilstepx;
    ilstepx = temp;
    temp = xlstepy;
    xlstepy = ilstepy;
    ilstepy = temp;
  }
  ilstepx *= inline_step;
  ilstepy *= inline_step;
  xlstepx *= xline_step;
  xlstepy *= xline_step;

  const NRLib::SegyGeometry *geometry;
  if (geometry_in == NULL)
    geometry = new NRLib::SegyGeometry(data.GetXMin(), data.GetYMin(), dx, dy,
    nx, ny, inline_start - 0.5, xline_start - 0.5, ilstepx, ilstepy, xlstepx, xlstepy, rot);
  else {
    geometry = new NRLib::SegyGeometry(geometry_in);
  }

  ///---Check that the precision asked for is not too high.---
  double x_max_coord = 0;
  double y_max_coord = 0;
  double cos_rot = geometry->GetCosRot();
  double sin_rot = geometry->GetSinRot();
  double x0_coord = geometry->GetX0();
  double y0_coord = geometry->GetY0();
  double xt_coord, yt_coord, x_coord, y_coord;
  for (i = 0; i < nx; i = i + (nx - 1)) {
    for (j = 0; j < ny; j = j + (ny - 1)) {
      xt_coord = double((i + 0.5) * dx);
      yt_coord = double((j + 0.5) * dy);
      x_coord = double(x0_coord + xt_coord * cos_rot - yt_coord * sin_rot);
      y_coord = double(y0_coord + yt_coord * cos_rot + xt_coord * sin_rot);
      if (x_coord > x_max_coord) {
        x_max_coord = x_coord;
      }
      if (y_coord > y_max_coord) {
        y_max_coord = y_coord;
      }
    }
  }
  double utmx_test, utmy_test;
  if (scalco < 0) {
    utmx_test = x_max_coord * scalco * -1;
    utmy_test = y_max_coord * scalco * -1;
  } else {
    utmx_test = x_max_coord / scalco;
    utmy_test = y_max_coord / scalco;
  }
  double max_int_value = 2147483647;
  if (utmx_test > max_int_value || utmy_test > max_int_value) {
    printf("Required precision of UTM coordinates not possible. No Segy file written. Try a lower precision.\n");
    return;
  }
  //-----------------------

  double z_min = data.GetZMin();
  double z_max = data.GetZMax();
  double dz     = float(data.GetLZ() / data.GetNK());
  double z0     = 0.0;
  int    nz     = static_cast<int>(ceil(data.GetZMax()) / dz);
  if (window_specified) {
    if (top_window > z_max) {
      printf("Top window is below grid. No Segy file written.\n");
      return;
    }
    if (bot_window < z_min) {
      printf("Bottom window is above grid. No Segy file written.\n");
      return;
    }


    if (top_window != -9999) {
      z0 = top_window;
    }
    if (bot_window != -9999) {
      nz = static_cast<int>(ceil(bot_window - z0) / dz);
    }
  } else if (nz < 0) {
    printf("Maximum depth is negative. No Segy file written.\n");
    return;
  }

  NRLib::TraceHeaderFormat thf(2);
  NRLib::SegY segyout(fileName, z0, nz, dz, header, thf);
  segyout.SetGeometry(geometry);
  segyout.SetDelayRecTime(z0);
  // geometry->WriteGeometry();
  // geometry->WriteILXL();
  std::vector<float> datavec;
  datavec.resize(nz);
  bool above_zero = false;
  double x, y, xt, yt, z;
  for (size_t j = 0; j < ny; j++) {
    for (size_t i = 0; i < nx; i++) {
      geometry->FindXYFromIJ(i, j, x, y);

      double zbot = data.GetBotSurface().GetZ(x, y);
      double ztop = data.GetTopSurface().GetZ(x, y);
      int firstData = static_cast<int>(floor((ztop) / dz));
      int endData = static_cast<int>(floor((zbot) / dz));
      int windowTop = static_cast<int>(floor((top_window) / dz));
      int windowBot = static_cast<int>(floor((bot_window) / dz));

      if (window_specified) {
        if (windowTop < firstData) {
          for (k = windowTop; k < firstData; k++) {
            datavec[k - windowTop] = 0.0;
          }
        } else {
          firstData = windowTop;
        }
        if (windowBot < endData) {
          endData = windowBot;
        }
        for (k = firstData; k < endData; k++) {
          z = k * dz;
          datavec[k - windowTop] = float(data.GetValueZInterpolated(x, y, z));
        }
        if (windowBot > endData) {
          for (k = endData; k < windowBot; k++) {
            datavec[k - windowTop] = 0.0;
          }
        }
      } else {
        if (firstData < 0) {
          //if(above_zero == false)
          //printf("Internal warning: SEGY-grid starts at 0. Truncating data at top.\n");
          above_zero = true;
          firstData = 0;
        }
        if (endData > nz) {
          //printf("Internal warning: SEGY-grid too small (%d, %d needed). Truncating data.\n", static_cast<int>(nz), static_cast<int>(endData));
          endData = nz;
        }
        for (k = 0; k < firstData; k++) {
          datavec[k] = 0.0;
        }
        for (k = firstData; k < endData; k++) {
          z = z0 + k * dz;
          datavec[k] = float(data.GetValueZInterpolated(x, y, z));
        }
        for (k = endData; k < nz; k++) {
          datavec[k] = 0.0;
        }
        //segyout.StoreTrace(x,y,datavec,NULL);
      }
      segyout.StoreTrace(x,y, datavec, NULL);
      //segyout.WriteTrace(x, y, datavec, NULL, 0.0, 0.0, scalco);
    }
  }
  segyout.WriteAllTracesToFile(scalco);
  delete geometry;
}