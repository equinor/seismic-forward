#ifndef SEGY_WRITER_HPP
#define SEGY_WRITER_HPP

#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/segy/segygeometry.hpp"
#include <string>

class SEGY {
  public:
    static void writeSegy(NRLib::StormContGrid &data,
                          std::string fileName,
                          int inline_start,
                          int xline_start,
                          bool xline_x_axis,
                          int inline_step,
                          int xline_step,
                          const NRLib::SegyGeometry *geometry_in,
                          short scalco,
                          double top_window,
                          double bot_window,
                          bool window_specified);

};

#endif
