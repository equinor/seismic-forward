#ifndef STORM_WRITER_HPP
#define STORM_WRITER_HPP

#include "nrlib/stormgrid/stormcontgrid.hpp"
#include <string>

class STORM {
  public:
    static void WriteStorm(NRLib::StormContGrid &grid,
            std::string &filename,
            double top_window,
            double bot_window,
            bool window_specified,
            bool keep_grid = false);

};

#endif
