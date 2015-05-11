#ifndef SEISMIC_REGRIDDING_HPP
#define SEISMIC_REGRIDDING_HPP

#include <seismic_parameters.hpp>
#include <nrlib/grid/grid.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/eclipsegrid/eclipsegeometry.hpp>
#include <nrlib/geometry/point.hpp>


class SeismicRegridding {
  public:
    static void seismicRegridding(SeismicParameters &seismic_parameters);
    static void addNoiseToReflections(unsigned long seed, double std_dev, std::vector<NRLib::StormContGrid> &grid_vec);
    static void addNoiseToReflectionsPos(unsigned long seed, double std_dev, std::vector<std::vector<double> > &refl);

  private:
    static void findZValues(SeismicParameters &seismic_parameters);
    static void findVpAndR(SeismicParameters &seismic_parameters);
   
    static void findVrms(SeismicParameters &seismic_parameters);
    
    static void findTWT(NRLib::StormContGrid &vpgrid, NRLib::StormContGrid &vsgrid, NRLib::StormContGrid &twtgrid, NRLib::StormContGrid &zgrid, NRLib::RegularSurface<double> &toptime, NRLib::RegularSurface<double> &bottime, bool ps_seismic);
//    static void findTWT(NRLib::StormContGrid &vpgrid, NRLib::StormContGrid &twtgrid, NRLib::StormContGrid &zgrid,NRLib::RegularSurface<double> &toptime, NRLib::RegularSurface<double> &bottime);

    static void findPointZValue(size_t i, size_t j, size_t k,
                                NRLib::Point &point,
                                const NRLib::EclipseGeometry &geometry,
                                const NRLib::Grid<double> &grid,
                                const NRLib::Grid2D<double> &value_above,
                                double &default_value,
                                double &zlimit);
};

#endif
