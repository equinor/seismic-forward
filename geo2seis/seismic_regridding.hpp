#ifndef SEISMIC_REGRIDDING_HPP
#define SEISMIC_REGRIDDING_HPP

#include <seismic_parameters.hpp>
#include <nrlib/grid/grid.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/eclipsegrid/eclipsegeometry.hpp>
#include <nrlib/geometry/point.hpp>
#include <nrlib/geometry/triangle.hpp>


class SeismicRegridding {
  public:
    static void MakeSeismicRegridding(SeismicParameters &seismic_parameters);
    static void AddNoiseToReflectionsPos(unsigned long         seed,
                                         double                std_dev,
                                         NRLib::Grid2D<double> &refl);

  private:
    static void FindZValues(SeismicParameters &seismic_parameters);

    static void FindVp(SeismicParameters &seismic_parameters);

    static void FindVrms(SeismicParameters          &seismic_parameters,
                         const NRLib::StormContGrid &vgrid,
                         const NRLib::StormContGrid &twtgrid);

    static void FindTWT(SeismicParameters &seismic_parameters,
                        NRLib::RegularSurface<double> &toptime,
                        NRLib::RegularSurface<double> &bottime);

    static void FindPointZValue(size_t i, size_t j, size_t k,
                                NRLib::Point &point,
                                const NRLib::EclipseGeometry &geometry,
                                const NRLib::Grid<double> &grid,
                                const NRLib::Grid2D<double> &value_above,
                                double &default_value,
                                double &zlimit);

    static void FindVpEdges(const NRLib::EclipseGeometry        &geometry,
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
                            bool top, bool bot, bool right, bool left);

    static void FindVpCorners(const NRLib::EclipseGeometry        &geometry,
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
                              std::vector<NRLib::Point>           &pt_vp);

    static void SetElasticTriangles(std::vector<NRLib::Point>               & pt_vp,
                                    std::vector<NRLib::Point>               & pt_vs,
                                    std::vector<NRLib::Point>               & pt_rho,
                                    std::vector<std::vector<NRLib::Point> > & pt_extra_param,
                                    bool                                      triangulate_124,
                                    std::vector<NRLib::Triangle>            & triangles_elastic,
                                    std::vector<NRLib::Triangle>            & triangles_extra_param);

    static bool Is124Triangulate(std::vector<NRLib::Point> pt_vp);

    static void GetCornerPointDir(std::vector<size_t> & a,
                                  std::vector<size_t> & b,
                                  std::vector<size_t> & c,
                                  bool left,
                                  bool right,
                                  bool bot,
                                  bool top);

    static bool FindTopCell(const NRLib::EclipseGeometry &geometry,
                            size_t  i,
                            size_t &j);
    static bool FindBotCell(const NRLib::EclipseGeometry &geometry,
                            size_t  nj,
                            size_t  i,
                            size_t &j);
    static bool FindLeftCell(const NRLib::EclipseGeometry &geometry,
                             size_t  ni,
                             size_t &i,
                             size_t  j);
    static bool FindRightCell(const NRLib::EclipseGeometry &geometry,
                             size_t &i,
                             size_t  j);
    static void FindCornerCellPoints(const NRLib::EclipseGeometry &geometry,
                                     std::vector<NRLib::Point>    &vp_point,
                                     size_t                        i,
                                     size_t                        j,
                                     size_t                        k,
                                     size_t                        botk);
};

#endif
