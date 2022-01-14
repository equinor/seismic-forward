#ifndef SEISMIC_REGRIDDING_HPP
#define SEISMIC_REGRIDDING_HPP

#include <seismic_parameters.hpp>
#include <nrlib/grid/grid.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include <nrlib/eclipsegrid/eclipsegeometry.hpp>
#include <nrlib/geometry/point.hpp>
#include <nrlib/geometry/triangle.hpp>
#include "utils/resampl_trace.hpp"
#include "utils/resampl_output.hpp"
#include "utils/gen_resampl_param.hpp"

class SeismicRegridding {
  public:
    static void MakeSeismicRegridding(SeismicParameters &seismic_parameters);

  private:
    static void WriteElasticParametersSegy(SeismicParameters &seismic_parameters,
                                           size_t             n_threads,
                                           bool               time);

    static void WriteExtraParametersSegy(SeismicParameters &seismic_parameters,
                                         size_t n_threads,
                                         bool time);

    static void WriteParametersTimeSegy(SeismicParameters                  &seismic_parameters,
                                        size_t                              n_threads,
                                        std::vector<NRLib::StormContGrid*>  input_grid,
                                        std::vector<std::string>            filenames);

    static void WriteParametersDepthSegy(SeismicParameters                  &seismic_parameters,
                                         size_t                              n_threads,
                                         std::vector<NRLib::StormContGrid*>  input_grid,
                                         std::vector<std::string>            filenames);

    static void WriteParametersSegyInParallel(SeismicParameters                 &seismic_parameters,
                                              size_t                             n_threads,
                                              std::vector<NRLib::StormContGrid*> input_grid,
                                              std::vector<std::string>           filenames,
                                              std::vector<double>               &time_or_depth_vec_reg,
                                              NRLib::StormContGrid              &time_or_depth_grid,
                                              bool                               time);

    static void GenerateParameterGridForOutput(GenResamplParam *params,
                                               Trace           *trace,
                                               ResamplOutput   *resampl_output);

    static void GenerateParameterGridForOutputQueue(GenResamplParam *params,
                                                    ResamplOutput   *resampl_output);

    static void WriteResampledParameter(GenResamplParam *params,
                                        ResamplOutput   *resampl_output);

    static size_t FindCellIndex(size_t                i,
                                size_t                j,
                                double                target_k,
                                NRLib::StormContGrid &grid);

    static void FindZValues(SeismicParameters &seismic_parameters,
                            size_t             n_threads);

    static void FindVp(SeismicParameters &seismic_parameters,
                       size_t             n_threads);

    static void VpPostProcess(SeismicParameters &seismic_parameters);

    static void FindVrms(SeismicParameters          &seismic_parameters,
                         const NRLib::StormContGrid &vgrid,
                         const NRLib::StormContGrid &twtgrid);

    static void FindTWT(SeismicParameters             &seismic_parameters,
                        NRLib::RegularSurface<double> &toptime,
                        NRLib::RegularSurface<double> &bottime,
                        size_t                         n_threads);


    static void FillInGridValues(const NRLib::EclipseGeometry &geometry,
                                 NRLib::Grid<double>          &grid_copy,
                                 double                        default_value,
                                 double                        zlimit,
                                 double                        default_top,
                                 size_t                        ni,
                                 size_t                        nj,
                                 size_t                        topk,
                                 size_t                        botk);

    static void FindVpEdges(const NRLib::EclipseGeometry        &geometry,
                            size_t                               n_extra_param,
                            SeismicParameters                   &seismic_parameters,
                            const NRLib::Grid<double>           &vp_grid,
                            const NRLib::Grid<double>           &vs_grid,
                            const NRLib::Grid<double>           &rho_grid,
                            std::vector<NRLib::Grid<double> >   &parameter_grid_from_eclipse,
                            size_t i, size_t j, size_t k,
                            bool top, bool bot, bool right, bool left);

    static void FindVpCorners(const NRLib::EclipseGeometry        &geometry,
                              size_t                              n_extra_param,
                              SeismicParameters                   &seismic_parameters,
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
