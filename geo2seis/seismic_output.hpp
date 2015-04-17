#ifndef SEISMIC_OUTPUT_HPP
#define SEISMIC_OUTPUT_HPP

#include <stdio.h>
#include <string>
#include <vector>
#include "modelsettings.hpp"
#include "seismic_parameters.hpp"


#include <nrlib/surface/regularsurface.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>

class SeismicParameters;

class SeismicOutput {
  public:
    SeismicOutput(ModelSettings *model_settings);

    void writeDepthSurfaces(const NRLib::RegularSurface<double> &top_eclipse, const NRLib::RegularSurface<double> &bottom_eclipse);

    void writeReflections(SeismicParameters &seismic_parameters, bool noise_added, bool use_window);

    void writeTimeSurfaces(SeismicParameters &seismic_parameters);
    void writeElasticParametersTimeSegy(SeismicParameters &seismic_parameters);
    void writeElasticParametersDepthSegy(SeismicParameters &seismic_parameters);
    void writeExtraParametersTimeSegy(SeismicParameters &seismic_parameters);
    void writeExtraParametersDepthSegy(SeismicParameters &seismic_parameters);

    void writeVpVsRho(SeismicParameters &seismic_parameters, bool use_window);
    void writeZValues(SeismicParameters &seismic_parameters, bool use_window);
    void writeTwt(SeismicParameters &seismic_parameters, bool use_window);

    void writeSeismicTimeSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec);
    void writeSeismicTimeStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec);
    void writeSeismicTimeshiftSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec);
    void writeSeismicTimeshiftStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec);
    void writeSeismicDepthSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec);
    void writeSeismicDepthStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec);

    void writeSeismicStackTime(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec);
    void writeSeismicStackTimeshift(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec);
    void writeSeismicStackDepth(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec);

  private:
    void generateParameterGridForOutput(NRLib::StormContGrid &input_grid, NRLib::StormContGrid &time_or_depth_grid, NRLib::StormContGrid &output_grid, double delta_time_or_depth, double zero_time_or_depth, NRLib::RegularSurface<double> &toptime);
    size_t findCellIndex(size_t i, size_t j, double target_k, NRLib::StormContGrid &grid);

    double top_time_window;
    double bot_time_window;
    bool time_window;
    bool depth_window;
    double top_depth_window;
    double bot_depth_window;

    std::string prefix;
    std::string suffix;

    std::vector<std::string> extra_parameter_names;

    int inline_start;
    int xline_start;
    std::string inline_direction;
    short scalco;
    bool xline_x_axis;
    int inline_step;
    int xline_step;
};

#endif
