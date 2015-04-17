#include <seismic_output.hpp>

#include <utils/storm_writer.hpp>
#include "utils/segy_writer.hpp"
#include "seismic_parameters.hpp"
#include "seismic_geometry.hpp"

SeismicOutput::SeismicOutput(ModelSettings *model_settings) {
    top_time_window = model_settings->GetTopTimeWindow();
    bot_time_window = model_settings->GetBotTimeWindow();
    time_window = model_settings->GetTimeWindowSpecified();
    depth_window = model_settings->GetDepthWindowSpecified();
    top_depth_window = model_settings->GetTopDepthWindow();
    bot_depth_window = model_settings->GetBotDepthWindow();

    prefix = "";
    if (model_settings->GetPrefix() != "") {
        prefix = model_settings->GetPrefix() + "_";
    }

    suffix = "";
    if (model_settings->GetSuffix() != "") {
        suffix = "_" + model_settings->GetSuffix();
    }

    extra_parameter_names = model_settings->GetExtraParameterNames();

    //-----------------Get segy indexes for print----------------------------
    inline_start = model_settings->GetSegyInlineStart();
    xline_start = model_settings->GetSegyXlineStart();
    inline_direction = model_settings->GetSegyInlineDirection();

    //-----------------UTM precision in segy header--------------------------
    scalco = model_settings->GetUtmPrecision();


    xline_x_axis = true;
    if (NRLib::Uppercase(inline_direction) == "X") {
        xline_x_axis = false;
    } else if (NRLib::Uppercase(inline_direction) == "Y") {
        xline_x_axis = true;
    }
    inline_step = model_settings->GetSegyInlineStep();
    xline_step = model_settings->GetSegyXlineStep();
}

void SeismicOutput::writeDepthSurfaces(const NRLib::RegularSurface<double> &top_eclipse, const NRLib::RegularSurface<double> &bottom_eclipse) {
    printf("Write depth surfaces on Storm format\n");
    std::string filename = prefix + "topeclipse" + suffix + ".storm";
    top_eclipse.WriteToFile(filename);
    filename = prefix + "boteclipse" + suffix + ".storm";
    bottom_eclipse.WriteToFile(filename);
}

void SeismicOutput::writeReflections(SeismicParameters &seismic_parameters, bool noise_added, bool use_window) {
    std::vector<NRLib::StormContGrid> &rgridvec = seismic_parameters.rGrids();
    double theta_0 = seismic_parameters.theta0();
    double dtheta = seismic_parameters.dTheta();
    size_t ntheta = seismic_parameters.nTheta();

    std::string reflection_string = "reflections_";
    if (noise_added) {
        printf("Write reflections with noise on Storm format");
        reflection_string = reflection_string + "noise_";
    } else {
        printf("Write reflections on Storm format");
    }

    if(use_window) {
        printf(" using window\n");

    } else {
        printf(" not using window\n");
    }

    double theta = theta_0;
    for (size_t i = 0; i < ntheta; i++) {
        std::string filename = prefix + reflection_string + NRLib::ToString(NRLib::Radian * theta) + suffix + ".storm";
        if(use_window) {
            STORM::writeStorm((rgridvec)[i], filename, top_depth_window, bot_depth_window, depth_window, true);
        } else {
            rgridvec[i].WriteToFile(filename);
        }

        theta += dtheta;
    }
}

void SeismicOutput::writeElasticParametersTimeSegy(SeismicParameters &seismic_parameters) {
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nt = seismic_parameters.seismicGeometry()->nt();

    double dt = seismic_parameters.seismicGeometry()->dt();

    NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();

    NRLib::Volume volume_time = seismic_parameters.seismicGeometry()->createTimeVolume();
    NRLib::StormContGrid vp_time_grid(volume_time, nx, ny, nt);
    NRLib::StormContGrid vs_time_grid(volume_time, nx, ny, nt);
    NRLib::StormContGrid rho_time_grid(volume_time, nx, ny, nt);

    NRLib::StormContGrid &vpgrid = seismic_parameters.vpGrid();
    NRLib::StormContGrid &vsgrid = seismic_parameters.vsGrid();
    NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();
    NRLib::StormContGrid &twtgrid = seismic_parameters.twtGrid();

    double t_min = toptime.Min();
    generateParameterGridForOutput(vpgrid, twtgrid, vp_time_grid, dt, t_min, toptime);
    generateParameterGridForOutput(vsgrid, twtgrid, vs_time_grid, dt, t_min, toptime);
    generateParameterGridForOutput(rhogrid, twtgrid, rho_time_grid, dt, t_min, toptime);

    printf("Write vp in time on Segy format\n");
    SEGY::writeSegy(vp_time_grid, prefix + "vp_time" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);

    printf("Write vs in time on Segy format\n");
    SEGY::writeSegy(vs_time_grid, prefix + "vs_time" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);

    printf("Write rho in time on Segy format\n");
    SEGY::writeSegy(rho_time_grid, prefix + "rho_time" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);

}

void SeismicOutput::writeExtraParametersTimeSegy(SeismicParameters &seismic_parameters) {
    printf("Write extra parameters in time on Segy format\n");
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nt = seismic_parameters.seismicGeometry()->nt();

    double dt = seismic_parameters.seismicGeometry()->dt();

    NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();

    NRLib::Volume volume_time = seismic_parameters.seismicGeometry()->createTimeVolume();
    NRLib::StormContGrid &twtgrid = seismic_parameters.twtGrid();
    std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();

    double tmin = toptime.Min();
    std::vector<NRLib::StormContGrid> extra_parameter_time_grid;
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        NRLib::StormContGrid extra_parameter_time_grid_temp(volume_time, nx, ny, nt);
        generateParameterGridForOutput((extra_parameter_grid)[i], twtgrid, extra_parameter_time_grid_temp, dt, tmin, toptime);
        extra_parameter_time_grid.push_back(extra_parameter_time_grid_temp);
    }
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        SEGY::writeSegy(extra_parameter_time_grid[i], prefix + extra_parameter_names[i] + "_time" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);
    }
}

void SeismicOutput::writeTimeSurfaces(SeismicParameters &seismic_parameters) {
    NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
    NRLib::RegularSurface<double> &bottime = seismic_parameters.bottomTime();

    printf("Write time surfaces on Storm format\n");
    bottime.WriteToFile(prefix + "bottime" + suffix + ".storm");
    toptime.WriteToFile(prefix + "toptime" + suffix + ".storm");
}

void SeismicOutput::writeElasticParametersDepthSegy(SeismicParameters &seismic_parameters) {
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nz = seismic_parameters.seismicGeometry()->nz();

    double dz = seismic_parameters.seismicGeometry()->dz();
    double z0 = seismic_parameters.seismicGeometry()->z0();

    NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();

    NRLib::Volume volume = seismic_parameters.seismicGeometry()->createDepthVolume();

    NRLib::StormContGrid &vpgrid = seismic_parameters.vpGrid();
    NRLib::StormContGrid &vsgrid = seismic_parameters.vsGrid();
    NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();
    NRLib::StormContGrid &zgrid = seismic_parameters.zGrid();

    NRLib::StormContGrid vp_depth_grid(volume, nx, ny, nz);
    NRLib::StormContGrid vs_depth_grid(volume, nx, ny, nz);
    NRLib::StormContGrid rho_depth_grid(volume, nx, ny, nz);

    generateParameterGridForOutput(vpgrid, zgrid, vp_depth_grid, dz, z0, toptime);
    generateParameterGridForOutput(vsgrid, zgrid, vs_depth_grid, dz, z0, toptime);
    generateParameterGridForOutput(rhogrid, zgrid, rho_depth_grid, dz, z0, toptime);
    printf("Write vp in depth on Segy format\n");
    SEGY::writeSegy(vp_depth_grid, prefix + "vp_depth" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, depth_window);

    printf("Write vs in depth on Segy format\n");
    SEGY::writeSegy(vs_depth_grid, prefix + "vs_depth" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, depth_window);

    printf("Write rho in depth on Segy format\n");
    SEGY::writeSegy(rho_depth_grid, prefix + "rho_depth" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, depth_window);
}

void SeismicOutput::writeExtraParametersDepthSegy(SeismicParameters &seismic_parameters) {
    printf("Write extra parameters in depth on Segy format\n");
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nz = seismic_parameters.seismicGeometry()->nz();
    double dz = seismic_parameters.seismicGeometry()->dz();
    NRLib::RegularSurface<double> &toptime = seismic_parameters.topTime();

    NRLib::Volume volume = seismic_parameters.seismicGeometry()->createDepthVolume();
    NRLib::StormContGrid &zgrid = seismic_parameters.zGrid();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();

    std::vector<NRLib::StormContGrid> &extra_parameter_grid = seismic_parameters.extraParametersGrids();
    std::vector<NRLib::StormContGrid> extra_parameter_depth_grid;

    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        NRLib::StormContGrid extra_parameter_depth_grid_temp(volume, nx, ny, nz);
        generateParameterGridForOutput((extra_parameter_grid)[i], zgrid, extra_parameter_depth_grid_temp, dz, toptime.Min(), toptime);
        extra_parameter_depth_grid.push_back(extra_parameter_depth_grid_temp);
    }
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        SEGY::writeSegy(extra_parameter_depth_grid[i], prefix + extra_parameter_names[i] + "_depth" + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, depth_window);
    }
}

void SeismicOutput::writeVpVsRho(SeismicParameters &seismic_parameters, bool use_window) {
    NRLib::StormContGrid &vpgrid = seismic_parameters.vpGrid();
    NRLib::StormContGrid &vsgrid = seismic_parameters.vsGrid();
    NRLib::StormContGrid &rhogrid = seismic_parameters.rhoGrid();

    if(use_window) {
        printf("Write elastic parameters on Storm format using window\n");
        std::string filename = prefix + "vp" + suffix + ".storm";
        STORM::writeStorm(vpgrid, filename, top_depth_window, bot_depth_window, depth_window);
        filename = prefix + "vs" + suffix + ".storm";
        STORM::writeStorm(vsgrid, filename, top_depth_window, bot_depth_window, depth_window);
        filename = prefix + "rho" + suffix + ".storm";
        STORM::writeStorm(rhogrid, filename, top_depth_window, bot_depth_window, depth_window);
    } else {
        printf("Write elastic parameters on Storm format not using window\n");
        std::string filename = prefix + "vp" + suffix + ".storm";
        vpgrid.WriteToFile(filename);
        filename = prefix + "vs" + suffix + ".storm";
        vsgrid.WriteToFile(filename);
        filename = prefix + "rho" + suffix + ".storm";
        rhogrid.WriteToFile(filename);
    }
}

void SeismicOutput::writeZValues(SeismicParameters &seismic_parameters, bool use_window) {
    NRLib::StormContGrid &zgrid = seismic_parameters.zGrid();
    std::string filename = prefix + "zgrid" + suffix + ".storm";


    if(use_window) {
        printf("Write zvalues on Storm format using window\n");
        STORM::writeStorm(zgrid, filename, top_depth_window, bot_depth_window, depth_window);
    } else {
        printf("Write zvalues on Storm format not using window\n");
        zgrid.WriteToFile(filename);
    }

}

void SeismicOutput::writeTwt(SeismicParameters &seismic_parameters, bool use_window) {
    NRLib::StormContGrid &twtgrid = seismic_parameters.twtGrid();
    std::string filename = prefix + "twt" + suffix + ".storm";

    if(use_window) {
        printf("Write two way time on Storm format using window\n");
        STORM::writeStorm(twtgrid, filename, top_depth_window, bot_depth_window, depth_window);
    } else {
        printf("Write two way time on Storm format not using window\n");
        twtgrid.WriteToFile(filename);
    }
}





void SeismicOutput::generateParameterGridForOutput(NRLib::StormContGrid &input_grid, NRLib::StormContGrid &time_or_depth_grid, NRLib::StormContGrid &output_grid, double delta_time_or_depth, double zero_time_or_depth, NRLib::RegularSurface<double> &toptime) {
    for (size_t i = 0; i < output_grid.GetNI(); i++) {
        for (size_t j = 0; j < output_grid.GetNJ(); j++) {
            double x, y, z;
            input_grid.FindCenterOfCell(i, j, 0, x, y, z);

            double topt = toptime.GetZ(x, y);
            if (!toptime.IsMissing(topt)) { //check whether there are values in input_grid in this pillar - if not, cells in output_grid will be zero
                double location = zero_time_or_depth + 0.5 * delta_time_or_depth;
                for (size_t k = 0; k < output_grid.GetNK(); k++) {
                    //find cell index in time or depth grid
                    size_t location_index = findCellIndex(i, j, location, time_or_depth_grid);
                    if (location_index == 999999) {          //if location is above all values in pillar of time_or_depth_grid,
                        location_index = input_grid.GetNK() - 1;    //output_grid is given the value of the bottom cell of input_grid
                    }
                    output_grid(i, j, k) = input_grid(i, j, location_index);
                    location += delta_time_or_depth;
                }
            } else {
                for (size_t k = 0; k < output_grid.GetNK(); k++) {
                    output_grid(i, j, k) = 0.0;
                }
            }
        }
    }
}

size_t SeismicOutput::findCellIndex(size_t i, size_t j, double target_k, NRLib::StormContGrid &grid) {
    size_t found_k = 999999;
    size_t nz = grid.GetNK();
    for (size_t k = 0; k < nz; k++) {
        if (grid(i, j, k) > target_k) {
            found_k = k;
            break;
        }
    }
    return found_k;
}

void SeismicOutput::writeSeismicTimeSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec) {
    printf("Write seismic in time on Segy format\n");
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    for (size_t i = 0; i < ntheta; i++) {
        SEGY::writeSegy(timegridvec[i], prefix + "seismic_time_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicTimeStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec) {
    printf("Write seismic in time on Storm format\n");
    ModelSettings *model_settings = seismic_parameters.model_settings;
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    for (size_t i = 0; i < ntheta; i++) {
        std::string filename = prefix + "seismic_time_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".storm";
        if (model_settings->GetOutputSeismicStackTimeStorm() || model_settings->GetOutputSeismicStackTimeSegy()) {
            STORM::writeStorm(timegridvec[i], filename, top_time_window, bot_time_window, time_window, true);
        } else {
            STORM::writeStorm(timegridvec[i], filename, top_time_window, bot_time_window, time_window);
        }
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicTimeshiftSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec) {
    printf("Write seismic shifted in time on Segy format\n");
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
    for (size_t i = 0; i < ntheta; i++) {
        SEGY::writeSegy(timeshiftgridvec[i], prefix + "seismic_timeshift_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicTimeshiftStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec) {
    printf("Write seismic shifted in time on Storm format\n");
    ModelSettings *model_settings = seismic_parameters.model_settings;
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    for (size_t i = 0; i < ntheta; i++) {
        std::string filename = prefix + "seismic_timeshift_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".storm";
        if (model_settings->GetOutputSeismicStackTimeShiftStorm() || model_settings->GetOutputSeismicStackTimeShiftSegy()) {
            STORM::writeStorm(timeshiftgridvec[i], filename, top_time_window, bot_time_window, time_window, true);
        } else {
            STORM::writeStorm(timeshiftgridvec[i], filename, top_time_window, bot_time_window, time_window);
        }
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicDepthSegy(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec) {
    printf("Write seismic in depth on Segy format\n");
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
    for (size_t i = 0; i < ntheta; i++) {
        SEGY::writeSegy(depthgridvec[i], prefix + "seismic_depth_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".segy", inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, depth_window);
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicDepthStorm(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec) {
    printf("Write seismic in depth on Storm format\n");
    ModelSettings *model_settings = seismic_parameters.model_settings;
    double theta = seismic_parameters.theta0();
    double ntheta = seismic_parameters.nTheta();
    double dtheta = seismic_parameters.dTheta();
    for (size_t i = 0; i < ntheta; i++) {
        std::string filename = prefix + "seismic_depth_" + NRLib::ToString(NRLib::Radian * theta) + suffix + ".storm";
        if (model_settings->GetOutputSeismicStackDepthStorm() || model_settings->GetOutputSeismicStackDepthSegy()) {
            STORM::writeStorm(depthgridvec[i], filename, top_depth_window, bot_depth_window, depth_window, true);
        } else {
            STORM::writeStorm(depthgridvec[i], filename, top_depth_window, bot_depth_window, depth_window);
        }
        theta = theta + dtheta;
    }
}

void SeismicOutput::writeSeismicStackTime(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timegridvec) {
    ModelSettings *model_settings = seismic_parameters.model_settings;
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nt = seismic_parameters.seismicGeometry()->nt();
    NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();
    size_t ntheta = seismic_parameters.nTheta();

    NRLib::StormContGrid stack_timegridvec(volume_t, nx, ny, nt);
    float ntheta_inv = static_cast<float>(1.0 / ntheta);
    for (size_t angle = 0; angle < ntheta; ++angle) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nt; ++k) {
                    stack_timegridvec(i, j, k) += ntheta_inv * timegridvec[angle](i, j, k);
                }
            }
        }
        timegridvec[angle] = NRLib::StormContGrid(0, 0, 0);
    }

    //print out seismic stack in time, in storm and segy
    if (model_settings->GetOutputSeismicStackTimeSegy()) {
        printf("Write seismic stack in time on Segy format\n");
        NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
        std::string filename = prefix + "seismic_time_stack" + suffix + ".segy";
        SEGY::writeSegy(stack_timegridvec, filename, inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);
    }
    if (model_settings->GetOutputSeismicStackTimeStorm()) {
        printf("Write seismic stack in time on Storm format\n");
        std::string filename = prefix + "seismic_time_stack" + suffix + ".storm";
        STORM::writeStorm(stack_timegridvec, filename, top_time_window, bot_time_window, time_window);
    }
}

void SeismicOutput::writeSeismicStackTimeshift(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &timeshiftgridvec) {
    ModelSettings *model_settings = seismic_parameters.model_settings;
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nt = seismic_parameters.seismicGeometry()->nt();
    NRLib::Volume volume_t = seismic_parameters.seismicGeometry()->createTimeVolume();
    size_t ntheta = seismic_parameters.nTheta();

    NRLib::StormContGrid stack_timeshiftgridvec(volume_t, nx, ny, nt);
    float ntheta_inv = static_cast<float>(1.0 / ntheta);
    for (size_t angle = 0; angle < ntheta; ++angle) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nt; ++k) {
                    stack_timeshiftgridvec(i, j, k) += ntheta_inv * timeshiftgridvec[angle](i, j, k);
                }
            }
        }
        timeshiftgridvec[angle] = NRLib::StormContGrid(0, 0, 0);
    }

    //print out seismic stack shifted in time, in storm and segy
    if (model_settings->GetOutputSeismicStackTimeShiftSegy()) {
        printf("Write seismic stack shifted in time on Segy format\n");
        NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
        std::string filename = prefix + "seismic_timeshift_stack" + suffix + ".segy";
        SEGY::writeSegy(stack_timeshiftgridvec, filename, inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_time_window, bot_time_window, time_window);
    }
    if (model_settings->GetOutputSeismicStackTimeShiftStorm()) {
        printf("Write seismic stack shifted in time on Storm format\n");
        std::string filename = prefix + "seismic_timeshift_stack" + suffix + ".storm";
        STORM::writeStorm(stack_timeshiftgridvec, filename, top_time_window, bot_time_window, time_window);
    }
}

void SeismicOutput::writeSeismicStackDepth(SeismicParameters &seismic_parameters, std::vector<NRLib::StormContGrid> &depthgridvec) {
    ModelSettings *model_settings = seismic_parameters.model_settings;
    size_t nx = seismic_parameters.seismicGeometry()->nx();
    size_t ny = seismic_parameters.seismicGeometry()->ny();
    size_t nz = seismic_parameters.seismicGeometry()->nz();
    NRLib::Volume volume = seismic_parameters.seismicGeometry()->createDepthVolume();
    size_t ntheta = seismic_parameters.nTheta();

    NRLib::StormContGrid stack_depthgrid(volume, nx, ny, nz);
    float ntheta_inv = static_cast<float>(1.0 / ntheta);
    for (size_t angle = 0; angle < ntheta; ++angle) {
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nz; ++k) {
                    stack_depthgrid(i, j, k) += ntheta_inv * depthgridvec[angle](i, j, k);
                }
            }
        }
        depthgridvec[angle] = NRLib::StormContGrid(0, 0, 0);
    }

    //print out seismic stack in depth, in storm and segy
    if (model_settings->GetOutputSeismicStackDepthSegy()) {
        printf("Write seismic stack in depth on Segy format\n");
        NRLib::SegyGeometry *segy_geometry = seismic_parameters.segyGeometry();
        std::string filename = prefix + "seismic_depth_stack" + suffix + ".segy";
        SEGY::writeSegy(stack_depthgrid, filename, inline_start, xline_start, xline_x_axis, inline_step, xline_step, segy_geometry, scalco, top_depth_window, bot_depth_window, time_window);
    }
    if (model_settings->GetOutputSeismicStackDepthStorm()) {
        printf("Write seismic stack in depth on Storm format\n");
        std::string filename = prefix + "seismic_depth_stack" + suffix + ".storm";
        STORM::writeStorm(stack_depthgrid, filename, top_depth_window, bot_depth_window, depth_window);
    }
}
