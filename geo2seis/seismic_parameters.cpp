#include <seismic_parameters.hpp>

#include <nrlib/math/constants.hpp>
#include <physics/wavelet.hpp>
#include <nrlib/eclipsegrid/eclipsegrid.hpp>
#include <nrlib/surface/regularsurfacerotated.hpp>
#include <nrlib/segy/segygeometry.hpp>
#include <nrlib/stormgrid/stormcontgrid.hpp>
#include "seismic_geometry.hpp"


SeismicParameters::SeismicParameters(ModelSettings *model_settings) {
    this->model_settings = model_settings;

    seismic_geometry = new SeismicGeometry();
    seismic_output = new SeismicOutput(model_settings);

    calculateAngleSpan();

    segy_geometry = NULL;

    setupWavelet();
    readEclipseGrid();
    findGeometry();
    findSurfaceGeometry();

    createGrids();
}

void SeismicParameters::calculateAngleSpan() {
    theta_0 = model_settings->GetTheta0();
    dtheta = model_settings->GetDTheta();
    theta_max = model_settings->GetThetaMax();

    if (dtheta == 0) {
        ntheta = 1;
    } else {
        ntheta = size_t((theta_max - theta_0) / dtheta) + 1;
        dtheta = (theta_max - theta_0) / (ntheta - 1);
    }
}


void SeismicParameters::setupWavelet() {
    if (model_settings->GetRicker()) {
        double peakF = model_settings->GetPeakFrequency();
        _wavelet = new Wavelet(peakF);
    } else {
        std::string wavelet_file_format = model_settings->GetWaveletFileFormat();
        std::string wavelet_file_name = model_settings->GetWaveletFileName();
        _wavelet = new Wavelet(wavelet_file_name, wavelet_file_format);
    }
    _wavelet_scale = model_settings->GetWaveletScale();
}

void SeismicParameters::readEclipseGrid() {
    std::string filename = model_settings->GetEclipseFileName();

    printf("Start reading Eclipsegrid from file\n");
    eclipse_grid = new NRLib::EclipseGrid(filename);
    printf("Eclipsegrid read.\n");

    std::vector<std::string> names = model_settings->GetParameterNames();
    if (!eclipse_grid->HasParameter(names[0])) {
        std::cout << "Parameter " + names[0] + " is not found in Eclipse grid\n";
        exit(0);
    }
    if (!eclipse_grid->HasParameter(names[1])) {
        std::cout << "Parameter " + names[1] + " is not found in Eclipse grid\n";
        exit(0);
    }
    if (!eclipse_grid->HasParameter(names[2])) {
        std::cout << "Parameter " + names[2] + " is not found in Eclipse grid\n";
        exit(0);
    }
    std::vector<std::string> extra_parameter_names = model_settings->GetExtraParameterNames();
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        if (!eclipse_grid->HasParameter(extra_parameter_names[i])) {
            std::cout << "Parameter " + extra_parameter_names[i] + " is not found in Eclipse grid\n";
            exit(0);
        }
    }
}

void SeismicParameters::findGeometry() {
    seismic_geometry->setDxDy(model_settings->GetDx(), model_settings->GetDy());
    seismic_geometry->setDz(model_settings->GetDz());
    seismic_geometry->setDt(model_settings->GetDt());

    const NRLib::EclipseGeometry &geometry = eclipse_grid->GetGeometry();

    if (model_settings->GetAreaGiven()) {
        double x0 = model_settings->GetX0();
        double y0 = model_settings->GetY0();
        double lx = model_settings->GetLx();
        double ly = model_settings->GetLy();
        double angle = model_settings->GetAngle();
        seismic_geometry->setGeometry(x0, y0, lx, ly, angle);

    } else if (model_settings->GetAreaFromSurface() != "") {
        NRLib::RegularSurfaceRotated<double> toptime_rot = NRLib::RegularSurfaceRotated<double>(model_settings->GetAreaFromSurface());
        double x0 = toptime_rot.GetXRef();
        double y0 = toptime_rot.GetYRef();
        double lx = toptime_rot.GetLengthX();
        double ly = toptime_rot.GetLengthY();
        double angle = toptime_rot.GetAngle();
        seismic_geometry->setGeometry(x0, y0, lx, ly, angle);

    } else if (model_settings->GetAreaFromSegy() != "") {
        int scalcoloc = 71;
        NRLib::TraceHeaderFormat::coordSys_t coord = NRLib::TraceHeaderFormat::UTM;
        NRLib::TraceHeaderFormat *thf = new NRLib::TraceHeaderFormat(scalcoloc, model_settings->GetUtmxIn(), model_settings->GetUtmyIn(), model_settings->GetIL0In(), model_settings->GetXL0In(), coord);
        double z0 = 0.0;
        NRLib::Volume *volume = NULL;
        std::vector<NRLib::TraceHeaderFormat *> thfvec;
        thfvec.push_back(thf);

        NRLib::SegY segy(model_settings->GetAreaFromSegy(), static_cast<float>(z0), thfvec);
        segy.ReadAllTraces(volume, z0);
        segy.CreateRegularGrid();
        const NRLib::SegyGeometry *temp_segy_geometry = segy.GetGeometry();

        segy_geometry = new NRLib::SegyGeometry(temp_segy_geometry);
        segy_geometry->WriteGeometry();
        segy_geometry->WriteILXL();

        double x0 = temp_segy_geometry->GetX0();
        double y0 = temp_segy_geometry->GetY0();
        double lx = temp_segy_geometry->Getlx();
        double ly = temp_segy_geometry->Getly();
        double angle = temp_segy_geometry->GetAngle();
        double dx = temp_segy_geometry->GetDx();
        double dy = temp_segy_geometry->GetDy();

        seismic_geometry->setGeometry(x0, y0, lx, ly, angle);
        seismic_geometry->setDxDy(dx, dy);

    } else {
        double x0, y0, lx, ly, angle;
        geometry.FindEnclosingVolume(x0, y0, lx, ly, angle);
        seismic_geometry->setGeometry(x0, y0, lx, ly, angle);
    }

    if (model_settings->GetAreaGiven() && model_settings->GetAreaFromSurface() != "") {
        printf("WARNING! Area defined in two different ways. The area specified by the area command is used.\n");
    }
}

void SeismicParameters::findSurfaceGeometry() {
    const NRLib::EclipseGeometry &geometry = eclipse_grid->GetGeometry();

    double dx = seismic_geometry->dx();
    double dy = seismic_geometry->dy();

    double lxsurf = seismic_geometry->xsurfacelength();
    double lysurf = seismic_geometry->ysurfacelength();

    double xmin = seismic_geometry->xmin();
    double ymin = seismic_geometry->ymin();

    size_t nxsurfec = seismic_geometry->nxsurfaceeclipse();
    size_t nysurfec = seismic_geometry->nysurfaceeclipse();

    bool const_top_given = true;
    if (model_settings->GetTopTimeSurfaceFile() != "") {
        NRLib::RegularSurfaceRotated<double> top_time_rotated = NRLib::RegularSurfaceRotated<double>(model_settings->GetTopTimeSurfaceFile());
        double topmin = top_time_rotated.Min();
        top_time = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, topmin);
        top_time.SetMissingValue(top_time_rotated.GetMissingValue());
        for (size_t i = 0; i < top_time.GetNI(); i++) {
            for (size_t j = 0; j < top_time.GetNJ(); j++) {
                double x, y;
                top_time.GetXY(i, j, x, y);
                double value = top_time_rotated.GetZ(x, y);
                top_time(i, j) = value;
            }
        }

        bot_time = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, top_time.Max());
        const_top_given = false;
    } else {
        double t1 = model_settings->GetTopTimeConstant();
        top_time = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);
        bot_time = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, t1);

    }

    topeclipse = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);
    boteclipse = NRLib::RegularSurface<double>(xmin - dx, ymin - dy, lxsurf + 2 * dx, lysurf + 2 * dy, nxsurfec + 2, nysurfec + 2, -999.0);

    top_k = geometry.FindTopLayer();
    bottom_k = geometry.FindBottomLayer();

    seismic_geometry->setZReflectorCount(static_cast<size_t>(bottom_k + 2 - top_k));

    NRLib::Grid2D<double> values(nxsurfec + 2, nysurfec + 2, 0.0);
    //what? no interpolation on top?
    geometry.FindLayerSurface(values, top_k, 0, topeclipse.GetDX(), topeclipse.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
    for (size_t i = 0; i < topeclipse.GetNI(); i++) {
        for (size_t j = 0; j < topeclipse.GetNJ(); j++) {
            topeclipse(i, j) = values(i, j);
        }
    }

    if (model_settings->GetUseCornerpointInterpol()) {
        geometry.FindLayerSurfaceCornerpoint(values, bottom_k, 1, boteclipse.GetDX(), boteclipse.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
    } else {
        geometry.FindLayerSurface(values, bottom_k, 1, boteclipse.GetDX(), boteclipse.GetDY(), xmin - dx, ymin - dy, 0.0, 0);
    }

    for (size_t i = 0; i < boteclipse.GetNI(); i++) {
        for (size_t j = 0; j < boteclipse.GetNJ(); j++) {
            boteclipse(i, j) = values(i, j);
        }
    }

    if (model_settings->GetOutputDepthSurfaces()) {
        seismic_output->writeDepthSurfaces(topeclipse, boteclipse);
    }

    double d1 = topeclipse.Min();
    double d2 = boteclipse.Max();

    if (const_top_given) {
        double t1 = model_settings->GetTopTimeConstant();
        std::vector<double> const_vp = model_settings->GetConstVp();
        for (size_t i = 0; i < top_time.GetNI(); i++)
            for (size_t j = 0; j < top_time.GetNJ(); j++) {
                top_time(i, j) = t1 + 2000.0 * (topeclipse(i, j) - d1) / const_vp[0];
                bot_time(i, j) = top_time(i, j);
            }
    }

    topeclipse.Add(-1 * _wavelet->GetDepthAdjustmentFactor()); // add one wavelet length to bot and subtract from top
    boteclipse.Add(_wavelet->GetDepthAdjustmentFactor());
    d1 = topeclipse.Min();
    d2 = boteclipse.Max();

    seismic_geometry->setZRange(d1, d2);
}

void SeismicParameters::createGrids() {
    size_t nx = seismic_geometry->nx();
    size_t ny = seismic_geometry->ny();
    size_t nzrefl = seismic_geometry->zreflectorcount();

    NRLib::Volume volume = seismic_geometry->createDepthVolume();

    zgrid = new NRLib::StormContGrid(volume, nx, ny, nzrefl);
    vpgrid = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
    vsgrid = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
    rhogrid = new NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
    twtgrid = new NRLib::StormContGrid(volume, nx, ny, nzrefl);

    rgridvec = new std::vector<NRLib::StormContGrid>(ntheta);
    NRLib::StormContGrid rgrid(volume, nx, ny, nzrefl);

    std::vector<std::string> extra_parameter_names = model_settings->GetExtraParameterNames();

    std::vector<double> extra_parameter_default_values = model_settings->GetExtraParameterDefaultValues();
    extra_parameter_grid = new std::vector<NRLib::StormContGrid>(extra_parameter_names.size());
    for (size_t i = 0; i < extra_parameter_names.size(); ++i) {
        (*extra_parameter_grid)[i] = NRLib::StormContGrid(volume, nx, ny, nzrefl + 1);
    }

    std::vector<double> const_vp = model_settings->GetConstVp();
    std::vector<double> const_vs = model_settings->GetConstVs();
    std::vector<double> const_rho = model_settings->GetConstRho();

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t k = 0; k < nzrefl; k++) {
                (*zgrid)(i, j, k) = 0.0;
                (*vpgrid)(i, j, k) = static_cast<float>(const_vp[1]);
                (*vsgrid)(i, j, k) = static_cast<float>(const_vs[1]);
                (*rhogrid)(i, j, k) = static_cast<float>(const_rho[1]);
                (*twtgrid)(i, j, k) = 0.0;
                rgrid(i, j, k) = 0.0;

                for (size_t epi = 0; epi < extra_parameter_names.size(); ++epi) {
                    (*extra_parameter_grid)[epi](i, j, k) = static_cast<float>(extra_parameter_default_values[epi]);
                }
            }

            (*vpgrid)(i, j, nzrefl) = static_cast<float>(const_vp[2]);
            (*vsgrid)(i, j, nzrefl) = static_cast<float>(const_vs[2]);
            (*rhogrid)(i, j, nzrefl) = static_cast<float>(const_rho[2]);
            for (size_t epi = 0; epi < extra_parameter_names.size(); ++epi) {
                (*extra_parameter_grid)[epi](i, j, nzrefl) = 0.0;
            }
        }
    }
    for (size_t i = 0; i < ntheta; i++) {
        (*rgridvec)[i] = rgrid;
    }
}
