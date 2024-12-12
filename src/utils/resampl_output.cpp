#include <utils/resampl_output.hpp>

#include <seismic_parameters.hpp>
#include <seismic_geometry.hpp>

ResamplOutput::ResamplOutput(SeismicParameters & seismic_parameters,
                             bool                time,
                             size_t              n_samples)
{
  n_samples_ = n_samples;
  time_      = time;

  SeismicOutput       * seismic_output   = seismic_parameters.GetSeismicOutput();
  SeismicGeometry     * seismic_geometry = seismic_parameters.GetSeismicGeometry();
  NRLib::SegyGeometry * segy_geometry    = seismic_parameters.GetSegyGeometry();

  size_t nx  = seismic_geometry->nx();
  size_t ny  = seismic_geometry->ny();

  NRLib::Volume volume;

  if (time)
    volume = seismic_geometry->createTimeVolume();
  else
    volume = seismic_geometry->createDepthVolume();

  if (segy_geometry == NULL){
    NRLib::SegyGeometry * geometry = seismic_output->CreateSegyGeometry(volume, nx, ny);
    seismic_parameters.SetSegyGeometry(geometry);
    delete geometry;
  }

  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nSegy geometry:\n");
  segy_geometry->WriteGeometry();
  segy_geometry->WriteILXL();

  segy_ok_ = seismic_output->CheckUTMPrecision(segy_geometry, volume, nx, ny);

  segy_files_.push_back(&segy_1_);
  segy_files_.push_back(&segy_2_);
  segy_files_.push_back(&segy_3_);
  segy_files_.push_back(&segy_4_);
  segy_files_.push_back(&segy_5_);
  segy_files_.push_back(&segy_6_);
  segy_files_.push_back(&segy_7_);
  segy_files_.push_back(&segy_8_);
  segy_files_.push_back(&segy_9_);
  segy_files_.push_back(&segy_10_);
}

void ResamplOutput::AddResampleCase(std::string            filename,
                                    NRLib::StormContGrid & input_grid,
                                    bool                   time,
                                    std::vector<double>  & time_or_depth_vec_reg,
                                    SeismicParameters    & seismic_parameters)
{
  std::vector<double> dummy_vec(1);
  dummy_vec[0] = 0;
  size_t case_number = traces_.size();
  if (segy_ok_)
    segy_files_ok_.push_back(seismic_parameters.GetSeismicOutput()->PrepareSegy(*(segy_files_[case_number]),
                                                                                time_or_depth_vec_reg,
                                                                                time_or_depth_vec_reg.size(),
                                                                                filename,
                                                                                seismic_parameters,
                                                                                dummy_vec,
                                                                                1,
                                                                                time,
                                                                                false));
  else
    segy_files_ok_.push_back(false);

  NRLib::Grid2D<double> new_trace(n_samples_, 1, 0);
  traces_.push_back(new_trace);

  input_grid_.push_back(&input_grid);
}

//--------------------------------------------------------------------------------------------
void ResamplOutput::AddTrace(SeismicOutput                            * seismic_output,
                             std::vector<double>                      & time_or_depth_vec_reg,
                             const std::vector<NRLib::Grid2D<double>> & traces,
                             const double                               x,
                             const double                               y,
                             const bool                                 time)
//--------------------------------------------------------------------------------------------
{
  std::vector<short> zero_vec(1, 0);
  for (size_t l = 0; l < traces.size(); ++l) {
    if (segy_files_ok_[l]) {
      seismic_output->WriteSegyGather(traces[l],
                                      *(segy_files_[l]),
                                      time_or_depth_vec_reg,
                                      zero_vec,
                                      time,
                                      x,
                                      y,
                                      false);
    }
  }
}
