#include <utils/resampl_output.hpp>

#include <seismic_parameters.hpp>
#include <seismic_geometry.hpp>


ResamplOutput::ResamplOutput(SeismicParameters &seismic_parameters, bool time, size_t n_samples)
{
  size_t nx = seismic_parameters.GetSeismicGeometry()->nx();
  size_t ny = seismic_parameters.GetSeismicGeometry()->ny();
  n_samples_ = n_samples;
  time_ = time;

  NRLib::Volume volume;
  if (time)
    volume = seismic_parameters.GetSeismicGeometry()->createTimeVolume();
  else
    volume = seismic_parameters.GetSeismicGeometry()->createDepthVolume();
  seismic_parameters.GetSeismicOutput()->SetSegyGeometry(seismic_parameters, volume, nx, ny);
  segy_ok_ = seismic_parameters.GetSeismicOutput()->CheckUTMPrecision(seismic_parameters, volume, nx, ny);

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
    segy_files_ok_.push_back(seismic_parameters.GetSeismicOutput()->PrepareSegy(*(segy_files_[case_number]), time_or_depth_vec_reg, time_or_depth_vec_reg.size(), filename, seismic_parameters, dummy_vec, 1, time, false));
  else
    segy_files_ok_.push_back(false);

  NRLib::Grid2D<double> new_trace(n_samples_, 1, 0);
  traces_.push_back(new_trace);

  input_grid_.push_back(&input_grid);
}

//--------------------------------------------------------------------------------------------
void ResamplOutput::AddTrace(SeismicParameters                        & seismic_parameters,
                             std::vector<double>                      & time_or_depth_vec_reg,
                             const std::vector<NRLib::Grid2D<double>> & traces,
                             double                                     x,
                             double                                     y)
//--------------------------------------------------------------------------------------------
{
  std::vector<short> zero_vec(1, 0);
  for (size_t l = 0; l < traces.size(); ++l) {
    if (segy_files_ok_[l]) {
      seismic_parameters.GetSeismicOutput()->WriteSegyGather(traces[l], *(segy_files_[l]), time_or_depth_vec_reg, zero_vec, time_, x, y, false);
    }
  }
}
