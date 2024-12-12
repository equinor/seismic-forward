#include <utils/resampl_output.hpp>

#include <seismic_parameters.hpp>
#include <seismic_geometry.hpp>

//--------------------------------------------------------------------------------------
ResamplOutput::ResamplOutput(const SeismicParameters            & seismic_parameters,
                             const std::vector<double>          & time_or_depth_vec_reg,
                             const std::vector<std::string>     & filenames,
                             const bool                           segy_ok,
                             const bool                           time,
                             const size_t                         n_samples)
//--------------------------------------------------------------------------------------
  : segy_ok_  ( segy_ok   ),
    n_samples_( n_samples ),
    time_     ( time      )
{
  for (int i = 0 ; i < 10 ; i++) {
    NRLib::SegY * segy = new NRLib::SegY();
    segy_files_.push_back(segy);
  }

  for (size_t i = 0; i < filenames.size(); ++i) {
    std::vector<double> dummy_vec(1, 0.0);
    if (segy_ok_) {
      segy_files_ok_.push_back(seismic_parameters.GetSeismicOutput()->PrepareSegy(*(segy_files_[i]),
                                                                                  time_or_depth_vec_reg,
                                                                                  time_or_depth_vec_reg.size(),
                                                                                  filenames[i],
                                                                                  seismic_parameters.GetSegyGeometry(),
                                                                                  dummy_vec,
                                                                                  1,
                                                                                  time,
                                                                                  false));
    }
    else {
      segy_files_ok_.push_back(false);
    }
  }
}

//---------------------------------
ResamplOutput::~ResamplOutput(void)
//---------------------------------
{
  for (int i = 0 ; i < 10 ; i++)
    delete segy_files_[i];
}

//--------------------------------------------------------------------------------------------
void ResamplOutput::AddTrace(SeismicOutput                            * seismic_output,
                             const std::vector<double>                & time_or_depth_vec_reg,
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
