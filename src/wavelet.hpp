#ifndef WAVELET_HPP
#define WAVELET_HPP

#include <string>
#include <vector>

class Wavelet
{
public:
  Wavelet(const double        peakF,
          const double        dt,
          const double        length,
          const double        length_factor,
          const bool          write_wavelet,
          const std::string & prefix);

  Wavelet(const std::string & filename,
          const std::string & file_format,
          const double        dt,
          const double        length,
          const double        length_factor,
          const bool          write_wavelet,
          const std::string & prefix,
          bool              & error);

  ~Wavelet(void);

  double FindWaveletPoint(double t);

  void   FindTwtLength(const std::vector<double> & wavelet,
                       const double                time_sampling_in_ms,
                       int                       & sample_number_for_zero_time,
                       double                    & twt_length);

  double GetTwtLength(void) { return twt_length_ ;}

private:

  void   FillTimeVector(std::vector<double> & time_vector,
                        const int             sample_number_peak,
                        const double          time_sampling_in_ms);

  void   FillWaveletVector(std::vector<double> & wavelet,
                           const double          time_sampling_in_ms,
                           const int             n);

  void   ReadLandMarkWavelet(const std::string   & filename,
                             std::vector<double> & wavelet,
                             double              & time_sampling_in_ms,
                             int                 & sample_number_for_zero_time);

  void   WriteLandMarkWavelet(const std::vector<double> & wavelet,
                              const int                   sample_number_zero_time,
                              const double                time_sampling,
                              const std::string         & filename) const;

  void   ResampleTrace(std::vector<double> & wavelet,
                       size_t                scale_factor);

  std::vector<double> wavelet_;
  std::vector<double> time_vector_;

  std::string         file_format_;

  double              time_sampling_in_ms_;
  double              twt_length_;               ///< Half a wavelet!

  double              peak_frequency_;

  bool                is_ricker_;
};

#endif
