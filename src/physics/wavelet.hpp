#ifndef WAVELET_HPP
#define WAVELET_HPP

#include <string>
#include <vector>

class Wavelet {
  public:
    //Constructor
    Wavelet(std::string filename, std::string file_format);

    Wavelet(double peakF);

    //Destructor
    ~Wavelet();

    double FindWaveletPoint(double t);

    void FindTwtWavelet();

    double GetTwtWavelet() {
      return twt_wavelet_;
    }

  private:

    double FindAbsMaxOfVector(std::vector<double> vector);

    void ResampleTrace(std::vector<double> &wavelet, std::vector<double> &wavelet_out, size_t scale_factor);

    std::vector<double> wavelet_;
    std::vector<double> time_vector_;
    std::string file_format_;
    int         sample_number_for_zero_time_;
    double      time_sampling_in_ms_;
    bool        is_ricker_;
    double      peak_frequency_;
    double      twt_wavelet_;
};

#endif
