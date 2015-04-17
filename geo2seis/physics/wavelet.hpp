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

    double GetPeakFrequency() {
        return peak_frequency_;
    }

    double GetDepthAdjustmentFactor() {
        return depth_adjustment_factor_;
    }

  private:

    double FindPeakFrequency(std::vector<double> wavelet, int sample_number_for_zero_time);

    double FindAbsMaxOfVector(std::vector<double> vector);

    double FindDepthAdjustmentFactor(std::vector<double> wavelet, double time_sampling_in_ms);

    void ResampleTrace(std::vector<double> &wavelet, std::vector<double> &wavelet_out, size_t scale_factor);

    std::string file_format_;
    int sample_number_for_zero_time_;
    double time_sampling_in_ms_;
    std::vector<double> wavelet_;
    std::vector<double> time_vector_;
    bool is_ricker_;
    double peak_frequency_;
    double depth_adjustment_factor_;

};

#endif
