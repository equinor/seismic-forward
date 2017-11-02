#include "wavelet.hpp"

#include <iostream>
#include <fstream>
#include <list>
#include <complex>

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/fft/fft.hpp"

Wavelet::Wavelet(std::string filename, std::string file_format) : is_ricker_(false)
{

  if (NRLib::Uppercase(file_format) == "LANDMARK" || NRLib::Uppercase(file_format) == "LANDMARK ASCII WAVELET") {
    std::ifstream file;
    NRLib::OpenRead(file, filename);
    if (!file.bad()) {
      double number_of_samples;
      std::getline(file, file_format_);
      file >> number_of_samples >> sample_number_for_zero_time_ >> time_sampling_in_ms_;
      sample_number_for_zero_time_ = sample_number_for_zero_time_ - 1;
      size_t i = 0;
      double value;
      while (i < number_of_samples) {
        file >> value;
        wavelet_.push_back(value);
        ++i;
      }
    }
    FindTwtWavelet();

    std::vector<double> wavelet_out;
    size_t scale_factor = static_cast<size_t>(time_sampling_in_ms_);
    if (scale_factor < time_sampling_in_ms_) {
      scale_factor += 1;
    }
    ResampleTrace(wavelet_, wavelet_out, scale_factor);
    wavelet_ = wavelet_out;
    time_sampling_in_ms_ /= scale_factor;
    sample_number_for_zero_time_ *= static_cast<int>(scale_factor);


    time_vector_.reserve(wavelet_.size());
    double time;
    for (int j = 0; j < static_cast<int>(wavelet_.size()); ++j) {
      time = time_sampling_in_ms_ * (j - sample_number_for_zero_time_);
      time_vector_.push_back(time);
    }

    ////print wavelet
    //std::ofstream file_out;
    //filename = "wavelet_.txt";
    //NRLib::OpenWrite(file_out, filename);
    //for (size_t i = 0; i < wavelet_.size(); ++i){
    //  file_out << wavelet_[i];
    //  file_out << "\n";
    //}
  }
}

Wavelet::Wavelet(double peakF) :
file_format_(""),
  sample_number_for_zero_time_(0),
  time_sampling_in_ms_(0),
  wavelet_(std::vector<double>(0)),
  time_vector_(std::vector<double>(0)),
  is_ricker_(true),
  peak_frequency_(peakF) {

  twt_wavelet_ = 1000 / peak_frequency_;
}


Wavelet::~Wavelet() {
}


double Wavelet::FindAbsMaxOfVector(std::vector<double> vector) {
  double max_value = 0.0;
  for (size_t i = 0; i < vector.size(); ++i) {
    if (max_value < std::abs(vector[i])) {
      max_value = std::abs(vector[i]);
    }
  }
  return max_value;
}

void Wavelet::FindTwtWavelet()
{
  size_t start = 0;
  size_t end = wavelet_.size() - 1;
  double wavelet_max = FindAbsMaxOfVector(wavelet_);
  for (size_t i = 0; i < wavelet_.size(); ++i) {
    if (std::abs(wavelet_[i]) > wavelet_max * 0.01) {
      start = i;
      break;
    }
  }
  for (size_t i = wavelet_.size() - 1; i >= 0; --i) {
    if (std::abs(wavelet_[i]) > wavelet_max * 0.01) {
      end = i;
      break;
    }
  }
  double w1    = (sample_number_for_zero_time_ - start) * time_sampling_in_ms_;
  double w2    = (end -   sample_number_for_zero_time_) * time_sampling_in_ms_;
  std::cout << "w1, w2, start_i, end_i, sample_zero = " << w1 << " " << w2 << " " << start << " " << end << " " << sample_number_for_zero_time_ << "\n";
  twt_wavelet_ = std::max(w1, w2);
}


double Wavelet::FindWaveletPoint(double t)
{
  if (is_ricker_) {
    double rickerConst = NRLib::Pi * NRLib::Pi * peak_frequency_ * peak_frequency_ * 1e-6;
    double c = rickerConst * t * t;
    return (1 - 2 * c) * exp(-c);
  }
  else {
    if (wavelet_.size() > 0 && time_sampling_in_ms_ > 0) {

      size_t i;
      if (t < time_vector_[0]) {
        i = 0;
      }
      else {
        double start = (t - time_vector_[0]) / time_sampling_in_ms_;
        i = static_cast<size_t>(start);
        if (i < wavelet_.size() - 1 && t > time_vector_[i]) {
          ++i;
        }
      }
      if (i > wavelet_.size() - 1) {
        return 0;
      }
      else if (i > 0) {
        double a = (time_vector_[i] - t) / (time_vector_[i] - time_vector_[i - 1]);
        return a * wavelet_[i - 1] + (1 - a) * wavelet_[i];
      }
      else {
        return 0;
      }
    }
    else {
      return 0;
    }
  }
}


void Wavelet::ResampleTrace(std::vector<double> &wavelet, std::vector<double> &wavelet_out, size_t scale_factor)
{
  //
  // Transform to Fourier domain
  //
  std::vector<std::complex<double> > data_fft;
  NRLib::ComputeFFT1D(wavelet, data_fft, false, 0);

  //
  // Fill fine-sampled grid
  //
  std::vector<std::complex<double> > fine_data_fft(data_fft.size() * scale_factor);
  for (size_t i = 0; i < data_fft.size() / 2 + 1; i++) {
    fine_data_fft[i] = data_fft[i];
  }

  //Pad with zeros

  for (size_t i = data_fft.size() / 2 + 1; i < data_fft.size() * scale_factor; i++) {
    fine_data_fft[i] = 0;
  }

  //
  // Fine-sampled grid: Fourier --> Time
  //
  wavelet_out.resize(fine_data_fft.size());
  NRLib::ComputeFFTInv1D(fine_data_fft, wavelet_out, false);

  // Scale wavelet out according to change of length

  wavelet_out.resize(fine_data_fft.size() - scale_factor + 1);
  for (size_t i = 0; i < wavelet_out.size(); ++i) {
    wavelet_out[i] *= scale_factor;
  }
}


