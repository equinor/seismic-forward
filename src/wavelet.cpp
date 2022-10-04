#include "nrlib/surface/regularsurface.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/fileio.hpp"

#include "nrlib/fft/fft.hpp"

#include "wavelet.hpp"

#include <iostream>
#include <fstream>
#include <complex>
#include <list>


//--------------------------------------------------
Wavelet::Wavelet(const double         peakF,
                 const double         dt,
                 const double         length,
                 const bool           write_wavelet,
                 const std::string  & prefix)
//--------------------------------------------------
  : wavelet_            ( std::vector<double>(0) ),
    time_vector_        ( std::vector<double>(0) ),
    file_format_        ( ""                     ),
    time_sampling_in_ms_( -999.0                 ),
    twt_length_         ( -999.0                 ),
    peak_frequency_     ( peakF                  ),
    is_ricker_          ( true                   )
{
  //
  // EVERYTHING below is CURRENTLY for wavelet export only. Wavelet is calculated on the fly!
  //
  double      dummy_length       = 1000.0/peakF;
  size_t      n                  = static_cast<size_t>(floor(dummy_length));
  int         sample_number_peak = n;
  std::string filename           = "Wavelet_" + NRLib::ToString(peakF) + "Hz_as_used.txt";
  if (prefix != "") {
    filename = prefix + "_" + filename;
  }

  FillWaveletVector(wavelet_,       // This is not the wavelet used in final calculations
                    dt,
                    n);

  if (length == -999.0)
    FindTwtLength(wavelet_,
                  dt,
                  sample_number_peak, // Also find peak position
                  twt_length_);       // Size of a half wavelet
  else
    twt_length_ = 0.5*length;

  if (write_wavelet) {
    WriteLandMarkWavelet(wavelet_,
                         sample_number_peak,
                         dt,
                         filename);
  }
}

//-----------------------------------------------
Wavelet::Wavelet(const std::string & filename,
                 const std::string & file_format,
                 const double        dt,
                 const double        length,
                 const bool          write_wavelet,
                 const std::string & prefix,
                 bool              & error)
//-----------------------------------------------
  : wavelet_            ( std::vector<double>(0) ),
    time_vector_        ( std::vector<double>(0) ),
    time_sampling_in_ms_( -999.0                 ),
    twt_length_         ( -999.0                 ),
    peak_frequency_     ( -999.0                 ),
    is_ricker_          ( false                  )
{
  if (NRLib::Uppercase(file_format) == "LANDMARK" || NRLib::Uppercase(file_format) == "LANDMARK ASCII WAVELET") {

    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nLandmark format assumed for wavelet\n");

    int sample_number_peak = -999;

    ReadLandMarkWavelet(filename,
                        wavelet_,
                        time_sampling_in_ms_,
                        sample_number_peak);

    if (length == -999.0)
      FindTwtLength(wavelet_,
                    time_sampling_in_ms_,
                    sample_number_peak,
                    twt_length_);         // Size of a half wavelet
    else
      twt_length_ = 0.5*length;

    //
    // Export wavelet as read from file - for debugging
    //
    //WriteLandMarkWavelet(wavelet_,
    //                     sample_number_peak,
    //                     time_sampling_in_ms_,
    //                     "WAVELET_from_file_as_read.txt");

    //
    // Scale wavelet if needed
    //
    size_t scale_factor = static_cast<size_t>(time_sampling_in_ms_);
    if (scale_factor < time_sampling_in_ms_) {
      scale_factor += 1;
    }
    time_sampling_in_ms_ /= scale_factor;

    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nResampling wavelet from %.1f ms to %.1f ms\n",
                                time_sampling_in_ms_*scale_factor,time_sampling_in_ms_);


    ResampleTrace(wavelet_, scale_factor);

    sample_number_peak *= static_cast<int>(scale_factor);

    //
    // Export resampled wavelet
    //
    if (write_wavelet) {
      std::string filename2 = std::string("Wavelet_resampled_to_") + NRLib::ToString(time_sampling_in_ms_, 1) + "ms.txt";
      if (prefix != "") {
        filename2 = prefix + "_" + filename2;
      }
      WriteLandMarkWavelet(wavelet_,
                           sample_number_peak,
                           time_sampling_in_ms_,
                           filename2);
    }

    //for (size_t ii = 0 ; ii < wavelet_.size() ; ii++) {
    //  printf("wavelet: %10.5f    time: %10.5f\n",wavelet_[ii],time_vector_[ii]);
    //}

    //
    // Fill time vector
    //
    time_vector_.clear();
    time_vector_.resize(wavelet_.size());

    FillTimeVector(time_vector_,
                   sample_number_peak,
                   time_sampling_in_ms_);
  }
  else {
    error = true;
  }
}

//---------------------
Wavelet::~Wavelet(void)
//---------------------
{
}

//---------------------------------------------------------------------
void Wavelet::FillTimeVector(std::vector<double> & time_vector,
                             const int             sample_number_peak,
                             const double          dt)
//---------------------------------------------------------------------
{
  double time;
  for (int i = 0 ; i < time_vector.size() ; ++i) {
    time = dt * (i - sample_number_peak);
    time_vector[i] = time;
  }
}

//------------------------------------------------------------------------
void Wavelet::FillWaveletVector(std::vector<double> & wavelet,
                                const double          time_sampling_in_ms,
                                const int             n)
//------------------------------------------------------------------------
{
  wavelet_.resize(2*n + 1);
  wavelet[n] = FindWaveletPoint(0.0);
  for (size_t i = 1 ; i < n + 1 ; i++) {
    double t   = static_cast<double>(i)*time_sampling_in_ms;
    double amp = FindWaveletPoint(t);
    wavelet[n - i] = amp;
    wavelet[n + i] = amp;
  }
}

//----------------------------------------------------------------------------------
void Wavelet::ReadLandMarkWavelet(const std::string   & filename,
                                  std::vector<double> & wavelet,
                                  double              & time_sampling_in_ms,
                                  int                 & sample_number_peak)
//----------------------------------------------------------------------------------
{
  std::ifstream file;
  NRLib::OpenRead(file, filename);
  if (!file.bad()) {
    double number_of_samples;
    std::getline(file, file_format_);
    file >> number_of_samples >> sample_number_peak >> time_sampling_in_ms;
    sample_number_peak -= 1;
    size_t i = 0;
    double value;
    while (i < number_of_samples) {
      file >> value;
      wavelet.push_back(value);
      ++i;
    }
  }
}

//-----------------------------------------------------------------------------------------
void Wavelet::WriteLandMarkWavelet(const std::vector<double> & wavelet,
                                   const int                   sample_number_peak,
                                   const double                time_sampling_in_ms,
                                   const std::string         & filename) const
//-----------------------------------------------------------------------------------------
{
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nExporting Landmark ASCII wavelet \'%s\'",filename.c_str());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Wavelet size   : %3d",wavelet.size());
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Peak at sample : %3d",sample_number_peak);
  NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n  Density in ms  : %.1f\n",time_sampling_in_ms);

  std::ofstream file_out;
  NRLib::OpenWrite(file_out, filename);
  file_out << "Landmark ASCII Wavelet\n"
           << wavelet.size()              << "\n"
           << sample_number_peak << "\n";
  file_out << std::fixed
           << std::setprecision(2)
           << time_sampling_in_ms         << "\n";
  for (size_t i = 0 ; i < wavelet.size() ; ++i){
    file_out << std::setprecision(6) << wavelet[i] << "\n";
  }
}

//----------------------------------------------------------------------------------
void Wavelet::FindTwtLength(const std::vector<double> & wavelet,
                            const double                time_sampling_in_ms,
                            int                       & sample_number_peak,
                            double                    & twt_length)
//----------------------------------------------------------------------------------
{
  size_t n = wavelet.size();

  double factor = 0.001;  // Fraction of peak amplitude where wavelet is considered zero
  double w_max  = 0.0;
  int    i_max  =  -1;

  for (size_t i = 0 ; i < n ; ++i) {
    if (w_max < std::abs(wavelet[i])) {
      w_max = std::abs(wavelet[i]);
      i_max = i;
    }
  }

  if (sample_number_peak != -999 && i_max != sample_number_peak) {
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\nWARNING: Sample incorrect sample number for peak amplitude in file?");
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n   Peak sample number from file   : %3d", sample_number_peak);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n   Peak sample number calculated  : %3d", i_max);
    NRLib::LogKit::LogFormatted(NRLib::LogKit::Low, "\n   Calculated sample number will be used\n");
    sample_number_peak = i_max;
  }

  double threshold = w_max*factor;  // Amplitude that defines length of wavelet
  size_t start     = 0;
  size_t end       = n - 1;

  for (size_t i = 0 ; i < n ; ++i) {
    if (std::abs(wavelet[i]) > threshold) {
      start = i;
      break;
    }
  }
  for (int i = n - 1; i >= 0; --i) {
    if (std::abs(wavelet[i]) > threshold) {
      end = i;
      break;
    }
  }

  double w1  = (i_max - start) * time_sampling_in_ms;
  double w2  = (end   - i_max) * time_sampling_in_ms;

 twt_length = std::max(w1, w2);
}

//----------------------------------------
double Wavelet::FindWaveletPoint(double t)
//----------------------------------------
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

//-------------------------------------------------------------
void Wavelet::ResampleTrace(std::vector<double> & wavelet,
                            size_t                scale_factor)
//-------------------------------------------------------------
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
  std::vector<double> wavelet_out(fine_data_fft.size());
  NRLib::ComputeFFTInv1D(fine_data_fft, wavelet_out, false);

  // Scale wavelet out according to change of length

  wavelet_out.resize(fine_data_fft.size() - scale_factor + 1);
  for (size_t i = 0; i < wavelet_out.size(); ++i) {
    wavelet_out[i] *= scale_factor;
  }

  std::swap(wavelet, wavelet_out);
}
