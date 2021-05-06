// $Id: waveletio.hpp 1380 2016-10-31 13:13:39Z perroe $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/// \file Parsers for several wavelet formats.

#ifndef NRLIB_WAVELET_WAVELETIO_HPP
#define NRLIB_WAVELET_WAVELETIO_HPP

#include <fstream>
#include <string>

namespace NRLib {

class Wavelet;

enum WaveletFileType {
  UnknownWavelet,
  JasonAscii,
  SyntoolOW,
  Strata,
  LandmarkAscii    ///< Landmark ASCII format, used by RokDoc and Syntool.
};

/// Find file type.
WaveletFileType FindWaveletFileType(const std::string& filename);

bool CheckFileIsJasonAscii(   const std::string& filename);
bool CheckFileIsSyntoolOW(    const std::string& filename);
bool CheckFileIsStrata(       const std::string& filename);
bool CheckFileIsLandmarkAscii(const std::string& filename);

void ReadJasonAsciiFile(   const std::string& filename, Wavelet& wavelet);
void ReadSyntoolOWFile(    const std::string& filename, Wavelet& wavelet);
void ReadStrataFile(       const std::string& filename, Wavelet& wavelet);
void ReadLandmarkAsciiFile(const std::string& filename, Wavelet& wavelet);

void WriteJasonAsciiFile(const std::string&              filename,
                         const Wavelet&                  wavelet,
                         const std::vector<std::string>& text_header);

// Helper functions, only for use in this file.
namespace NRLibPrivate {
  std::string JasonGetFirstNonCommentLine(std::ifstream& infile, int& line_number);
  void JasonParseUnits(const std::string & line,
                       std::string       & data_type,
                       std::string       & data_unit,
                       std::string       & depth_type,
                       std::string       & depth_unit);

  /// Get wavelet length that is large enough so that it can hold all the
  /// data in the file with center-index at nw/2. Always returns an odd number.
  /// \param[in] file_wavelet_length  Length of wavelet in file.
  size_t FindWaveletLength(size_t file_wavelet_length, size_t center_index);

  /// Find shift used to map from index in file to index in wavelet.
  /// \param[in] wavelet_length  Length of wavelet datastructure.
  size_t FindWaveletShift(size_t wavelet_length, size_t center_index);
} // namespace NRLibPrivate

} // namespace NRLib

#endif // NRLIB_WAVELET_WAVELETIO_HPP
