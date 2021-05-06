// $Id: waveletio.cpp 1380 2016-10-31 13:13:39Z perroe $

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

#include "wavelet.hpp"
#include "waveletio.hpp"

#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"

#include <algorithm>
#include <fstream>
#include <string>

namespace NRLib {

WaveletFileType FindWaveletFileType(const std::string& filename)
{
  if (CheckFileIsJasonAscii(filename))
    return JasonAscii;
  if (CheckFileIsSyntoolOW(filename))
    return SyntoolOW;
  if (CheckFileIsStrata(filename))
    return Strata;
  if (CheckFileIsLandmarkAscii(filename))
    return LandmarkAscii;
  return UnknownWavelet;
}


bool CheckFileIsJasonAscii(const std::string& filename)
{
  std::ifstream file;
  OpenRead(file, filename);

  try {
    int line_number;
    std::string line = NRLibPrivate::JasonGetFirstNonCommentLine(file, line_number);
    std::string data_type, data_unit, depth_type, depth_unit;
    NRLibPrivate::JasonParseUnits(line, data_type, data_unit, depth_type, depth_unit);

    // Seems like a Jason ASCII file.
    return true;
  }
  catch (Exception& ) {
    // Format did not match, not Jason ASCII file.
    return false;
  }
}


bool CheckFileIsSyntoolOW(const std::string& filename)
{
  std::ifstream file;
  OpenRead(file, filename);

  // Line 3 should look like: "Wavelet Data:"
  std::string line3;
  std::getline(file, line3);
  std::getline(file, line3);
  std::getline(file, line3);
  if (line3.find("Wavelet Data") != std::string::npos)
    return true;
  return false;
}


bool CheckFileIsStrata(const std::string& filename)
{
  std::ifstream file;
  OpenRead(file, filename);

  // Line 1 should be #STRATA_WHIST
  std::string line1;
  std::getline(file, line1);
  if (line1.find("#STRATA_WHIST") != std::string::npos)
    return true;
  return false;
}


bool CheckFileIsLandmarkAscii(const std::string& filename)
{
  std::ifstream file;
  OpenRead(file, filename);

  // Line 1 should be "Landmark ASCII Wavelet"
  std::string line1;
  std::getline(file, line1);
  if (line1.find("Landmark ASCII Wavelet") != std::string::npos)
    return true;
  return false;
}


void ReadJasonAsciiFile(const std::string& filename, Wavelet& wavelet)
{
  int line_number = 0;

  std::ifstream infile;
  OpenRead(infile, filename);

  try {
    std::string line = NRLibPrivate::JasonGetFirstNonCommentLine(infile, line_number);

    std::string data_type, data_unit, depth_type, depth_unit;
    try {
      NRLibPrivate::JasonParseUnits(line, data_type, data_unit, depth_type, depth_unit);
    } catch (FileFormatError&) {
      std::string msg = "Error parsing " + filename + " , line "
                        + ToString(line_number) + ": Units line expected.";
      throw FileFormatError(msg);
    }

    double start  = ReadNext<double>(infile, line_number);
    double dt     = ReadNext<double>(infile, line_number);
    size_t center = static_cast<size_t>(-start/dt + 0.5);
    size_t nw     = ReadNext<size_t>(infile, line_number);

    wavelet.SetSampleInterval(dt);
    size_t data_length = NRLibPrivate::FindWaveletLength(nw, center);
    wavelet.Resize(data_length);
    size_t shift = NRLibPrivate::FindWaveletShift(data_length, center);

    for (size_t i = 0; i < nw; ++i) {
      wavelet[i + shift] = ReadNext<double>(infile, line_number);
    }

  } catch (EndOfFile&) {
    throw FileFormatError("Unexpected end of file.");
  } catch (Exception& e) {
    std::string msg = "Error parsing " + filename + ", line number "
                      + ToString(line_number) + ": " + e.what();
    throw FileFormatError(msg);
  }

  try {
    /// \todo Maybe just ignore this.
    std::string line = NRLibPrivate::JasonGetFirstNonCommentLine(infile, line_number);
    throw FileFormatError("Unexpected data at end of file.");
  } catch (EndOfFile&) {
    // OK, finished with file.
  }
}


void WriteJasonAsciiFile(const std::string&              filename,
                         const Wavelet&                  wavelet,
                         const std::vector<std::string>& text_header)
{
  std::ofstream out;
  NRLib::OpenWrite(out, filename);

  for (size_t i = 0; i < text_header.size(); ++i) {
    out << "\"* " + text_header[i] + "\n";
  }
  out << "\"*\n";
  out << "\"* File format:\n"
         "\"* - N lines starting with * are comment (such as this line)\n"
         "\"* - 1 line with four fields (data type, data unit, depth type, depth unit)\n"
         "\"* - 1 line with start time/depth\n"
         "\"* - 1 line with sample interval\n"
         "\"* - 1 line with number of data lines\n"
         "\"* - N lines with trace data\n"
         "\"* Data values are represented as floating point numbers,\n"
         "\"* except for lithology data, which is represented as strings\n"
         "\"wavelet\", \"none\", \"time\", \"ms\"\n";

  out << NRLib::ToString(wavelet.GetTimeForIndex(0)) + "\n";
  out << NRLib::ToString(wavelet.GetSampleInterval()) + "\n";
  out << NRLib::ToString(wavelet.GetNSamples()) + "\n";

  NRLib::WriteAsciiArray(out, wavelet.begin(), wavelet.end(), 1);
}


void ReadSyntoolOWFile(const std::string& filename, Wavelet& wavelet)
{
  int line_number = 0;
  std::ifstream infile;

  OpenRead(infile, filename);

  // Setting default values so we know if we didn't find the values we
  // were looking for.
  double n_samples = -1;
  double time_zero_sample = -1;
  double sample_interval = -1;
  bool end_of_header = false;
  std::string line;
  // Parsing header
  while(!end_of_header && std::getline(infile, line)) {
    ++line_number;
    if (line.find("Total Samples:") != std::string::npos) {
      if (n_samples != -1) {
        std::string msg = "Error parsing line " + ToString(line_number)
                          + ": Total Samples is set several times.";
        throw FileFormatError(msg);
      }
      std::string::size_type pos = line.find(':');
      n_samples = NRLib::ParseType<double>(NRLib::Chomp(line.substr(pos + 1)));
    }
    if (line.find("Sample Interval:") != std::string::npos) {
      if (sample_interval != -1) {
        std::string msg = "Error parsing  line " + ToString(line_number)
                          + ": Sample Interval is set several times.";
        throw FileFormatError(msg);
      }
      std::string::size_type pos = line.find(':');
      sample_interval = NRLib::ParseType<double>(NRLib::Chomp(line.substr(pos + 1)));
    }
    if (line.find("Time Zero Sample:") != std::string::npos) {
      if (time_zero_sample != -1) {
        std::string msg = "Error parsing  line " + ToString(line_number)
                          + ": Time Zero Sample is set several times.";
        throw FileFormatError(msg);
      }
      std::string::size_type pos = line.find(':');
      time_zero_sample = NRLib::ParseType<double>(NRLib::Chomp(line.substr(pos + 1)));
    }
    if (line.find("--------") != std::string::npos) {
       // End of header
      end_of_header = true;
    }
  }

  if (!end_of_header) {
    std::string msg = "Unexpected end of file.";
    throw FileFormatError(msg);
  }

  if (n_samples == -1) {
    std::string msg = "\"Total Samples\" line missing.";
    throw FileFormatError(msg);
  }

  if (sample_interval == -1) {
    std::string msg = "\"Sample Interval\" line missing.";
    throw FileFormatError(msg);
  }

  if (time_zero_sample == -1) {
    std::string msg = "\"Time Zero Sample\" line missing.";
    throw FileFormatError(msg);
  }

  size_t center = static_cast<size_t>(time_zero_sample);
  wavelet.SetSampleInterval(sample_interval);
  size_t nw     = static_cast<size_t>(n_samples);

  size_t data_length = NRLibPrivate::FindWaveletLength(nw, center);
  wavelet.Resize(data_length);
  size_t shift = NRLibPrivate::FindWaveletShift(data_length, center);

  try {
    for (size_t i = 0; i < nw; ++i) {
      double val;
      GetNextNonEmptyLine(infile, line_number, line);
      if (infile.eof())
        throw EndOfFile();
      std::vector<std::string> tokens = GetTokens(line);
      if (tokens.size() == 2) {
        try {
          ParseType<double>(tokens[0]); // z, not used.
          val = ParseType<double>(tokens[1]);
          wavelet[i + shift] = val;
        } catch (Exception&) {
          std::string msg = "Error parsing line " + ToString(line_number)
            + ": Line with z-value and wavelet value expected.";
          throw FileFormatError(msg);
        }
      }
      else {
        std::string msg = "Error parsing line " + ToString(line_number)
          + ": Line with z-value and wavelet value expected.";
        throw FileFormatError(msg);
      }
    }
  } catch (EndOfFile&) {
    std::string msg = "Unexpected end of file, to few samples.";
    throw FileFormatError(msg);
  }
}


void ReadStrataFile(const std::string& filename, Wavelet& wavelet)
{
  int line_number = 0;
  std::ifstream infile;
  OpenRead(infile, filename);

  bool end_of_header = false;
  std::string line;
  // Skip textual header
  while(!end_of_header && std::getline(infile, line)) {
    ++line_number;
    if (line.find("#STRATA_WPARAMS") != std::string::npos) {
      end_of_header = true;
    }
  }

  if (!end_of_header) {
    throw FileFormatError("Unexpected end of file.");
  }

  end_of_header = false;

  // Setting default values so we know if we didn't find the values we
  // were looking for.
  int n_samples = -1;
  int time_zero_sample = -1;
  double sample_interval = -1;

  // Parse wavelet parameters
  while(!end_of_header && std::getline(infile, line)) {
    ++line_number;
    if (line[0] != '~') {
      end_of_header = true;
    }
    else {
      std::string code = line.substr(1, 2);
      if (code == "SR") {
        std::string val = line.substr(3);
        sample_interval = ParseType<double>(val);
      }
      else if (code == "TZ") {
        std::string val = line.substr(3);
        time_zero_sample = ParseType<int>(val);
      }
      else if (code == "NS") {
        std::string val = line.substr(3);
        n_samples = ParseType<int>(val);
      }
      else if (code == "PR") {
        // Ignore
      }
      else {
        // Ignore
        /// @todo Print warning.
      }
    }
  }

  if (n_samples == -1) {
    throw FileFormatError("\"Total Samples\" line missing.");
  }

  if (sample_interval == -1) {
    throw FileFormatError("\"Sample Interval\" line missing.");
  }

  if (time_zero_sample == -1) {
    throw FileFormatError("\"Time Zero Sample\" line missing.");
  }

  size_t center = time_zero_sample - 1;
  wavelet.SetSampleInterval(sample_interval);

  size_t data_length = NRLibPrivate::FindWaveletLength(n_samples, center);
  wavelet.Resize(data_length);
  size_t shift = NRLibPrivate::FindWaveletShift(data_length, center);

  wavelet[shift] = ParseType<double>(line);

  for (int i = 1; i < n_samples; ++i) {
    wavelet[i + shift] = ReadNext<double>(infile, line_number);
  }

  if (std::getline(infile, line)) {
    if (line.find_first_not_of(Whitespace()) != std::string::npos) {
      throw FileFormatError("Unexpected data at end of file.");
    }
  }
}


void ReadLandmarkAsciiFile(const std::string& filename, Wavelet& wavelet)
{
  int line_number = 0;
  try {
    std::ifstream infile;
    OpenRead(infile, filename);

    std::string line;
    std::getline(infile, line);
    line_number++;

    if (line.find("Landmark ASCII Wavelet") == std::string::npos) {
      throw FileFormatError("File is not a Landmark ASCII Wavelet file.");
    }

    size_t nw     = ReadNext<size_t>(infile, line_number);
    size_t center = ReadNext<size_t>(infile, line_number);
    double dt     = ReadNext<double>(infile, line_number);

    wavelet.SetSampleInterval(dt);

    size_t data_length = NRLibPrivate::FindWaveletLength(nw, center);
    wavelet.Resize(data_length);
    size_t shift = NRLibPrivate::FindWaveletShift(data_length, center);

    for (size_t i = 0; i < nw; ++i) {
      wavelet[i + shift] = ReadNext<double>(infile, line_number);
    }
  }
  catch (EndOfFile&) {
    throw FileFormatError("Error reading Landmark ASCII Wavelet file " + filename
                          + ", Unexpected end of file.");
  }
  catch (IOError& e) {
    throw e;
  }
  catch (Exception& e) {
    throw FileFormatError("Error reading Landmark ASCII Wavelet file " + filename
                          + ", line " + ToString(line_number) + ": " + e.what());
  }
}


namespace NRLibPrivate {

std::string JasonGetFirstNonCommentLine(std::ifstream & infile,
                                        int           & line_number)
{
  std::string line;
  while(std::getline(infile, line)) {
    ++line_number;
    // Skip whitespace.
    std::string::size_type start = line.find_first_not_of(Whitespace());
    if (start == std::string::npos) {
      // Empty line
      continue;
    }
    if (line[start] == '\"') {
      std::string::size_type p
        = line.find_first_not_of(Whitespace(), start + 1);
      if (p != std::string::npos && line[p] == '*') {
        // Comment line
        continue;
      }
    }

    std::string::size_type end
      = line.find_last_not_of(Whitespace());

    std::string strippedLine = line.substr(start, end - start + 1);
    return strippedLine;
  }
  // End of file reached, or something bad happend while reading
  // from file.
  throw EndOfFile();
}


void JasonParseUnits(const std::string & line,
                     std::string       & data_type,
                     std::string       & data_unit,
                     std::string       & depth_type,
                     std::string       & depth_unit)
{
  // Data Type
  std::string::size_type pos1 = line.find(',');
  if (pos1 == std::string::npos) {
    throw FileFormatError();
  }
  std::string::size_type start = line.find('\"', 0);
  std::string::size_type end = line.rfind('\"', pos1);
  if (start == std::string::npos
      || end == std::string::npos
      || start >= end) {
    throw FileFormatError();
  }

  data_type = line.substr(start + 1, end - start - 1);

  // Data Unit
  std::string::size_type pos2 = line.find(',', pos1 + 1);
  if (pos2 == std::string::npos) {
    throw FileFormatError();
  }
  start = line.find('\"', pos1);
  end = line.rfind('\"', pos2);
  if (start == std::string::npos
      || end == std::string::npos
      || start >= end) {
    throw FileFormatError();
  }

  data_unit = line.substr(start + 1, end - start - 1);

  // Depth Type
  pos1 = pos2;
  pos2 = line.find(',', pos1 + 1);
  if (pos2 == std::string::npos) {
    throw FileFormatError();
  }
  start = line.find('\"', pos1);
  end = line.rfind('\"', pos2);
  if (start == std::string::npos
      || end == std::string::npos
      || start >= end) {
    throw FileFormatError();
  }

  depth_type = line.substr(start + 1, end - start - 1);

  // Depth Unit
  pos1 = pos2;

  start = line.find('\"', pos1);
  end = line.rfind('\"', std::string::npos);
  if (start == std::string::npos
      || end == std::string::npos
      || start >= end) {
    throw FileFormatError();
  }

  depth_unit = line.substr(start + 1, end - start - 1);
}

size_t FindWaveletLength(size_t file_wavelet_length, size_t center_index)
{
  assert(file_wavelet_length > center_index);

  size_t n_before = center_index;
  size_t n_after  = file_wavelet_length - center_index - 1;
  size_t length = 2 * std::max(n_before, n_after) + 1;

  return length;
}

size_t FindWaveletShift(size_t wavelet_length, size_t center_index)
{
  size_t new_center = wavelet_length/2;
  assert(new_center >= center_index);

  return new_center - center_index;
}


} // namespace NRLibPrivate

} // namespace NRLib
