// $Id: commonheaders.hpp 883 2011-09-26 09:17:05Z perroe $

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

#ifndef COMMONHEADER_HPP
#define COMMONHEADER_HPP

#include <cstdio>
#include <iostream>
#include <fstream>
#include "segy.hpp"

namespace NRLib {

class SegyGeometry;

  /// The textual header of a SEG Y file.
  /// Also known as the EBCDIC header.
class TextualHeader
{
public:
  /// Default constructor.
  TextualHeader();

  /// The old standard header.
  static TextualHeader standardHeader();

  /// Set a line in the textual header.
  /// \param[in] lineNo Line number in interval [0, 39].
  /// \param[in] text   Text for the given line. ASCII encoded. Max 75 char.
  /// \return 0 on success, -1 if error in input, 1 if the text was truncated.
  int SetLine(int lineNo, const std::string& text);

  /// Write header to file.
  /// \param[in] file  Output file.
  void Write(std::ostream& file) const;

private:
  /// output buffer in EBCDIC encoding.
  std::string buffer_;
  //const char * buffer_;
};


class BinaryHeader
{
public:
  /// Constructor
  BinaryHeader(std::istream& file);
  // constructor for writing
  BinaryHeader();
  /// Update variables
  void Update(std::istream& file);
  void Write(std::ostream& file, SegyGeometry *geometry, double dz, size_t nz) const;
  short GetFormat() {return(format_);}
  int GetLino() {return(lino_);}
  short GetHns() {return(hns_);}
  short GetHdt() {return(hdt_);}
private:
  short format_;
  int lino_;
  short hns_;
  short hdt_;
};

} // namespace NRLib

#endif
