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

#ifndef NRLIB_SERIALIZE_HPP
#define NRLIB_SERIALIZE_HPP

#include <ostream>
#include <vector>

#include "nrlib/iotools/fileio.hpp"

namespace NRLib {
  template<typename T>
  void SerializeVector(const std::vector<T> & vec,
                       std::ostream & ostream);

  template<typename T>
  void Serialize2dVector(const std::vector<std::vector<T> > & vec,
                         std::ostream & ostream);

  /// Does not work for std::string containing whitespace characters
  template<typename T>
  void DeserializeVector(std::istream & istream, std::vector<T> & vec);

  template<typename T>
  void Deserialize2dVector(std::istream & istream, std::vector<std::vector<T> > & vec, char end_char);
}


/**** Implementations ****/

template<typename T>
void NRLib::SerializeVector(const std::vector<T> & vec,
                            std::ostream & ostream)
{
  NRLib::WriteAsciiArray(ostream,
                         vec.begin(),
                         vec.end(),
                         static_cast<int>(vec.size()));
}


template <typename T>
void NRLib::Serialize2dVector(const std::vector<std::vector<T> >& vec,
                              std::ostream& ostream)
{
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i].size() == 0)
      ostream << ' ' << '\n';
    else
      NRLib::WriteAsciiArray(ostream,
                             vec[i].begin(),
                             vec[i].end(),
                             static_cast<int>(vec[i].size()));
  }
}

template <typename T>
void NRLib::DeserializeVector(std::istream   & istream,
                              std::vector<T> & vec)
{
  std::string line;
  std::getline(istream, line, '\n');
  if (line.size() == 0)
    return;
  std::vector<std::string> tokens = GetTokens(line);
  for (size_t i = 0; i < tokens.size(); ++i) {
    vec.push_back(NRLib::ParseType<T>(tokens[i]));
  }
}

template <typename T>
void NRLib::Deserialize2dVector(std::istream& istream,
                                std::vector<std::vector<T> > & vec,
                                char end_char)
{
  while (istream.good()) {
    std::string line;
    std::getline(istream, line, '\n');
    if (line.size() == 0)
      continue;
    if (line[0] == end_char)
      break;
    std::vector<std::string> tokens = GetTokens(line);
    vec.push_back(std::vector<size_t>());
    for (size_t i = 0; i < tokens.size(); ++i) {
      vec[vec.size() - 1].push_back(NRLib::ParseType<T>(tokens[i]));
    }
  }
}


#endif
