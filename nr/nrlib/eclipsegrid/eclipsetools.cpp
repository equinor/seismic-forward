// $Id: eclipsetools.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "eclipsetools.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cctype>

using namespace std;

namespace NRLib {
void SkipComments( std::ifstream& in_file )
{
  char c;
  std::string line_removed;

  bool comment = true;

  while (comment == true){
    in_file.get(c);
    if ( isspace(c) ){

    }
    else if (c == '-'){
      in_file.get(c);
      if (c == '-'){
        std::getline(in_file, line_removed);
      }
      else{
        in_file.unget();
        in_file.unget();
        comment = false;
      }
    }
    else{
      in_file.unget();
      comment = false;
    }
  }
}


std::string ReadParameterBuffer(std::ifstream& in_file)
{
  std::streamoff start_pos = in_file.tellg();
  const size_t K = 1000000; //how much the buffer size increases
  std::string buffer;
  size_t pos = buffer.find("/");
  size_t counter = 0;
  while ( pos == buffer.npos )
  {
    counter+=K;
    buffer.resize(counter);
    in_file.read(&(buffer[counter-K]), K);
    pos = buffer.find("/", counter-K);
  }
  if ( in_file.eof() ){
    in_file.clear();
  }

  in_file.seekg(start_pos + pos + 1);
  buffer.resize(pos);
  return buffer;
}


} // namespace NRLib

