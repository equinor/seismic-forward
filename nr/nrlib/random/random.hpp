// $Id: random.hpp 1189 2013-07-02 12:24:28Z anner $

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

#ifndef NRLIB_RANDOM_H
#define NRLIB_RANDOM_H

#include <string>

#include "dSFMT.h"

namespace NRLib {

/// Random generator class based on the Mersenne-Twister random
/// number generator.
/// Always initialize before use!
class Random {
public:
  ///Initializes with current time
  static void Initialize();

  static void Initialize(unsigned long seed);

  static void Initialize(const std::string& seed_file_);

  /// \return uniform number in [0,1)
  static double Unif01()             { return dsfmt_gv_genrand_close_open(); }

  /// \return uniform number in (0,1)
  static double Unif01Open()             { return dsfmt_gv_genrand_open_open(); }

  /// \return unsigned 32-bit integer betwen 0 and 0xFFFFFFFF
  static unsigned long DrawUint32()  { return dsfmt_gv_genrand_uint32(); }

  /// Marsaglia-Bray's method, see Ripley, p. 84.
  static double Norm01();

  /// Get start seed.
  static unsigned long GetStartSeed();

  /// Writes seed to file if seed-file is used.
  static void WriteSeedToFile();

private:
  /// Support function for Norm01
  static double g(double x);

  static void InitializeMT(unsigned long seed);

  static unsigned long start_seed_;

  static bool is_initialized_;

  static bool use_seed_file_;

  static std::string seed_file_;
};

}

#endif

