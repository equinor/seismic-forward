// $Id: randomgenerator.hpp 1189 2013-07-02 12:24:28Z anner $

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

#ifndef NRLIB_RANDOMGENERATOR_H
#define NRLIB_RANDOMGENERATOR_H

#include "dSFMT.h"

namespace NRLib {

class RandomGenerator {
public:
  RandomGenerator();

  ~RandomGenerator();
  ///Initializes with current time
  void Initialize();

  void Initialize(unsigned long seed);

  /// \return unsigned 32-bit integer betwen 0 and 0xFFFFFFFF
  unsigned long DrawUint32()  { return dsfmt_genrand_uint32(&dsfmt); }

  /// \return uniform number in [0,1)
  double Unif01()             { return dsfmt_genrand_close_open(&dsfmt); }

  /// \return uniform number in (0,1)
  double Unif01Open()             { return dsfmt_genrand_open_open(&dsfmt); }

  /// Marsaglia-Bray's method, see Ripley, p. 84.
  double Norm01();

  /// Get start seed.
  unsigned long GetStartSeed();

private:
  /// RNG state
  dsfmt_t dsfmt;

  /// Support function for Norm01
  double g(double x);

  void InitializeMT(unsigned long seed);

  unsigned long start_seed_;

  bool is_initialized_;
};

}

#endif // NRLIB_RANDOMGENERATOR_H
