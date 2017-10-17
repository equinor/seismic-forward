// $Id: posteriormultinormal.hpp 1241 2014-02-18 12:17:26Z anner $

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

#ifndef NRLIB_STATISTICS_POSTERIORMULTINORMAL_HPP
#define NRLIB_STATISTICS_POSTERIORMULTINORMAL_HPP

#include <cstdlib>
#include <vector>


#include "../flens/nrlib_flens.hpp"
#include "../variogram/variogram.hpp"

namespace NRLib {

  class Point;

  bool CondDistrMultiNormal(const std::vector<double> &x_known,
                            const std::vector<double> &known_values,
                            const std::vector<double> &x_unknown,
                            const std::vector<double> &exp_known,
                            const std::vector<double> &exp_unknown,
                            const std::vector<double> &std_known,
                            const std::vector<double> &std_unknown,
                            const Variogram           &vario,
                            Vector                    &exp,
                            SymmetricMatrix           &cov);

  bool SecondayCondDistrMultiNormal2D(const std::vector<NRLib::Point> &pos_known,
                                      const std::vector<double>       &known_values,
                                      const std::vector<double>       &local_corr,   //Allows secondary observations. Must have same correlation with unknown. A value of 1 indicates primary observation, lower value seconday.
                                      const std::vector<NRLib::Point> &pos_unknown,
                                      const std::vector<double>       &exp_known,
                                      const std::vector<double>       &exp_unknown,
                                      const std::vector<double>       &std_known,
                                      const std::vector<double>       &std_unknown,
                                      const Variogram                 &vario,
                                      Vector                          &exp,
                                      SymmetricMatrix                 &cov);

  void EliminateLargeCorr(std::vector<NRLib::Point> & pos_known,
                          Vector                    & delta_known,
                          std::vector<double>       & std_known,
                          const std::vector<double> & local_corr,
                          const Variogram           & vario);

  void IdentifyLargeCorr(const std::vector<NRLib::Point> & pos_known,
                         const std::vector<double>       & local_corr,
                         const Variogram                 & vario,
                         std::vector<int>                & obs_chain_map,
                         std::vector<std::vector<int> >  & obs_chains);


  void Posterior1DNormal(double exp_prior, double cov_prior,
                       double exp_obs, double cov_obs,
                       double corr, double obs,
                       double & exp_post, double & cov_post);



  bool CondDistrMultiNormalCore(const SymmetricMatrix & sigma_11,
                              const Matrix          & sigma_12,
                              SymmetricMatrix       & sigma_22,
                              const Vector          & mu_unknown,
                              const Vector          & delta_known,
                              Vector                & exp,
                              SymmetricMatrix       & cov);

}

#endif
