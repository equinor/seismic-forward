// $Id: posteriormultinormal.cpp 1320 2016-03-06 13:01:47Z perroe $

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

#include <vector>

#include "posteriormultinormal.hpp"
#include "../math/constants.hpp"
#include "../geometry/point.hpp"


using namespace NRLib;

//
// Calculates conditional mean (exp) and covariance matrix (cov)
//
//

bool NRLib::CondDistrMultiNormal(const std::vector<double> & x_known,
                                 const std::vector<double> & known_values,
                                 const std::vector<double> & x_unknown,
                                 const std::vector<double> & exp_known,
                                 const std::vector<double> & exp_unknown,
                                 const std::vector<double> & std_known,
                                 const std::vector<double> & std_unknown,
                                 const Variogram           & vario,
                                 Vector                    & exp,
                                 SymmetricMatrix           & cov)
{

  int n_known = static_cast<int>(x_known.size());
  int n_unknown = static_cast<int>(x_unknown.size());
  assert(static_cast<int>(known_values.size()) == n_known);
  assert(static_cast<int>(exp_known.size()) == n_known);
  assert(static_cast<int>(exp_unknown.size()) == n_unknown);

  exp = Vector(n_unknown);
  cov = SymmetricMatrix(n_unknown);
  if(n_known == 0) { // no observations, return unconditional mean and covariance
    for(int i = 0; i < n_unknown; i++) {
      exp(i) = exp_unknown[i];
      for(int j = 0; j < n_unknown; j++) {
        cov(i, j) = std_unknown[i]*std_unknown[j]*vario.GetCorr(x_unknown[i]-x_unknown[j]);
      }
    }
    return true;
  }

  NRLib::Vector delta_known(n_known);
  NRLib::Vector mu_unknown(n_unknown);

  for(int i = 0; i < n_known; i++) {
    delta_known(i) =  known_values[i] - exp_known[i];
  }
  for(int i = 0; i < n_unknown; i++)
    mu_unknown(i) = exp_unknown[i];

  SymmetricMatrix sigma_11, sigma_22;
  Matrix sigma_12;
  sigma_11 = SymmetricMatrix(n_unknown);
  sigma_22 = SymmetricMatrix(n_known);
  sigma_12 = Matrix(n_unknown, n_known);
  Matrix sigma_21 = Matrix(n_known, n_unknown);
  for(int i = 0; i < n_unknown; i++) {
    for(int j = 0; j <= i; j++) {
      sigma_11(j,i) = std_unknown[i]*std_unknown[j]*vario.GetCorr(x_unknown[i]-x_unknown[j]);
    }
  }

  for(int i = 0; i < n_known; i++) {
    for(int j = 0; j <= i; j++)
      sigma_22(j,i) = std_known[i]*std_known[j]*vario.GetCorr(x_known[i]-x_known[j]);

    for(int j = 0; j < n_unknown; j++)
      sigma_12(j,i) = std_known[i]*std_unknown[j]*vario.GetCorr(x_known[i]-x_unknown[j]);
  }

  CondDistrMultiNormalCore(sigma_11,
                           sigma_12,
                           sigma_22,
                           mu_unknown,
                           delta_known,
                           exp,
                           cov);

  return true;
}




bool NRLib::SecondayCondDistrMultiNormal2D(const std::vector<NRLib::Point> & pos_known,
                                           const std::vector<double>       & known_values,
                                           const std::vector<double>       & local_corr,
                                           const std::vector<NRLib::Point> & pos_unknown,
                                           const std::vector<double>       & exp_known,
                                           const std::vector<double>       & exp_unknown,
                                           const std::vector<double>       & std_known,
                                           const std::vector<double>       & std_unknown,
                                           const Variogram                 & vario,
                                           Vector                          & exp,
                                           SymmetricMatrix                 & cov)
{

  int n_known = static_cast<int>(pos_known.size());
  int n_unknown = static_cast<int>(pos_unknown.size());
  assert(static_cast<int>(known_values.size()) == n_known);
  assert(static_cast<int>(exp_known.size()) == n_known);
  assert(static_cast<int>(exp_unknown.size()) == n_unknown);

  exp = Vector(n_unknown);
  cov = SymmetricMatrix(n_unknown);
  if(n_known == 0) { // no observations, return unconditional mean and covariance
    for(int i = 0; i < n_unknown; i++) {
      exp(i) = exp_unknown[i];
      for(int j = 0; j < n_unknown; j++) {
        cov(i, j) = std_unknown[i]*std_unknown[j]*vario.GetCorr(pos_unknown[i].x-pos_unknown[j].x,
                                                                pos_unknown[i].y-pos_unknown[j].y);
      }
    }
    return true;
  }


  NRLib::Vector delta_known(n_known);
  NRLib::Vector mu_unknown(n_unknown);

  for(int i = 0; i < n_known; i++) {
    delta_known(i) =  known_values[i] - exp_known[i];
  }
  for(int i = 0; i < n_unknown; i++)
    mu_unknown(i) = exp_unknown[i];

  SymmetricMatrix sigma_11, sigma_22;
  Matrix sigma_12;
  sigma_11 = SymmetricMatrix(n_unknown);
  Matrix sigma_21 = Matrix(n_known, n_unknown);
  for(int i = 0; i < n_unknown; i++) {
    for(int j = 0; j <= i; j++) {
      sigma_11(j,i) = std_unknown[i]*std_unknown[j]*vario.GetCorr(pos_unknown[i].x-pos_unknown[j].x,
                                                                  pos_unknown[i].y-pos_unknown[j].y);
    }
  }

  std::vector<NRLib::Point> pos_known_used = pos_known;
  std::vector<double>       std_known_used = std_known;
  EliminateLargeCorr(pos_known_used, delta_known, std_known_used, local_corr, vario);
  int n_used = static_cast<int>(pos_known_used.size());

  sigma_22 = SymmetricMatrix(n_used);
  sigma_12 = Matrix(n_unknown, n_used);

  for(int i = 0; i < n_used; i++) {
    for(int j = 0; j <= i; j++) {
      double corr;
      if(local_corr[i] > local_corr[j])
        corr = local_corr[j];
      else if(local_corr[j] > local_corr[i])
        corr = local_corr[i];
      else
        corr = 1.0;

      corr *= vario.GetCorr(pos_known_used[i].x-pos_known_used[j].x, pos_known_used[i].y-pos_known_used[j].y);

      sigma_22(j,i) = std_known_used[i]*std_known_used[j]*corr;
    }

    for(int j = 0; j < n_unknown; j++)
      sigma_12(j,i) = std_known_used[i]*std_unknown[j]*local_corr[i]*vario.GetCorr(pos_known_used[i].x-pos_unknown[j].x,
                                                                              pos_known_used[i].y-pos_unknown[j].y);
  }



  CondDistrMultiNormalCore(sigma_11,
                           sigma_12,
                           sigma_22,
                           mu_unknown,
                           delta_known,
                           exp,
                           cov);

  return true;
}


bool NRLib::CondDistrMultiNormalCore(const SymmetricMatrix & sigma_11,
                                     const Matrix          & sigma_12,
                                     SymmetricMatrix       & sigma_22,
                                     const Vector          & mu_unknown,
                                     const Vector          & delta_known,
                                     Vector                & exp,
                                     SymmetricMatrix       & cov)
{
  CholeskyInvert(sigma_22);

  Matrix sigma_21 = NRLib::transpose(sigma_12);

  NRLib::Vector sigma_22_delta_known = sigma_22*delta_known;
  Matrix sigma_22_sigma_21 = sigma_22*sigma_21;
  exp = mu_unknown + sigma_12*sigma_22_delta_known;
  Matrix help = sigma_12*sigma_22_sigma_21;
  SymmetricMatrix symmetric_help = help.upper();
  cov = sigma_11 - symmetric_help;

  return true;

}

void
NRLib::EliminateLargeCorr(std::vector<NRLib::Point> & pos_known,
                          Vector                    & delta_known,
                          std::vector<double>       & std_known,
                          const std::vector<double> & local_corr,
                          const Variogram           & vario)
{
  std::vector<int>                obs_chain_map;
  std::vector<std::vector<int> >  obs_chains;
  IdentifyLargeCorr(pos_known, local_corr, vario, obs_chain_map, obs_chains);

  std::vector<NRLib::Point> pos_known_used(obs_chains.size());
  Vector                    delta_known_used(static_cast<int>(obs_chains.size()));
  std::vector<double>       std_known_used(obs_chains.size());

  for(size_t i=0;i<obs_chains.size();i++) {
    //For simplicity, use values from first obs, instead of averaging.
    pos_known_used[i]   = pos_known[obs_chains[i][0]];
    delta_known_used(static_cast<int>(i)) = delta_known(static_cast<int>(obs_chains[i][0]));
    std_known_used[i]   = std_known[obs_chains[i][0]];
  }
  pos_known   = pos_known_used;
  delta_known = delta_known_used;
  std_known   = std_known_used;
}

void
NRLib::IdentifyLargeCorr(const std::vector<NRLib::Point> & pos_known,
                         const std::vector<double>       & local_corr,
                         const Variogram                 & vario,
                         std::vector<int>                & obs_chain_map,
                         std::vector<std::vector<int> >  & obs_chains)
{
  int n_known = static_cast<int>(pos_known.size());
  obs_chain_map.resize(n_known);
  for(int i=0;i<n_known;i++)
    obs_chain_map[i] = -1;

  int chain_no = 0;
  for(int base = 0; base < n_known; base++) {
    std::vector<int> chain;
    if(obs_chain_map[base] < 0) {
      chain.push_back(base);
      obs_chain_map[base] = chain_no;
      for(size_t cur = 0; cur < chain.size(); cur++) {
        int i = chain[cur];
        for(int j=i+1;j<n_known;j++) {
          if(obs_chain_map[j] < 0) {
            double corr;
            if(local_corr[i] > local_corr[j])
              corr = local_corr[j];
            else if(local_corr[j] > local_corr[i])
              corr = local_corr[j];
            else
              corr = 1.0;
            corr *= vario.GetCorr(pos_known[i].x-pos_known[j].x, pos_known[i].y-pos_known[j].y);
            if(corr > 0.999) {
              chain.push_back(j);
              obs_chain_map[j] = chain_no;
            }
          }
        }
      }
      obs_chains.push_back(chain);
      chain_no++;
    }
  }
}



// Calculates posterior distribution in binormal distribution when one of the variables is observed

void NRLib::Posterior1DNormal(double   exp_prior,
                              double   cov_prior,
                              double   exp_obs,
                              double   cov_obs,
                              double   corr,
                              double   obs,
                              double & exp_post,
                              double & cov_post)
{
  exp_post = exp_prior + (corr*sqrt(cov_prior/cov_obs))*(obs - exp_obs);
  cov_post = (1-corr*corr)*cov_prior;

  return;
}
