// $Id: multinormal.cpp 1242 2014-02-18 12:40:58Z anner $

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


#include "multinormal.hpp"
#include "../random/normal.hpp"

namespace NRLib {

MultiNormal::MultiNormal()
{

}

MultiNormal::MultiNormal(const Vector    & exp,
                         SymmetricMatrix   cov)
                         : exp_(exp), cov_(cov)
{

  CholeskyFactorize(cov);
  chol_fac_cov_ = cov;

}


MultiNormal::~MultiNormal()
{}



Vector
MultiNormal::Draw() const
{
  Vector result = exp_;
  Vector rand_01(static_cast<int>(exp_.length()));
  Normal normal;
  for (int i = 0; i < rand_01.length(); i++)
    rand_01(i) = normal.Draw();

  Vector help = chol_fac_cov_*rand_01;
  result +=help;

  return result;

}


Vector
MultiNormal::Draw(RandomGenerator & rg) const
{
  Vector result = exp_;
  Vector rand_01(static_cast<int>(exp_.length()));
  Normal normal;
  for (int i = 0; i < rand_01.length(); i++)
    rand_01(i) = normal.Draw(rg);

  Vector help = chol_fac_cov_*rand_01;
  result += help;

  return result;

}


double
MultiNormal::Potential(const Vector & x)  const
{
  Matrix e_mat;
  Vector e_vec;
  ComputeEigenVectorsSymmetric(cov_, e_vec, e_mat);

  double det(1.0);
  int dim = cov_.dim();
  for (int i = 0; i < dim; i++)
    det *= e_vec(i);
  assert(det > 0);
  Vector diff = x - exp_;
  Vector e_vec_inv = e_vec;
  for(int i = 0; i < dim; i++)
    e_vec_inv(i) = 1.0/e_vec(i);
  double sum(0.0);
  for (int i = 0; i < dim; i++) {
    //double r1 = x(i) - mu(i);
    double r1 = diff(i);
    for (int j = 0; j < dim; j++) {
      double q1 = e_mat(i,j);
      //double s = 1.0 / e_vec(j);
      double r1q1s = r1*q1*e_vec_inv(j);
      for (int k = 0; k < dim; k++) {
        //double r2 = x(k) - mu(k);
        double r2 = diff(k);
        double q2 = e_mat(k,j);
        //sum += q1 * q2 * r1 * r2 * s;
        sum += r1q1s*r2*q2;
      }
    }
  }

  double potential = 0.5 * (dim * NRLib::Log2Pi + std::log(det) + sum);

  return(potential);
}








}



