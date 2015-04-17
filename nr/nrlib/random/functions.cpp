// $Id: functions.cpp 882 2011-09-23 13:10:16Z perroe $

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

#include "functions.hpp"
#include "../exception/exception.hpp"

#include <cmath>

// The implementation of the functions in this file is based on the boost
// implemenation for 53-bits presision, corresponding to 64-bit double.

namespace NRLib {

// Calculation of error function and inverse error function is obtained from boost.
//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

double erf_imp(double z, bool invert)
{
  if(z < 0)
  {
    if(!invert)
      return -erf_imp(-z, invert);
    else if(z < -0.5)
      return 2 - erf_imp(-z, invert);
    else
      return 1 + erf_imp(-z, false);
  }

  double result;

  //
  // Big bunch of selection statements now to pick
  // which implementation to use,
  // try to put most likely options first:
  //
  if(z < 0.5)
  {
    //
    // We're going to calculate erf:
    //
    if(z == 0)
    {
      result = 0.0;
    }
    else if(z < 1e-10)
    {
      result = static_cast<double>(z * 1.125f + z * 0.003379167095512573896158903121545171688L);
    }
    else
    {
      // Maximum Deviation Found:                     1.561e-17
      // Expected Error Term:                         1.561e-17
      // Maximum Relative Change in Control Points:   1.155e-04
      // Max Error found at double precision =        2.961182e-17

      static const double Y = 1.044948577880859375;
      static const double P[] = {
        0.0834305892146531832907L,
        -0.338165134459360935041L,
        -0.0509990735146777432841L,
        -0.00772758345802133288487L,
        -0.000322780120964605683831L,
      };
      static const double Q[] = {
        1L,
        0.455004033050794024546L,
        0.0875222600142252549554L,
        0.00858571925074406212772L,
        0.000370900071787748000569L,
      };
      result = z * (Y + evaluate_polynomial(P, z * z) / evaluate_polynomial(Q, z * z));
    }
  }
  else if((z < 14) || ((z < 28) && invert))
  {
    //
    // We'll be calculating erfc:
    //
    invert = !invert;
    if(z < 1.5f)
    {
      // Maximum Deviation Found:                     3.702e-17
      // Expected Error Term:                         3.702e-17
      // Maximum Relative Change in Control Points:   2.845e-04
      // Max Error found at double precision =        4.841816e-17
      static const double Y = 0.405935764312744140625;
      static const double P[] = {
        -0.098090592216281240205L,
        0.178114665841120341155L,
        0.191003695796775433986L,
        0.0888900368967884466578L,
        0.0195049001251218801359L,
        0.00180424538297014223957L,
      };
      static const double Q[] = {
        1L,
        1.84759070983002217845L,
        1.42628004845511324508L,
        0.578052804889902404909L,
        0.12385097467900864233L,
        0.0113385233577001411017L,
        0.337511472483094676155e-5L,
      };
      result = Y + evaluate_polynomial(P, z - 0.5) / evaluate_polynomial(Q, z - 0.5);
      result *= std::exp(-z * z) / z;
    }
    else if(z < 2.5f)
    {
      // Max Error found at double precision =        6.599585e-18
      // Maximum Deviation Found:                     3.909e-18
      // Expected Error Term:                         3.909e-18
      // Maximum Relative Change in Control Points:   9.886e-05
      static const double Y = 0.50672817230224609375;
      static const double P[] = {
        -0.0243500476207698441272L,
        0.0386540375035707201728L,
        0.04394818964209516296L,
        0.0175679436311802092299L,
        0.00323962406290842133584L,
        0.000235839115596880717416L,
      };
      static const double Q[] = {
        1L,
        1.53991494948552447182L,
        0.982403709157920235114L,
        0.325732924782444448493L,
        0.0563921837420478160373L,
        0.00410369723978904575884L,
      };
      result = Y + evaluate_polynomial(P, z - 1.5) / evaluate_polynomial(Q, z - 1.5);
      result *= std::exp(-z * z) / z;
    }
    else if(z < 4.5f)
    {
      // Maximum Deviation Found:                     1.512e-17
      // Expected Error Term:                         1.512e-17
      // Maximum Relative Change in Control Points:   2.222e-04
      // Max Error found at double precision =        2.062515e-17
      static const double Y = 0.5405750274658203125;
      static const double P[] = {
        0.00295276716530971662634L,
        0.0137384425896355332126L,
        0.00840807615555585383007L,
        0.00212825620914618649141L,
        0.000250269961544794627958L,
        0.113212406648847561139e-4L,
      };
      static const double Q[] = {
        1L,
        1.04217814166938418171L,
        0.442597659481563127003L,
        0.0958492726301061423444L,
        0.0105982906484876531489L,
        0.000479411269521714493907L,
      };
      result = Y + evaluate_polynomial(P, z - 3.5) / evaluate_polynomial(Q, z - 3.5);
      result *= std::exp(-z * z) / z;
    }
    else
    {
      // Max Error found at double precision =        2.997958e-17
      // Maximum Deviation Found:                     2.860e-17
      // Expected Error Term:                         2.859e-17
      // Maximum Relative Change in Control Points:   1.357e-05
      static const double Y = 0.5579090118408203125;
      static const double P[] = {
        0.00628057170626964891937L,
        0.0175389834052493308818L,
        -0.212652252872804219852L,
        -0.687717681153649930619L,
        -2.5518551727311523996L,
        -3.22729451764143718517L,
        -2.8175401114513378771L,
      };
      static const double Q[] = {
        1L,
        2.79257750980575282228L,
        11.0567237927800161565L,
        15.930646027911794143L,
        22.9367376522880577224L,
        13.5064170191802889145L,
        5.48409182238641741584L,
      };
      result = Y + evaluate_polynomial(P, 1 / z) / evaluate_polynomial(Q, 1 / z);
      result *= std::exp(-z * z) / z;
    }
  }
  else
  {
    //
    // Any value of z larger than 28 will underflow to zero:
    //
    result = 0;
    invert = !invert;
  }

  if(invert)
  {
    result = 1 - result;
  }

  return result;
}


//
// The inverse erf and erfc functions share a common implementation,
// this version is for 80-bit long double's and smaller:
//
double erf_inv_imp(double p, double q)
{
  double result = 0;

  if(p <= 0.5)
  {
    //
    // Evaluate inverse erf using the rational approximation:
    //
    // x = p(p+10)(Y+R(p))
    //
    // Where Y is a constant, and R(p) is optimised for a low
    // absolute error compared to |Y|.
    //
    // double: Max error found: 2.001849e-18
    // long double: Max error found: 1.017064e-20
    // Maximum Deviation Found (actual error term at infinite precision) 8.030e-21
    //
    static const float Y = 0.0891314744949340820313f;
    static const double P[] = {
      -0.000508781949658280665617L,
      -0.00836874819741736770379L,
      0.0334806625409744615033L,
      -0.0126926147662974029034L,
      -0.0365637971411762664006L,
      0.0219878681111168899165L,
      0.00822687874676915743155L,
      -0.00538772965071242932965L
    };
    static const double Q[] = {
      1,
      -0.970005043303290640362L,
      -1.56574558234175846809L,
      1.56221558398423026363L,
      0.662328840472002992063L,
      -0.71228902341542847553L,
      -0.0527396382340099713954L,
      0.0795283687341571680018L,
      -0.00233393759374190016776L,
      0.000886216390456424707504L
    };
    double g = p * (p + 10);
    double r = evaluate_polynomial(P, p) / evaluate_polynomial(Q, p);
    result = g * Y + g * r;
  }
  else if(q >= 0.25)
  {
    //
    // Rational approximation for 0.5 > q >= 0.25
    //
    // x = sqrt(-2*log(q)) / (Y + R(q))
    //
    // Where Y is a constant, and R(q) is optimised for a low
    // absolute error compared to Y.
    //
    // double : Max error found: 7.403372e-17
    // long double : Max error found: 6.084616e-20
    // Maximum Deviation Found (error term) 4.811e-20
    //
    static const float Y = 2.249481201171875f;
    static const double P[] = {
      -0.202433508355938759655L,
      0.105264680699391713268L,
      8.37050328343119927838L,
      17.6447298408374015486L,
      -18.8510648058714251895L,
      -44.6382324441786960818L,
      17.445385985570866523L,
      21.1294655448340526258L,
      -3.67192254707729348546L
    };
    static const double Q[] = {
      1L,
      6.24264124854247537712L,
      3.9713437953343869095L,
      -28.6608180499800029974L,
      -20.1432634680485188801L,
      48.5609213108739935468L,
      10.8268667355460159008L,
      -22.6436933413139721736L,
      1.72114765761200282724L
    };
    double g = sqrt(-2 * log(q));
    double xs = q - 0.25;
    double r = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
    result = g / (Y + r);
  }
  else
  {
    //
    // For q < 0.25 we have a series of rational approximations all
    // of the general form:
    //
    // let: x = sqrt(-log(q))
    //
    // Then the result is given by:
    //
    // x(Y+R(x-B))
    //
    // where Y is a constant, B is the lowest value of x for which
    // the approximation is valid, and R(x-B) is optimised for a low
    // absolute error compared to Y.
    //
    // Note that almost all code will really go through the first
    // or maybe second approximation.  After than we're dealing with very
    // small input values indeed: 80 and 128 bit long double's go all the
    // way down to ~ 1e-5000 so the "tail" is rather long...
    //
    double x = sqrt(-log(q));
    if(x < 3)
    {
      // Max error found: 1.089051e-20
      static const float Y = 0.807220458984375f;
      static const double P[] = {
        -0.131102781679951906451L,
        -0.163794047193317060787L,
        0.117030156341995252019L,
        0.387079738972604337464L,
        0.337785538912035898924L,
        0.142869534408157156766L,
        0.0290157910005329060432L,
        0.00214558995388805277169L,
        -0.679465575181126350155e-6L,
        0.285225331782217055858e-7L,
        -0.681149956853776992068e-9L
      };
      static const double Q[] = {
        1,
        3.46625407242567245975L,
        5.38168345707006855425L,
        4.77846592945843778382L,
        2.59301921623620271374L,
        0.848854343457902036425L,
        0.152264338295331783612L,
        0.01105924229346489121L
      };
      double xs = x - 1.125;
      double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
      result = Y * x + R * x;
    }
    else if(x < 6)
    {
      // Max error found: 8.389174e-21
      static const float Y = 0.93995571136474609375f;
      static const double P[] = {
        -0.0350353787183177984712L,
        -0.00222426529213447927281L,
        0.0185573306514231072324L,
        0.00950804701325919603619L,
        0.00187123492819559223345L,
        0.000157544617424960554631L,
        0.460469890584317994083e-5L,
        -0.230404776911882601748e-9L,
        0.266339227425782031962e-11L
      };
      static const double Q[] = {
        1L,
        1.3653349817554063097L,
        0.762059164553623404043L,
        0.220091105764131249824L,
        0.0341589143670947727934L,
        0.00263861676657015992959L,
        0.764675292302794483503e-4L
      };
      double xs = x - 3;
      double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
      result = Y * x + R * x;
    }
    else if(x < 18)
    {
      // Max error found: 1.481312e-19
      static const float Y = 0.98362827301025390625f;
      static const double P[] = {
        -0.0167431005076633737133L,
        -0.00112951438745580278863L,
        0.00105628862152492910091L,
        0.000209386317487588078668L,
        0.149624783758342370182e-4L,
        0.449696789927706453732e-6L,
        0.462596163522878599135e-8L,
        -0.281128735628831791805e-13L,
        0.99055709973310326855e-16L
      };
      static const double Q[] = {
        1L,
        0.591429344886417493481L,
        0.138151865749083321638L,
        0.0160746087093676504695L,
        0.000964011807005165528527L,
        0.275335474764726041141e-4L,
        0.282243172016108031869e-6L
      };
      double xs = x - 6;
      double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
      result = Y * x + R * x;
    }
    else if(x < 44)
    {
      // Max error found: 5.697761e-20
      static const float Y = 0.99714565277099609375f;
      static const double P[] = {
        -0.0024978212791898131227L,
        -0.779190719229053954292e-5L,
        0.254723037413027451751e-4L,
        0.162397777342510920873e-5L,
        0.396341011304801168516e-7L,
        0.411632831190944208473e-9L,
        0.145596286718675035587e-11L,
        -0.116765012397184275695e-17L
      };
      static const double Q[] = {
        1L,
        0.207123112214422517181L,
        0.0169410838120975906478L,
        0.000690538265622684595676L,
        0.145007359818232637924e-4L,
        0.144437756628144157666e-6L,
        0.509761276599778486139e-9L
      };
      double xs = x - 18;
      double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
      result = Y * x + R * x;
    }
    else
    {
      // Max error found: 1.279746e-20
      static const float Y = 0.99941349029541015625f;
      static const double P[] = {
        -0.000539042911019078575891L,
        -0.28398759004727721098e-6L,
        0.899465114892291446442e-6L,
        0.229345859265920864296e-7L,
        0.225561444863500149219e-9L,
        0.947846627503022684216e-12L,
        0.135880130108924861008e-14L,
        -0.348890393399948882918e-21L
      };
      static const double Q[] = {
        1L,
        0.0845746234001899436914L,
        0.00282092984726264681981L,
        0.468292921940894236786e-4L,
        0.399968812193862100054e-6L,
        0.161809290887904476097e-8L,
        0.231558608310259605225e-11L
      };
      double xs = x - 44;
      double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
      result = Y * x + R * x;
    }
  }
  return result;
}


double erf(double z)
{
  return erf_imp(z, false);
}


double erfc(double z)
{
  return erf_imp(z, true);
}


double erf_inv(double z)
{
  if ((z <= -1) || (z >= 1))
    throw Exception("Parameter error in inverse error function.");
  if (z == 0)
    return 0;
  //
  // Normalise the input, so it's in the range [0,1], we will
  // negate the result if z is outside that range.  This is a simple
  // application of the erf reflection formula: erf(-z) = -erf(z)
  //
  double p, q, s;
  if (z < 0) {
    p = -z;
    q = 1 - p;
    s = -1;
  }
  else {
    p = z;
    q = 1 - z;
    s = 1;
  }
  return s * erf_inv_imp(p, q);
}


double erfc_inv(double z)
{
  if((z <= 0) || (z >= 2))
    throw Exception("Parameter error in inverse complementary error function.");

  //
  // Normalise the input, so it's in the range [0,1], we will
  // negate the result if z is outside that range.  This is a simple
  // application of the erfc reflection formula: erfc(-z) = 2 - erfc(z)
  //
  double p, q, s;
  if(z > 1) {
    q = 2 - z;
    p = 1 - q;
    s = -1;
  }
  else {
    p = 1 - z;
    q = z;
    s = 1;
  }
   return s * erf_inv_imp(p, q);
}

} // namespace NRLib
