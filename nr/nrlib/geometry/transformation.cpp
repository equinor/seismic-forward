// $Id: transformation.cpp 969 2012-02-07 15:03:56Z perroe $

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

#include "transformation.hpp"
#include "point.hpp"

#include <cmath>

using namespace NRLib;

Transformation::Transformation()
{
  Identity(transformation_);
}


Transformation& Transformation::TranslateX(double translation)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][3] = translation;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::TranslateY(double translation)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[1][3] = translation;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::TranslateZ(double translation)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[2][3] = translation;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::StretchX(double expansion_factor)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][0] = expansion_factor;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::StretchY(double expansion_factor)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[1][1] = expansion_factor;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::StretchZ(double expansion_factor)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[2][2] = expansion_factor;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::RotateX(double angle)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[1][1] = std::cos(angle);
  deform[1][2] = -std::sin(angle);
  deform[2][1] = std::sin(angle);
  deform[2][2] = std::cos(angle);
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::RotateY(double angle)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][0] = std::cos(angle);
  deform[0][2] = std::sin(angle);
  deform[2][0] = -std::sin(angle);
  deform[2][2] = std::cos(angle);
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::RotateZ(double angle)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][0] = std::cos(angle);
  deform[0][1] = -std::sin(angle);
  deform[1][0] = std::sin(angle);
  deform[1][1] = std::cos(angle);
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::ShearXY(double shear_x, double shear_y)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][2] = shear_x;
  deform[1][2] = shear_y;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::ShearXZ(double shear_x, double shear_z)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][1] = shear_x;
  deform[2][1] = shear_z;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::ShearYZ(double shear_y, double shear_z)
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[1][0] = shear_y;
  deform[2][0] = shear_z;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::MirrorXY()
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[2][2] = -1;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::MirrorXZ()
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[1][1] = -1;
  AddTransformation(deform);
  return *this;
}


Transformation& Transformation::MirrorYZ()
{
  Transformation::Matrix deform;
  Identity(deform);
  deform[0][0] = -1;
  AddTransformation(deform);
  return *this;
}


void Transformation::Transform(NRLib::Point& pt) const
{
  double xp[3];

  for (int i = 0; i < 3; ++i) {
    xp[i] =
      transformation_[i][0] * pt.x +
      transformation_[i][1] * pt.y +
      transformation_[i][2] * pt.z +
      transformation_[i][3];
  }

  pt.x = xp[0]; pt.y = xp[1]; pt.z = xp[2];
}


void Transformation::Identity(Transformation::Matrix& matrix)
{
  unsigned int i;
  for (i = 0; i < 4; i++) {
    unsigned int j;
    for (j = 0; j < 4; j++) {
      matrix[i][j] = (i == j ? 1 : 0 );
    }
  }
}


void Transformation::AddTransformation(const Transformation::Matrix& deform)
{
  Matrix tmp_transformation;

  unsigned int I;
  for (I = 0; I < 4; I++) {
    unsigned int J;
    for (J = 0; J < 4; J++) {
      tmp_transformation[I][J] = transformation_[I][J];
    }
  }

  for (I = 0; I < 4; I++) {
    unsigned int J;
    for (J = 0; J < 4; J++) {
      transformation_[I][J] = deform[I][0] * tmp_transformation[0][J];
      unsigned int K;
      for ( K = 1; K < 4; K++) {
        transformation_[I][J] +=
          deform[I][K] * tmp_transformation[K][J];
      }
    }
  }
}
