// $Id: transformation.hpp 969 2012-02-07 15:03:56Z perroe $

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

#ifndef NRLIB_GEOMETRY_TRANSFORMATION_HPP
#define NRLIB_GEOMETRY_TRANSFORMATION_HPP

namespace NRLib {

class Point;

/// Defines a transformation transforming a point.
/// Possible transformations include translation, rotation,
/// stretching, shearing, and mirroring.
class Transformation {
public:
  /// Default constructor. Initialises to the identity transformation.
  Transformation();


  /// Add a translation in z direction.
  /// \param translation
  Transformation& TranslateX(double translation);

  /// Add a translation in z direction.
  /// \param translation
  Transformation& TranslateY(double translation);

  /// Add a translation in z direction.
  /// \param translation
  Transformation& TranslateZ(double translation);

  /// Add a stretch in x direction.
  /// \param expansion_factor
  Transformation& StretchX(double expansion_factor);

  /// Add a stretch in y direction.
  /// \param expansion_factor
  Transformation& StretchY(double expansion_factor);

  /// Add a stretch in z direction.
  /// \param expansion_factor
  Transformation& StretchZ(double expansion_factor);

  /// Add a global rigid body rotation around the x axis.
  /// \param angle rotation angle in radians.
  Transformation& RotateX(double angle);

  /// Add a global rigid body rotation around the y axis.
  /// \param angle rotation angle in radians.
  Transformation& RotateY(double angle);

  /// Add a global rigid body rotation around the z axis.
  /// \param angle rotation angle in radians.
  Transformation& RotateZ(double angle);

  /// Add a deformation specified by the off-diagonal
  /// components of the shear tensor. This is a non-rigid body operation.
  /// \param shear_x shear component for x.
  /// \param shear_y shear component for y.
  Transformation& ShearXY(double shear_x, double shear_y);

  /// Add a deformation specified by the off-diagonal
  /// components of the shear tensor. This is a non-rigid body operation.
  /// \param shear_x shear component for x.
  /// \param shear_z shear component for z.
  Transformation& ShearXZ(double shear_x, double shear_z);

  /// Add a deformation specified by the off-diagonal
  /// components of the shear tensor. This is a non-rigid body operation.
  /// \param shear_y shear component for y.
  /// \param shear_z shear component for z.
  Transformation& ShearYZ(double shear_y, double shear_z);

  /// Mirror across the XY-plane.
  Transformation& MirrorXY();

  /// Mirror across the XZ-plane.
  Transformation& MirrorXZ();

  /// Mirror across the YZ-plane.
  Transformation& MirrorYZ();

  /// Apply transformation to a point.
  /// The point is moved according to the transformation matrix.
  /// \param pt  Point to move.
  void Transform(NRLib::Point& pt) const;

private:
  typedef double Matrix[4][4];

  /// Transformation matrix.
  Matrix transformation_;

  /// Initializes the input matrix to the identity matrix.
  /// \param matrix matrix to initialize.
  static void Identity(Transformation::Matrix& matrix);

  /// Add transformation to the total transformation.
  /// \param deform transformation to add.
  void AddTransformation(const Transformation::Matrix& deform);
};

}

#endif
