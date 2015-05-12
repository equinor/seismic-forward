// $Id: nrlib_flens.hpp 1282 2014-09-01 12:27:47Z enesvold $

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

#ifndef NRLIB_FLENS_HPP
#define NRLIB_FLENS_HPP

#include <complex>
#include <flens/flens.h>
#include "../exception/exception.hpp"

namespace NRLib {
  using flens::conjugateTranspose;
  using flens::transpose;
  using flens::_;

  typedef flens::GeMatrix<flens::FullStorage<double, flens::ColMajor> >               Matrix;
  typedef flens::GeMatrix<flens::FullStorage<float, flens::ColMajor> >                FloatMatrix;
  typedef flens::GeMatrix<flens::FullStorage<std::complex<double>, flens::ColMajor> > ComplexMatrix;

  typedef flens::DenseVector<flens::Array<double> >                Vector;
  typedef flens::DenseVector<flens::Array<float> >                 FloatVector;
  typedef flens::DenseVector<flens::Array<int> >                   IntegerVector;
  typedef flens::DenseVector<flens::Array<std::complex<double> > > ComplexVector;

  typedef flens::SyMatrix<flens::FullStorage<double, flens::ColMajor> > SymmetricMatrix;
  typedef flens::SyMatrix<flens::FullStorage<float,  flens::ColMajor> > SymmetricFloatMatrix;

  // Use Invert instead!!!
  template <typename FS> void invert(flens::GeMatrix<FS> &A) { return Invert(A); }
  template <typename FS> void Invert(flens::GeMatrix<FS> &A);

  template <typename FS> void TriangleFactorize(flens::GeMatrix<FS> & A,
                                                flens::DenseVector<flens::Array<int> > & p);

  template <typename FS> void TriangleInvert(flens::GeMatrix<FS> & A,
                                             flens::DenseVector<flens::Array<int> > & p);

  template <typename T>  void SetMatrixFrom2DArray(Matrix &  A,
                                                   T      ** F);

  template <typename T>  void Set2DArrayFromMatrix(const Matrix &  A,
                                                   T            ** F);

  void            InitializeMatrix(Matrix           & A,
                                   double             value);

  void            InitializeComplexMatrix(ComplexMatrix         & A,
                                          std::complex<double>    value);

  void            InitializeSymmetricMatrix(SymmetricMatrix     & A,
                                            double                value);

  double          FindLargestElement(const NRLib::Vector & v);

  Vector          ZeroVector(int n);
  Matrix          ZeroMatrix(int n);
  SymmetricMatrix SymmetricZeroMatrix(int n);
  Matrix          IdentityMatrix(int n);

  void            Sort3x3(Vector & Eval,
                          Matrix & Evec);

  /// \brief Solves the equation system Ax = b, where A is symmetric
  ///        and positive definite using Cholesky factorization.
  /// \throw Exception if A is not positive definite.
  void CholeskySolve(const SymmetricMatrix & A,
                     const Vector          & b,
                     Vector                & x);

  void CholeskySolve(const SymmetricMatrix & A,
                     Matrix                & B); // This is also where the solution is stored

  void CholeskySolveComplex(ComplexMatrix & A,
                            ComplexMatrix & B); // This is also where the solution is stored

  void CholeskyInvert(SymmetricMatrix & A);


  void CholeskyFactorize(SymmetricMatrix & A);


  void ComputeEigenVectors(Matrix & A,
                           Vector & eigen_values,
                           Matrix & eigen_vectors);

  void ComputeEigenVectorsComplex(ComplexMatrix        & A,
                                  ComplexVector        & eigen_values,
                                  ComplexMatrix        & eigen_vectors);

  void ComputeEigenVectorsSymmetric(const SymmetricMatrix & A,
                                    Vector                & eigen_values,
                                    Matrix                & eigen_vectors);

  /// \brief File format for input and output of marixes and vectors.
  enum LinalgFileFormat {
    MatrixAscii,    ///< Space-separated columns with one row for each line.
    MatrixBinary    ///< Binary, column major.
  };

  /// \brief Read matrix from file.
  Matrix ReadMatrixFromFile(const std::string & filename,
                            int                 nRows,
                            int                 nCols,
                            LinalgFileFormat    format = MatrixAscii);

  Vector ReadVectorFromFile(const std::string & filename,
                            int                 n,
                            LinalgFileFormat    format = MatrixAscii);

  void   Adjoint(const ComplexMatrix & in,
                 ComplexMatrix       & out);

  void   WriteComplexVector(const std::string   & header,
                            const ComplexVector & c);

  void   WriteComplexMatrix(const std::string   & header,
                            const ComplexMatrix & C);

  void   WriteMatrix(const std::string & header,
                     const Matrix      & C);

  /// \brief Write matrix to file.
  void WriteMatrixToFile(const std::string & filename,
                         const Matrix      & matrix,
                         LinalgFileFormat    format = MatrixAscii);

  /// \brief Write symmetric matrix to file.
  void WriteMatrixToFile(const std::string     & filename,
                         const SymmetricMatrix & matrix,
                         LinalgFileFormat        format = MatrixAscii);


  void WriteVectorToFile(const std::string & filename,
                         const Vector      & vector,
                         LinalgFileFormat    format = MatrixAscii);

// ================== TEMPLATE FUNCTION IMPLEMENTATION ======================

template <typename FS>
void Invert(flens::GeMatrix<FS> &A)
{
  assert(A.numRows()==A.numCols());

  flens::DenseVector<flens::Array<int> > p(A.numRows());

  int info = flens::trf(A, p);

  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
        << " of xxtrf call.";
    }
    else {  // info > 0
      oss << "Error in triangle factorization: The matrix is singular with its " << info
        << " diagonal element exactly equal to zero.";
    }
    throw Exception(oss.str());
  }

  info = flens::tri(A, p);

  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
        << " of xxtri call.";
    }
    else {  // info > 0
      oss << "Error in triangle inversion: The matrix is singular with its " << info
        << " diagonal element exactly equal to zero.";
    }
    throw Exception(oss.str());
  }

}

template <typename FS>
void TriangleFactorize(flens::GeMatrix<FS> &A, flens::DenseVector<flens::Array<int> > &p)
{
  assert(A.numRows()==A.numCols());

  int info = flens::trf(A, p);

  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
        << " of xxtrf call.";
    }
    else {  // info > 0
      oss << "Error in triangle factorization: The matrix is singular with its " << info
        << " diagonal element exactly equal to zero.";
    }
    throw Exception(oss.str());
  }

}

template <typename FS>
void TriangleInvert(flens::GeMatrix<FS> &A, flens::DenseVector<flens::Array<int> > &p)
{
  assert(A.numRows()==A.numCols());

  int info = flens::tri(A, p);

  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
        << " of xxtri call.";
    }
    else {  // info > 0
      oss << "Error in triangle inversion: The matrix is singular with its " << info
        << " diagonal element exactly equal to zero.";
    }
    throw Exception(oss.str());
  }

}

//------------------------------------
template <typename T>
void SetMatrixFrom2DArray(Matrix &  A,
                          T      ** F)
//------------------------------------
{
  for (int i=0 ; i < A.numRows() ; i++) {
    for (int j=0 ; j < A.numCols() ; j++) {
      A(i,j) = static_cast<double>( F[i][j] );
    }
  }
}

//------------------------------------------
template <typename T>
void Set2DArrayFromMatrix(const Matrix &  A,
                          T            ** F)
//------------------------------------------
{
  for (int i=0 ; i < A.numRows() ; i++) {
    for (int j=0 ; j < A.numCols() ; j++) {
      F[i][j] = static_cast<T>( A(i,j) );
    }
  }
}


} // namespace NRLib

#endif // NRLIB_FLENS_HPP

