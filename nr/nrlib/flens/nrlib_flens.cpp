// $Id: nrlib_flens.cpp 1283 2014-09-02 12:33:49Z perroe $

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

#include "nrlib_flens.hpp"

#include <sstream>

#include "../exception/exception.hpp"
#include "../iotools/fileio.hpp"
#include "../iotools/logkit.hpp"

using namespace NRLib;

void NRLib::CholeskySolve(const SymmetricMatrix & A,
                          const Vector          & b,
                          Vector                & x)
{
  SymmetricMatrix temp = A;
  Matrix B(b.length(), 1);
  B(flens::_, 0) = b;                   // First column of B

  int info = flens::posv(temp, B);
  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
          << " of posv call.";
    }
    else {  // info > 0
      oss << "Error in Cholesky: The leading minor of order " << info
          << " is not positive definite.";
    }
    throw Exception(oss.str());
  }

  x = B(flens::_, 0);
}

//--------------------------------------------------
void NRLib::CholeskySolve(const SymmetricMatrix & A,
                          Matrix                & B) // This is also where the solution is stored
//--------------------------------------------------
{
  SymmetricMatrix temp = A;

  int info = flens::posv(temp, B);
  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
          << " of posv call.";
    }
    else {  // info > 0
      oss << "Error in Cholesky: The leading minor of order " << info
          << " is not positive definite.";
    }
    throw Exception(oss.str());
  }
}

//---------------------------------------------------------
void NRLib::CholeskySolveComplex(ComplexMatrix & A,
                                 ComplexMatrix & B) // This is also where the solution is stored
//---------------------------------------------------------
{
  NRLib::ComplexMatrix::TriangularView::HermitianView H = A.lower().hermitian();

  int info = flens::posv(H, B);

  if (info != 0) {
    std::ostringstream oss;
    if (info < 0) {
      oss << "Internal FLENS/Lapack error: Error in argument " << -info
          << " of posv call.";
    }
    else {  // info > 0
      oss << "Error in Cholesky: The leading minor of order " << info
          << " is not positive definite.";
    }
    throw Exception(oss.str());
  }
}

//--------------------------------------------
void NRLib::CholeskyInvert(SymmetricMatrix& A)
//--------------------------------------------
{
  //NBNB legge til feilmld her, se metode ovenfor
  int infof = flens::potrf(A.upLo(), A.dim(), A.data(), A.leadingDimension());
  int infoi = flens::potri(A.upLo(), A.dim(), A.data(), A.leadingDimension());

  if(infof > 0 || infoi > 0)
    throw NRLib::Exception("Error in Cholesky inversion");
}

void NRLib::ComputeEigenVectors(Matrix &A,
                                Vector       &eigen_values,
                                Matrix       &eigen_vectors)
{
  Matrix dummy_mat(A.numRows(), A.numCols());
  Vector dummy_vec(A.numRows());
  flens::ev(false, true, A, eigen_values, dummy_vec, dummy_mat, eigen_vectors);
}

void NRLib::ComputeEigenVectorsComplex(ComplexMatrix        & A,
                                       ComplexVector        & eigen_values,
                                       ComplexMatrix        & eigen_vectors)
{
  ComplexMatrix dummy_mat(A.numRows(), A.numCols());
  flens::ev(false, true, A, eigen_values, dummy_mat, eigen_vectors);
}

void NRLib::CholeskyFactorize(SymmetricMatrix & A)
{
  int info = flens::potrf(A);

  if (info != 0) {
    std::ostringstream oss;
     oss << "Error in Cholesky: The leading minor of order " << info
          << " is not positive definite.";
     throw Exception(oss.str());
  }


}



void NRLib::ComputeEigenVectorsSymmetric(const SymmetricMatrix & A,
                                         Vector                & eigen_values,
                                         Matrix                & eigen_vectors)
{
  eigen_vectors = A;
  int info = flens::ev(true, eigen_vectors, eigen_values);
  assert (info >= 0);
  if (info > 0)
    throw NRLib::Exception("Error in eigenvaluecalculation: Algorithm failed to converge.");
}


//--------------------------------------------------------
double NRLib::FindLargestElement(const NRLib::Vector & v)
//-------------------------------------------------------
{
  double max_value = -std::numeric_limits<double>::infinity();
  for (int i=0 ; i < v.length() ; i++) {
    if (v(i) > max_value) {
      max_value = v(i);
    }
  }
  return max_value;
}

//--------------------------------------------------
void NRLib::InitializeMatrix(NRLib::Matrix & A,
                             double          value)
//--------------------------------------------------
{
  for (int i=0 ; i < A.numRows() ; i++) {
    for (int j=0 ; j < A.numCols() ; j++) {
      A(i,j) = value;
    }
  }
}

//--------------------------------------------------
void NRLib::InitializeComplexMatrix(NRLib::ComplexMatrix  & A,
                                    std::complex<double>    value)
//--------------------------------------------------
{
  for (int i=0 ; i < A.numRows() ; i++) {
    for (int j=0 ; j < A.numCols() ; j++) {
      A(i, j) = value;
    }
  }
}

//-------------------------------------------------------------------
void NRLib::InitializeSymmetricMatrix(NRLib::SymmetricMatrix & A,
                                      double                   value)
//-------------------------------------------------------------------
{
  for (int i=0 ; i < A.dim() ; i++) {
    for (int j=0 ; j <= i ; j++) {
      A(j,i) = value;
    }
  }
}

//------------------------------------
NRLib::Vector  NRLib::ZeroVector(int n)
//------------------------------------
{
  NRLib::Vector zero(n);
  zero = 0;
  return zero;
}

//------------------------------------
NRLib::Matrix NRLib::ZeroMatrix(int n)
//------------------------------------
{
  NRLib::Matrix Zero(n,n);
  NRLib::InitializeMatrix(Zero, 0.0);
  return Zero;
}

//------------------------------------------------------
NRLib::SymmetricMatrix NRLib::SymmetricZeroMatrix(int n)
//------------------------------------------------------
{
  NRLib::SymmetricMatrix Zero(n);
  NRLib::InitializeSymmetricMatrix(Zero, 0.0);
  return Zero;
}

//----------------------------------------
NRLib::Matrix NRLib::IdentityMatrix(int n)
//----------------------------------------
{
  NRLib::Matrix I = ZeroMatrix(n);
  for(int i=0 ; i<n ; i++) {
    I(i,i) = 1.0;
  }
  return I;
}

//---------------------------------------
void NRLib::Sort3x3(NRLib::Vector & Eval,
                    NRLib::Matrix & Evec)
//---------------------------------------
{
  NRLib::Vector h1(3);
  NRLib::Matrix H2(3,3);

  NRLib::IntegerVector index(3);
  for (int i=0 ; i < 3 ; i++) {
    index(i) = -1;
  }

  for (int j=0 ; j < 3 ; j++) {
    double max = -100000.0;
    for (int i=0 ; i < 3 ; i++) {
      if (j==0 || (j==1 && i!=index(0)) ||(j==2 && i!=index(0) && i!=index(1))) {
        if (Eval(i) > max) {
          max = Eval(i);
          index(j) = i;
        }
      }
    }
    h1(j) = max;

    for (int k=0 ; k < 3 ; k++) {
      H2(k,j) = Evec(k, index(j));
    }
  }

  for (int i=0 ; i < 3 ; i++) {
    Eval(i) = h1(i);
    for (int j=0 ; j < 3 ; j++)
      Evec(i,j) = H2(i,j);
  }
}


NRLib::Matrix NRLib::ReadMatrixFromFile(const std::string & filename,
                                        int                 nRows,
                                        int                 nCols,
                                        LinalgFileFormat    format)
{
  NRLib::Matrix m(nRows, nCols);
  std::ifstream file;

  if (format == MatrixAscii) {
    NRLib::OpenRead(file, filename);
    int line;
    for (int i = 0; i < nRows; ++i) {
      for (int j = 0; j < nCols; ++j) {
        m(i, j) = NRLib::ReadNext<double>(file, line);
      }
    }
  }
  else {
    NRLib::OpenRead(file, filename, std::ios::in | std::ios::binary);
    for (int i = 0; i < nRows; ++i) {
      for (int j = 0; j < nCols; ++j) {
        m(i, j) = NRLib::ReadBinaryDouble(file);
      }
    }
  }

  if (!CheckEndOfFile(file))
    throw NRLib::Exception("Too much data in file.");

  return m;
}


NRLib::Vector NRLib::ReadVectorFromFile(const std::string & filename,
                                        int                 n,
                                        LinalgFileFormat    format)
{
  NRLib::Vector v(n);
  std::ifstream file;

  if (format == MatrixAscii) {
    NRLib::OpenRead(file, filename);
    int line;
    for (int i = 0; i < n; ++i)
      v(i) = ReadNext<double>(file, line);
  }
  else if (format == MatrixBinary) {
    NRLib::OpenRead(file, filename, std::ios::in | std::ios::binary);
    for (int i = 0; i < n; ++i)
      v(i) = ReadBinaryDouble(file);
  }

  if (!CheckEndOfFile(file))
    throw NRLib::Exception("Too much data in file.");

  return v;
}

void NRLib::WriteComplexVector(const std::string   & header,
                               const ComplexVector & c)
{
  int n = c.length();
  LogKit::LogFormatted(LogKit::Error,"\n"+header+"\n");
  for (int i=0; i < n ; i++) {
    LogKit::LogFormatted(LogKit::Error,"(%12.8f, %12.8f)\n",c(i).real(),c(i).imag());
  }
  LogKit::LogFormatted(LogKit::Error,"\n");
}


void NRLib::WriteComplexMatrix(const std::string   & header,
                               const ComplexMatrix & C)
{
  int m = C.numRows();
  int n = C.numCols();
  LogKit::LogFormatted(LogKit::Error,"\n"+header+"\n");
  for (int i=0; i < m ; i++) {
    for (int j=0; j < n ; j++) {
      LogKit::LogFormatted(LogKit::Error,"(%12.8f, %12.8f) ",C(i,j).real(),C(i,j).imag());
    }
    LogKit::LogFormatted(LogKit::Error,"\n");
  }
  LogKit::LogFormatted(LogKit::Error,"\n");
}


void NRLib::Adjoint(const ComplexMatrix & in,
                    ComplexMatrix       & out)
{
  int m = out.numRows();
  int n = out.numCols();

  for (int i=0 ; i < m ; i++) {
    for (int j=0 ; j < n ; j++) {
      out(i,j) = std::conj(in(j,i));
    }
  }
}


void NRLib::WriteMatrix(const std::string & header,
                        const Matrix      & C)
{
  int m = C.numRows();
  int n = C.numCols();
  LogKit::LogFormatted(LogKit::Error,"\n"+header+"\n");
  for (int i=0; i < m ; i++) {
    for (int j=0; j < n ; j++) {
      LogKit::LogFormatted(LogKit::Error,"%20.8f ",C(i,j));
    }
    LogKit::LogFormatted(LogKit::Error,"\n");
  }
  LogKit::LogFormatted(LogKit::Error,"\n");
}


void NRLib::WriteMatrixToFile(const std::string   & filename,
                              const NRLib::Matrix & matrix,
                              LinalgFileFormat      format)
{
  std::ofstream file;

  if (format == MatrixAscii) {
    NRLib::OpenWrite(file, filename);
    file.precision(8);
    file.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < matrix.numRows(); ++i) {
      for (int j = 0; j < matrix.numCols(); ++j) {
        file.width(18);
        file << matrix(i, j);
      }
      file << "\n";
    }
  }
  else if (format == MatrixBinary) {
    NRLib::OpenWrite(file, filename, std::ios::out | std::ios::binary);
    for (int i = 0; i < matrix.numRows(); ++i) {
      for (int j = 0; j < matrix.numCols(); ++j) {
        NRLib::WriteBinaryDouble(file, matrix(i, j));
      }
    }
  }
}


void NRLib::WriteMatrixToFile(const std::string            & filename,
                              const NRLib::SymmetricMatrix & matrix,
                              LinalgFileFormat               format)
{
  std::ofstream file;

  if (format == MatrixAscii) {
    NRLib::OpenWrite(file, filename);
    file.precision(8);
    file.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < matrix.dim(); ++i) {
      for (int j = 0; j < matrix.dim(); ++j) {
        file.width(18);
        int I = (matrix.upLo()==flens::Upper) ? std::min(i,j) : std::max(i,j);
        int J = (matrix.upLo()==flens::Upper) ? std::max(i,j) : std::min(i,j);
        file << matrix(I, J);
      }
      file << "\n";
    }
  }
  else if (format == MatrixBinary) {
    NRLib::OpenWrite(file, filename, std::ios::out | std::ios::binary);
    for (int i = 0; i < matrix.dim(); ++i) {
      for (int j = 0; j < matrix.dim(); ++j) {
        int I = (matrix.upLo()==flens::Upper) ? std::min(i,j) : std::max(i,j);
        int J = (matrix.upLo()==flens::Upper) ? std::max(i,j) : std::min(i,j);
        NRLib::WriteBinaryDouble(file, matrix(I, J));
      }
    }
  }
}


void NRLib::WriteVectorToFile(const std::string   & filename,
                              const NRLib::Vector & vector,
                              LinalgFileFormat      format)
{
  std::ofstream file;

  if (format == MatrixAscii) {
    NRLib::OpenWrite(file, filename);
    file.precision(8);
    file.setf(std::ios_base::scientific, std::ios_base::floatfield);
    for (int i = 0; i < vector.length(); ++i) {
      file.width(18);
      file << vector(i);
    }
  }
  else if (format == MatrixBinary) {
    NRLib::OpenWrite(file, filename, std::ios::out | std::ios::binary);
    for (int i = 0; i < vector.length(); ++i) {
      NRLib::WriteBinaryDouble(file, vector(i));
    }
  }
}

