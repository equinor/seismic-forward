/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLENS_FULLSTORAGE_H
#define FLENS_FULLSTORAGE_H 1

#include <flens/array.h>
#include <flens/storage.h>

#ifndef FLENS_FIRST_INDEX
#    define FLENS_FIRST_INDEX 1
#endif

namespace flens {

// forward declarations
template <typename T, StorageOrder Order>
class FullStorageView;

template <typename T, StorageOrder Order>
class FullStorage;


template <typename T, StorageOrder Order>
class ConstFullStorageView
{
    public:
        typedef const T                        ElementType;
        typedef ConstFullStorageView           ConstView;
        typedef FullStorageView<T, Order>      View;
        typedef FullStorage<T, Order>          NoView;

        typedef ConstArrayView<T>              ConstVectorView;
        typedef ArrayView<T>                   VectorView;
        typedef Array<T>                       VectorNoView;

        ConstFullStorageView(const T *data,
                             int numRows, int numCols,
                             int leadingDimension,
                             int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

        ConstFullStorageView(const ConstFullStorageView<T, Order> &rhs);

        ~ConstFullStorageView();

        //-- operators ---------------------------------------------------------

        const T &
        operator()(int row, int col) const;

        //-- methods -----------------------------------------------------------

        int
        firstRow() const;

        int
        firstCol() const;

        int
        lastRow() const;

        int
        lastCol() const;

        int
        numRows() const;

        int
        numCols() const;

        int
        leadingDimension() const;

        int
        strideRow() const;

        int
        strideCol() const;

        const T *
        data() const;

        T *
        data();

        ConstFullStorageView<T, Order>
        view(int fromRow, int fromCol, int toRow, int toCol,
             int firstViewRow=FLENS_FIRST_INDEX, int firstViewCol=FLENS_FIRST_INDEX) const;

        ConstArrayView<T>
        viewRow(int row, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ConstArrayView<T>
        viewCol(int col, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ConstArrayView<T>
        viewRow(int row, int fromCol, int toCol, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ConstArrayView<T>
        viewCol(int col, int fromRow, int toRow, int firstViewIndex=FLENS_FIRST_INDEX) const;

    private:
        //const T * const *
        //leadingDimensionStorage() const;

        //const T **
        //leadingDimensionStorage();

        //void
        //_allocate(const T *data);

        //const T **_data;
        const T *_data;
        int _numRows, _numCols;
        int _leadingDimension;
        int _firstRow, _firstCol;
};

template <typename T, StorageOrder Order>
class FullStorageView
{
    public:
        typedef T                              ElementType;
        typedef ConstFullStorageView<T, Order> ConstView;
        typedef FullStorageView                View;
        typedef FullStorage<T, Order>          NoView;

        typedef ConstArrayView<T>              ConstVectorView;
        typedef ArrayView<T>                   VectorView;
        typedef Array<T>                       VectorNoView;


        FullStorageView(T *data,
                        int numRows, int numCols,
                        int leadingDimension,
                        int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

        FullStorageView(const FullStorageView<T, Order> &rhs);

        ~FullStorageView();

        //-- operators ---------------------------------------------------------

        FullStorageView<T, Order> &
        operator=(const FullStorage<T, Order> &rhs);

        FullStorageView<T, Order> &
        operator=(const FullStorageView<T, Order> &rhs);

        FullStorageView<T, Order> &
        operator=(const ConstFullStorageView<T, Order> &rhs);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        operator ConstView() const;

        //-- methods -----------------------------------------------------------

        int
        firstRow() const;

        int
        firstCol() const;

        int
        lastRow() const;

        int
        lastCol() const;

        int
        numRows() const;

        int
        numCols() const;

        int
        leadingDimension() const;

        int
        strideRow() const;

        int
        strideCol() const;

        const T *
        data() const;

        T *
        data();

        ConstFullStorageView<T, Order>
        view(int fromRow, int fromCol, int toRow, int toCol,
             int firstViewRow=FLENS_FIRST_INDEX, int firstViewCol=FLENS_FIRST_INDEX) const;

        FullStorageView<T, Order>
        view(int fromRow, int fromCol, int toRow, int toCol,
             int firstViewRow=FLENS_FIRST_INDEX, int firstViewCol=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewRow(int row, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewRow(int row, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewCol(int col, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewCol(int col, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewRow(int row, int fromCol, int toCol, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewRow(int row, int fromCol, int toCol, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewCol(int col, int fromRow, int toRow, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewCol(int col, int fromRow, int toRow, int firstViewIndex=FLENS_FIRST_INDEX);

                ArrayView<T>
                viewDiag(int d);

                ConstArrayView<T>
                viewDiag(int d) const;

        void
        resize(int numRows, int numCols, int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

    private:
        //const T *const *
        //leadingDimensionStorage() const;

        //T**
        //leadingDimensionStorage();

        //void
        //_allocate(T *data);

        //T **_data;
        T* _data;
        int _numRows, _numCols;
        int _leadingDimension;
        int _firstRow, _firstCol;
};

template <typename T, StorageOrder Order>
class FullStorage
{
    public:
        typedef T                              ElementType;
        typedef ConstFullStorageView<T, Order> ConstView;
        typedef FullStorageView<T, Order>      View;
        typedef FullStorage                    NoView;

        typedef ConstArrayView<T>              ConstVectorView;
        typedef ArrayView<T>                   VectorView;
        typedef Array<T>                       VectorNoView;

        FullStorage();

        FullStorage(int numRows, int numCols, int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

        FullStorage(const FullStorage<T, Order> &rhs);

        FullStorage(const FullStorageView<T, Order> &rhs);

        FullStorage(const ConstFullStorageView<T, Order> &rhs);

        ~FullStorage();

        //-- operators ---------------------------------------------------------

        FullStorage<T, Order> &
        operator=(const FullStorage<T, Order> &rhs);

        FullStorage<T, Order> &
        operator=(const FullStorageView<T, Order> &rhs);

        FullStorage<T, Order> &
        operator=(const ConstFullStorageView<T, Order> &rhs);

        const T &
        operator()(int row, int col) const;

        T &
        operator()(int row, int col);

        //-- methods -----------------------------------------------------------

        int
        firstRow() const;

        int
        firstCol() const;

        int
        lastRow() const;

        int
        lastCol() const;

        int
        numRows() const;

        int
        numCols() const;

        int
        leadingDimension() const;

        int
        strideRow() const;

        int
        strideCol() const;

        const T *
        data() const;

        T *
        data();

        void
        resize(int numRows, int numCols, int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

        void
        resizeOrClear(int numRows, int numCols, int firstRow=FLENS_FIRST_INDEX, int firstCol=FLENS_FIRST_INDEX);

        void
        initialize(T val = T());

        void
        shiftIndexTo(int firstRow, int firstCol);

        ConstFullStorageView<T, Order>
        view(int fromRow, int fromCol, int toRow, int toCol,
             int firstViewRow=FLENS_FIRST_INDEX, int firstViewCol=FLENS_FIRST_INDEX) const;

        FullStorageView<T, Order>
        view(int fromRow, int fromCol, int toRow, int toCol,
             int firstViewRow=FLENS_FIRST_INDEX, int firstViewCol=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewRow(int row, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewRow(int row, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewCol(int col, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewCol(int col, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewRow(int row, int fromCol, int toCol, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewRow(int row, int fromCol, int toCol, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewCol(int col, int fromRow, int toRow, int firstViewIndex=FLENS_FIRST_INDEX) const;

        ArrayView<T>
        viewCol(int col, int fromRow, int toRow, int firstViewIndex=FLENS_FIRST_INDEX);

        ArrayView<T>
        viewDiag(int d, int firstViewIndex=FLENS_FIRST_INDEX);

        ConstArrayView<T>
        viewDiag(int d, int firstViewIndex=FLENS_FIRST_INDEX) const;

    private:
        const T * const *
        leadingDimensionStorage() const;

        T**
        leadingDimensionStorage();

        void
        _allocate();

        void
        _release();

        int _numRows, _numCols, _firstRow, _firstCol;
        T **_data;
};


//------------------------------------------------------------------------------

template <typename T, StorageOrder Order>
struct StorageInfo<FullStorage<T, Order> >
{
    static const StorageOrder order = Order;
};

template <typename T, StorageOrder Order>
struct StorageInfo<FullStorageView<T, Order> >
{
    static const StorageOrder order = Order;
};

template <typename T, StorageOrder Order>
struct StorageInfo<ConstFullStorageView<T, Order> >
{
    static const StorageOrder order = Order;
};


} // namespace flens

#include <flens/fullstorage.tcc>

#endif // FLENS_FULLSTORAGE_H
