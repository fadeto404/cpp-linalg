//
// Created by rns on 6/13/21.
//

#ifndef CPP_LINALG_MATRIX_H
#define CPP_LINALG_MATRIX_H

#include "mem_blocks.h"

namespace cla
{
    template<int n, int m, class memBlockType = DenseMatrixBlock<n, m, double> >
    struct Matrix
    {
        const static int rows = n;
        const static int cols = m;
        memBlockType elements;

        // Constructors
        // TODO: Make constructor given same type
        // TODO: Constructor given memory block
        // TODO: Constructor given array formatted data

        // Assignment
        Matrix<n, m, memBlockType>& operator = (const Matrix<n, m, memBlockType>& mat);
        Matrix<n, m, memBlockType>& operator = (typename memBlockType::elementType data_array[n][m]);
        Matrix<n, m, memBlockType>& operator = (const typename memBlockType::elementType& value);

        // Dimensions
        static int getNumRows() const;
        static int getNumCols() const;
        static int getNumElements() const;

        // Element access
        memBlockType::elementType& operator () (int i, int j) const;

        // Matrix addition/subtraction
        Matrix<n, m, memBlockType>& operator + (const Matrix<n, m, memBlockType>& mat) const;
        Matrix<n, m, memBlockType>& operator - (const Matrix<n, m, memBlockType>& mat) const;

        Matrix<n, m, memBlockType>& operator += (const Matrix<n, m, memBlockType>& mat);
        Matrix<n, m, memBlockType>& operator -= (const Matrix<n, m, memBlockType>& mat);

        // Negation
        Matrix<n, m, memBlockType>& operator - () const;

        // Matrix multiplication
        template<int p, class otherMemBlockType>
        Matrix<n, m, memBlockType>& operator * (const Matrix<m, p, otherMemBlockType>& mat) const;

        // Only for square matrices (n == m) (n x m multiplied by m x p yields n x p)
        template<class otherMemBlockType>
        Matrix<n, n, memBlockType>& operator *= (const Matrix<n, n, otherMemBlockType>& mat);

        // Scalar operations
        Matrix<n, m, memBlockType>& operator + (const typename memBlockType::elementType& scalar) const;
        Matrix<n, m, memBlockType>& operator - (const typename memBlockType::elementType& scalar) const;
        Matrix<n, m, memBlockType>& operator * (const typename memBlockType::elementType& scalar) const;
        Matrix<n, m, memBlockType>& operator / (const typename memBlockType::elementType& scalar) const;
        Matrix<n, m, memBlockType>& operator += (const typename memBlockType::elementType& scalar);
        Matrix<n, m, memBlockType>& operator -= (const typename memBlockType::elementType& scalar);
        Matrix<n, m, memBlockType>& operator *= (const typename memBlockType::elementType& scalar);
        Matrix<n, m, memBlockType>& operator /= (const typename memBlockType::elementType& scalar);

    };

    // DIMENSIONS
    template<int n, int m, class memBlockType>
    int Matrix<n, m, memBlockType>::getNumRows() const
    {
        return rows;
    }

    template<int n, int m, class memBlockType>
    int Matrix<n, m, memBlockType>::getNumCols() const
    {
        return cols;
    }

    template<int n, int m, class memBlockType>
    int Matrix<n, m, memBlockType>::getNumElements() const
    {
        return rows*cols;
    }

    // ASSIGNMENT
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator = (typename memBlockType::elementType data_array[n][m])
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) = data_array[i][j];
            }
        }
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator = (const typename memBlockType::elementType &value)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) = value;
            }
        }
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator = (const Matrix<n, m, memBlockType> &mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) = mat(i, j);
            }
        }
    }

    // ELEMENT ACCESS
    template<int n, int m, class memBlockType>
    inline memBlockType::elementType& Matrix<n, m, memBlockType>::operator()(int i, int j) const
    {
        return elements(i, j);
    }

    // ARITHMETIC FUNCTIONS

    // ADDITION
    // TODO: Replace with general add function
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator + (const Matrix<n, m, memBlockType> &mat) const
    {
        Matrix<n, m, memBlockType> res;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                res(i, j) = elements(i, j) + mat(i, j);
            }
        }
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator += (const Matrix<n, m, memBlockType> &mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) += mat(i, j);
            }
        }
    }
}

#endif //CPP_LINALG_MATRIX_H
