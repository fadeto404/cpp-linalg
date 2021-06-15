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

        // Matrix addition
        Matrix<n, m, memBlockType>& operator + (const Matrix<n, m, memBlockType>& mat) const;
        Matrix<n, m, memBlockType>& operator += (const Matrix<n, m, memBlockType>& mat);

        // Matrix subtraction
        Matrix<n, m, memBlockType>& operator - (const Matrix<n, m, memBlockType>& mat) const;
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
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator =(typename memBlockType::elementType data_array[n][m])
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
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator =(const typename memBlockType::elementType &value)
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
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator =(const Matrix<n, m, memBlockType> &mat)
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

    /////////////////////////////// MATRIX ARITHMETIC OPERATIONS ///////////////////////////////
    // MATRIX ADDITION
    template<int n, int m, class memBlockType, class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +(const Matrix<n, m, otherMemBlockType> &mat) const
    {
        Matrix<n, m, memBlockType> res;
        add((*this), mat, res);
        return res;
    }

    template<int n, int m, class memBlockType, class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +=(const Matrix<n, m, otherMemBlockType> &mat)
    {
        add((*this), mat, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& add(const Matrix<n, m, memBlockType> &A,
                                    const Matrix<n, m, otherMemBlockType> &B,
                                    const Matrix<n, m, resultMemBlockType> &C)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }


    // MATRIX SUBTRACTION
    template<int n, int m, class memBlockType, class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -(const Matrix<n, m, otherMemBlockType> &mat) const
    {
        Matrix<n, m, memBlockType> res;
        subtract((*this), mat, res);
        return res;
    }

    template<int n, int m, class memBlockType, class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -=(const Matrix<n, m, otherMemBlockType> &mat)
    {
        subtract((*this), mat, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& subtract(const Matrix<n, m, memBlockType> &A,
                                    const Matrix<n, m, otherMemBlockType> &B,
                                    const Matrix<n, m, resultMemBlockType> &C)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                C(i, j) = A(i, j) - B(i, j);
            }
        }
        return C;
    }

    // NEGATION
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator -() const
    {
        Matrix<n, m, memBlockType> res;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                res(i, j) = -elements(i, j);
            }
        }
    }

    // MATRIX MULTIPLICATION
    template<int n, int m, class memBlockType>  // Template for the left side matrix
    template<int p, class otherMemBlockType>    // Template for the right side matrix
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator *(const Matrix<m, p, otherMemBlockType> &mat) const
    {
        Matrix<n, p, memBlockType> res;
        multiply((*this), mat, res);
        return res;
    }

    template<int n, int m, class memBlockType>  // Template for the left operand
    template<class otherMemBlockType>           // Template for the right operand
    Matrix<n, n, memBlockType>& Matrix<n, m, memBlockType>::operator *=(const Matrix<n, n, otherMemBlockType> &mat)
    {
        Matrix<n, n, memBlockType> res;
        multiply((*this), mat, res);
        (*this) = res;
        return (*this);
    }

    template<int n, int m, int p, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& multiply(const Matrix<n, m, memBlockType> &A,
                                         const Matrix<n, m, otherMemBlockType> &B,
                                         const Matrix<n, m, resultMemBlockType> &C)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < p; ++j)
            {
                C(i, j) = A(i, 0) * B(0, j);
                for (int k = 1; k < m; ++k)
                {
                    C(i, j) += A(i, k) * B(k, j);
                }
            }
        }
        return C;
    }
    // TODO: Elementwise matrix multiplication ( C(i,j) = A(i,j) * B(i,j) )

    /////////////////////////////// MATRIX-SCALAR OPERATIONS ///////////////////////////////
    // MATRIX-SCALAR ADDITION
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +(const elementType &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        add((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +=(const elementType &scalar)
    {
        add((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m , memBlockType>& add(Matrix<n, m, memBlockType> &mat,
                                     const elementType &scalar,
                                     Matrix<n, m, memBlockType> &ret_mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                ret_mat(i, j) = mat(i, j) + scalar;
            }
        }
        return mat;
    }

    // MATRIX-SCALAR SUBTRACTION
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -(const elementType &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        subtract((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -=(const elementType &scalar)
    {
        subtract((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m , memBlockType>& subtract(Matrix<n, m, memBlockType> &mat,
                                          const elementType &scalar,
                                          Matrix<n, m, memBlockType> &ret_mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                ret_mat(i, j) = mat(i, j) - scalar;
            }
        }
        return mat;
    }

    // MATRIX-SCALAR MULTIPLICATION
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator *(const elementType &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        multiply((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator *=(const elementType &scalar)
    {
        multiply((*this), scalar, (*this));
        return (*this);
    }

    Matrix<n, m, memBlockType>& multiply(Matrix<n, m, memBlockType> &mat,
                                         const elementType &scalar,
                                         Matrix<n, m, memBlockType> &ret_mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                ret_mat(i, j) = mat(i, j) * scalar;
            }
        }
        return mat;
    }

    // MATRIX-SCALAR DIVISION (ELEMENTWISE DIVISION)
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator /(const elementType &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        divide((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator /=(const elementType &scalar)
    {
        divide((*this), scalar, (*this));
        return (*this);
    }

    Matrix<n, m, memBlockType>& divide(Matrix<n, m, memBlockType> &mat,
                                         const elementType &scalar,
                                         Matrix<n, m, memBlockType> &ret_mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                ret_mat(i, j) = mat(i, j) / scalar;
            }
        }
        return mat;
    }
}

#endif //CPP_LINALG_MATRIX_H
