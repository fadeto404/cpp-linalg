#ifndef CPP_LINALG_CLA_MATRIX_CORE_H
#define CPP_LINALG_CLA_MATRIX_CORE_H

#include "mem_blocks.h"
#include <cmath>

namespace cla
{
    template<int n, int m, class memBlockType=DenseMatrixBlock<n, m, float> >
    struct Matrix
    {
        const static int rows = n;
        const static int cols = m;
        memBlockType elements;

        // Constructors
        Matrix<n, m, memBlockType>(){};
        Matrix<n, m, memBlockType>(memBlockType &mem_block) : elements(mem_block){};
        Matrix<n, m, memBlockType>(typename memBlockType::element_type data[n][m]){*this = data;};
        Matrix<n, m, memBlockType>(typename memBlockType::element_type data[n*m]){*this = data;};

        // Code written by tomstewart89
        template<typename ...ARGS> Matrix(ARGS... args) { FillRowMajor(args...); }

        // Assignment
        Matrix<n, m, memBlockType>& operator =(const Matrix<n, m, memBlockType>& mat);
        Matrix<n, m, memBlockType>& operator =(typename memBlockType::element_type data_array[n][m]);
        Matrix<n, m, memBlockType>& operator =(const typename memBlockType::element_type& value);

        // Code written by tomstewart89
        template<typename ...TAIL> void FillRowMajor(typename memBlockType::element_type head, TAIL... tail);
        void FillRowMajor() { }

        // Dimensions
        static const int getNumRows();
        static const int getNumCols();
        static const int getNumElements();

        // Element access
        typename memBlockType::element_type & operator ()(int i, int j) const;

        // Matrix addition
        template<class otherMemBlockType>
        Matrix<n, m, memBlockType> operator +(const Matrix<n, m, otherMemBlockType>& mat) const;

        template<class otherMemBlockType>
        Matrix<n, m, memBlockType>& operator +=(const Matrix<n, m, otherMemBlockType>& mat);

        // Matrix subtraction
        template<class otherMemBlockType>
        Matrix<n, m, memBlockType> operator -(const Matrix<n, m, otherMemBlockType>& mat) const;

        template<class otherMemBlockType>
        Matrix<n, m, memBlockType>& operator -=(const Matrix<n, m, otherMemBlockType>& mat);

        // Negation
        Matrix<n, m, memBlockType>& operator -() const;

        // Matrix multiplication
        template<int p, class otherMemBlockType>
        Matrix<n, p, DenseMatrixBlock<n, p, typename memBlockType::element_type> > operator *(const Matrix<m, p, otherMemBlockType>& mat) const;

        // Only for square matrices (n == m) (n x m multiplied by m x p yields n x p)
        template<class otherMemBlockType>
        Matrix<n, n, memBlockType>& operator *=(const Matrix<n, n, otherMemBlockType>& mat);

        // Scalar operations
        Matrix<n, m, memBlockType>& operator +(const typename memBlockType::element_type & scalar) const;
        Matrix<n, m, memBlockType>& operator -(const typename memBlockType::element_type & scalar) const;
        Matrix<n, m, memBlockType>& operator *(const typename memBlockType::element_type & scalar) const;
        Matrix<n, m, memBlockType>& operator /(const typename memBlockType::element_type & scalar) const;
        Matrix<n, m, memBlockType>& operator +=(const typename memBlockType::element_type & scalar);
        Matrix<n, m, memBlockType>& operator -=(const typename memBlockType::element_type & scalar);
        Matrix<n, m, memBlockType>& operator *=(const typename memBlockType::element_type & scalar);
        Matrix<n, m, memBlockType>& operator /=(const typename memBlockType::element_type & scalar);

        // Transpose
        Matrix<m, n, TransposeBlock<memBlockType> > T();
    };

    // Matrix arithmetic functions
    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& add(const Matrix<n, m, memBlockType> &A,
                                          const Matrix<n, m, otherMemBlockType> &B,
                                          Matrix<n, m, resultMemBlockType> &C);

    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& subtract(const Matrix<n, m, memBlockType> &A,
                                               const Matrix<n, m, otherMemBlockType> &B,
                                               Matrix<n, m, resultMemBlockType> &C);

    template<int n, int m, int p, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, p, resultMemBlockType>& multiply(const Matrix<n, m, memBlockType> &A,
                                               const Matrix<m, p, otherMemBlockType> &B,
                                               Matrix<n, p, resultMemBlockType> &C);

    template<int n, int m, class memBlockType>
    Matrix<m, n, TransposeBlock<memBlockType> > transpose(Matrix<n, m, memBlockType> &mat);

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> transpose(Matrix<m, n, TransposeBlock<memBlockType> > &transposed_mat);

    /////// IMPLEMENTATIONS

    // DIMENSIONS
    template<int n, int m, class memBlockType>
    const int Matrix<n, m, memBlockType>::getNumRows()
    {
        return rows;
    }

    template<int n, int m, class memBlockType>
    const int Matrix<n, m, memBlockType>::getNumCols()
    {
        return cols;
    }

    template<int n, int m, class memBlockType>
    const int Matrix<n, m, memBlockType>::getNumElements()
    {
        return rows*cols;
    }

    // ASSIGNMENT
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator =(typename memBlockType::element_type data_array[n][m])
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) = data_array[i][j];
            }
        }
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator =(const typename memBlockType::element_type &value)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                elements(i, j) = value;
            }
        }
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator =(const Matrix<n, m, memBlockType> &mat)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                (*this)(i, j) = mat(i, j);
            }
        }
        return (*this);
    }

    // Code written by tomstewart89
    template<int n, int m, class memBlockType>
    template<typename... TAIL>
    void Matrix<n, m, memBlockType>::FillRowMajor(typename memBlockType::element_type head, TAIL... tail) {
        static_assert(rows*cols > sizeof...(TAIL), "Too many arguments passed to FillRowMajor");
        (*this)((rows*cols - sizeof...(TAIL) - 1) / cols,(rows*cols - sizeof...(TAIL) - 1) % cols) = head;
        FillRowMajor(tail...);
    }

    // ELEMENT ACCESS
    template<int n, int m, class memBlockType>
    inline typename memBlockType::element_type& Matrix<n, m, memBlockType>::operator()(int i, int j) const
    {
        return elements(i, j);
    }

    /////////////////////////////// MATRIX ARITHMETIC OPERATIONS ///////////////////////////////
    // MATRIX ADDITION
    template<int n, int m, class memBlockType>
    template<class otherMemBlockType>
    Matrix<n, m, memBlockType> Matrix<n, m, memBlockType>::operator +(const Matrix<n, m, otherMemBlockType> &mat) const
    {
        Matrix<n, m, memBlockType> res;
        add((*this), mat, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    template<class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +=(const Matrix<n, m, otherMemBlockType> &mat)
    {
        add((*this), mat, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& add(const Matrix<n, m, memBlockType> &A,
                                          const Matrix<n, m, otherMemBlockType> &B,
                                          Matrix<n, m, resultMemBlockType> &C)
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
    template<int n, int m, class memBlockType>
    template<class otherMemBlockType>
    Matrix<n, m, memBlockType> Matrix<n, m, memBlockType>::operator -(const Matrix<n, m, otherMemBlockType> &mat) const
    {
        Matrix<n, m, memBlockType> res;
        subtract((*this), mat, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    template<class otherMemBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -=(const Matrix<n, m, otherMemBlockType> &mat)
    {
        subtract((*this), mat, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType, class otherMemBlockType, class resultMemBlockType>
    Matrix<n, m, resultMemBlockType>& subtract(const Matrix<n, m, memBlockType> &A,
                                               const Matrix<n, m, otherMemBlockType> &B,
                                               Matrix<n, m, resultMemBlockType> &C)
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
    Matrix<n, p, DenseMatrixBlock<n, p, typename memBlockType::element_type> > Matrix<n, m, memBlockType>::operator *(const Matrix<m, p, otherMemBlockType> &mat) const
    {
        Matrix<n, p, DenseMatrixBlock<n, p, typename memBlockType::element_type> > res;
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
    Matrix<n, p, resultMemBlockType>& multiply(const Matrix<n, m, memBlockType> &A,
                                               const Matrix<m, p, otherMemBlockType> &B,
                                               Matrix<n, p, resultMemBlockType> &C)
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
    // TODO: Elementwise matrix arithmetic ( C(i,j) = A(i,j) */+- B(i,j) )

    /////////////////////////////// MATRIX-SCALAR OPERATIONS ///////////////////////////////
    // MATRIX-SCALAR ADDITION
    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +(const typename memBlockType::element_type &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        add((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator +=(const typename memBlockType::element_type &scalar)
    {
        add((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m , memBlockType>& add(Matrix<n, m, memBlockType> &mat,
                                     const typename memBlockType::element_type &scalar,
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
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -(const typename memBlockType::element_type &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        subtract((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& Matrix<n, m, memBlockType>::operator -=(const typename memBlockType::element_type &scalar)
    {
        subtract((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m , memBlockType>& subtract(Matrix<n, m, memBlockType> &mat,
                                          const typename memBlockType::element_type &scalar,
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
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator *(const typename memBlockType::element_type &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        multiply((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator *=(const typename memBlockType::element_type &scalar)
    {
        multiply((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& multiply(Matrix<n, m, memBlockType> &mat,
                                         const typename memBlockType::element_type &scalar,
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
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator /(const typename memBlockType::element_type &scalar) const
    {
        Matrix<n, m, memBlockType> res;
        divide((*this), scalar, res);
        return res;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> &Matrix<n, m, memBlockType>::operator /=(const typename memBlockType::element_type &scalar)
    {
        divide((*this), scalar, (*this));
        return (*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType>& divide(Matrix<n, m, memBlockType> &mat,
                                       const typename memBlockType::element_type &scalar,
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

    ///////////////////////////////  TRANSPOSITION  ///////////////////////////////

    // Member function transpose
    template<int n, int m, class memBlockType>
    Matrix<m, n, TransposeBlock<memBlockType> > Matrix<n, m, memBlockType>::T()
    {
        return transpose(*this);
    }

    template<int n, int m, class memBlockType>
    Matrix<m, n, TransposeBlock<memBlockType> > transpose(Matrix<n, m, memBlockType> &mat)
    {
        TransposeBlock<memBlockType> t_blk(mat.elements);
        Matrix<m, n, TransposeBlock<memBlockType> > t(t_blk);
        return t;
    }

    template<int n, int m, class memBlockType>
    Matrix<n, m, memBlockType> transpose(Matrix<m, n, TransposeBlock<memBlockType> > &transposed_mat)
    {
        memBlockType orig_blk(transposed_mat.elements);
        Matrix<n, m, memBlockType> orig(orig_blk);
        return orig;
    }

    /////////////////////////////// EXTRA FUNCTIONS ///////////////////////////////
    // INNER PRODUCTS
    template<int n, class memBlockType, class otherMemBlockType>
    typename memBlockType::element_type dot(Matrix<n, 1, memBlockType> &a, Matrix<n, 1, otherMemBlockType> &b)
    {
        typename memBlockType::element_type result = 0;
        for (int i = 0; i < n; ++i)
        {
            result += a(i, 0)*b(i, 0);
        }
        return result;
    }

    // SAXPY, GAXPY, ETC.
    template<int n, class memBlockType, class otherMemBlockType = memBlockType>
    Matrix<n, 1, memBlockType>& saxpy(Matrix<n, 1, memBlockType> &y,
                                     typename memBlockType::element_type a,
                                     Matrix<n, 1, otherMemBlockType> &x)
    {
        for (int i = 0; i < n; ++i)
        {
            y(i,0) += a*x(i,0);
        }
        return y;
    }

    template<int n, int m, class memBlockType, class scalingMemBlockType, class otherMemBlockType>
    Matrix<n, 1, memBlockType>& gaxpy(Matrix<n, 1, memBlockType> &y,
                                      Matrix<m, n, scalingMemBlockType> &A,
                                      Matrix<m, 1, otherMemBlockType> &x)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                y(i,0) += A(i, j) * x(j,0);
            }
        }
        return y;
    }
}

#endif //CPP_LINALG_CLA_MATRIX_CORE_H
