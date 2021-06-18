//
// Created by rns on 6/16/21.
//

#ifndef CPP_LINALG_NORM_H
#define CPP_LINALG_NORM_H

#include "matrix_core.h"
#include "util.h"
#include <cmath>

namespace cla
{
    // TODO: Support other return types than double

    // Column vector p-norm
    template<int n, class memBlockType>
    double norm(Matrix<n, 1, memBlockType> &mat, int p = 2)
    {
        double l = 0;
        for (int i = 0; i < mat.rows; ++i)
        {
            l += pow(mat(i, 0), p);
        }
        return pow(l, 1./p);
    }

    // Row-vector p-norm
    template<int m, class memBlockType>
    double norm(Matrix<1, m, memBlockType> &mat, int p = 2)
    {
        double l = 0;
        for (int i = 0; i < mat.cols; ++i)
        {
            l += pow(mat(0, i), p);
        }
        return pow(l, 1./p);
    }

    // Frobenius norm
    template<int n, int m, class memBlockType>
    double f_norm(Matrix<n, m, memBlockType> &mat)
    {
        double l = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                l += pow(mat(i, j), 2);
            }
        }
        return sqrt(l);
    }

    // Generalizes Frobenius norm from 2-norm to p-norm
    template<int n, int m, class memBlockType>
    double entrywise_p_norm(Matrix<n, m, memBlockType> &mat, int p)
    {
        double l = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                l += pow(mat(i, j), p);
            }
        }
        return pow(l, 1./p);
    }

    // TODO: Implement induced matrix p-norm approximation
    // TODO: Implement induced 1-, 2-, and infinity-norm
    // Induced p-norm implementation (Yikes)
    template<class matrixType>
    double norm(matrixType &mat, int p)
    {
        return 0.;
    }
}

#endif //CPP_LINALG_NORM_H
