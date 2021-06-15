//
// Created by rns on 6/13/21.
//

#ifndef CPP_LINALG_MEM_BLOCKS_H
#define CPP_LINALG_MEM_BLOCKS_H

namespace cla
{
    // TODO: Formal testing
    template<int n, int m, class elementType>
    struct DenseMatrixBlock
    {
        typedef elementType elementType;
        elementType mem_block[n*m];

        inline elementType& operator ()(int i, int j)
        {
            return mem_block[i * m + j];
        }
    };

    // TODO: Test
    // TODO: Test whether getNumRows/Cols work, test whether matrix arithmetic work properly
    // TODO: Separate from memblocks
    template<class MatrixType>
    struct TransposeBlock
    {
        const MatrixType original;

        explicit TransposeBlock(const MatrixType& mat) : original(mat){};
        // TODO: Read about copy constructors
        TransposeBlock(const TransposeBlock<MatrixType>& mat) : original(mat.original){};

        inline elementType& operator ()(int i, int j)
        {
            return original(j, i);
        }
    };

    // TODO: Test
    template <int n, int m, int elementType>
    struct Zero
    {
        static elementType elem = 0;

        inline elementType& operator ()(int i, int j)
        {
            return elem;
        }
    };

    // TODO: Make memory block layout types for the following types of matrices
    // Square
    // Diagonal
    // Upper triangular
    // Lower triangular
    // Symmetric/Hermitian (A = A^T  /  A = A^*)
    // Unitary (For a real matrix: A A^T = A^T A = I, implying that A^T = A^(-1) )
    // Identity matrix (transform to diagonal when scaled)

    // Reference matrix
    // Minor (https://en.wikipedia.org/wiki/Minor_(linear_algebra))
    // Submatrix
    // Slice
    // Sparse matrix
    // Matrix augmentation blocks
    // Homogeneous transform augmentation (when augmenting a 3x3 matrix with 3x1 vector, automatically add [0 0 0 1])
}

#endif //CPP_LINALG_MEM_BLOCKS_H
