#ifndef CPP_LINALG_MEM_BLOCKS_H
#define CPP_LINALG_MEM_BLOCKS_H

namespace cla
{
    // TODO: Formal testing
    template<int n, int m, class elementType>
    struct DenseMatrixBlock
    {
        typedef elementType element_type;
        mutable element_type mem_block[n][m];

        inline element_type& operator()(int i, int j) const
        {
            return mem_block[i][j];
        }
    };

    template<class memBlockType>
    struct ReferenceBlock
    {
        typedef typename memBlockType::element_type element_type;
        const memBlockType& referenced;
        int row_offset, col_offset;

        ReferenceBlock(const memBlockType& block, int start_row=0, int start_col=0) : referenced(block),
                                                                                      row_offset(start_row),
                                                                                      col_offset(start_col){};
        explicit ReferenceBlock(const ReferenceBlock<memBlockType>& copied) : referenced(copied.referenced),
                                                                              row_offset(copied.row_offset),
                                                                              col_offset(copied.col_offset){};

        inline element_type& operator()(int i, int j) const
        {
            return referenced(i+row_offset, j+col_offset);
        }
    };

    // TODO: Test
    // TODO: Test whether getNumRows/Cols work, test whether matrix arithmetic work properly
    // TODO: Separate from memblocks?
    template<class memBlockType>
    struct TransposeBlock
    {
        typedef typename memBlockType::element_type element_type;
        const memBlockType original;


        explicit TransposeBlock(const memBlockType& block) : original(block){};
        TransposeBlock(const TransposeBlock<memBlockType>& copied) : original(copied.original){};

        inline element_type& operator()(int i, int j) const
        {
            return original(j, i);
        }
    };

    // TODO: Test
    template <int n, int m, class elementType>
    struct Zero
    {
        static elementType elem = 0;

        inline elementType& operator ()(int i, int j) const
        {
            return elem;
        }
    };

    // TODO: Make memory block layout types for the following types of matrices
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
