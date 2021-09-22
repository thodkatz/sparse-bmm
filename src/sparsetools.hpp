#pragma once

#include <algorithm>
#include "mmio.h"

/*
 * Source:
 * https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/coo.h
 *
 * Compute B = A for COO matrix A, CSR matrix B
 *
 * Input Arguments:
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 * Output Arguments:
 *   I Bp  - row pointer
 *   I Bj  - column indices
 *
 * Note:
 *   Output arrays Bp, Bj, and Bx must be preallocated
 *
 * Note:
 *   Matrix A assumed to be boolean, thus no values array needed.
 *
 * Note:
 *   Input:  row and column indices *are not* assumed to be ordered
 *
 *   Note: duplicate entries are carried over to the CSR represention
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 *
 */
template <class I>
void cooToCsr(const I n_row,
              const I n_col,
              const I nnz,
              const I Ai[],
              const I Aj[],
              I Bp[],
              I Bj[])
{
    // compute number of non-zero entries per row of A
    std::fill(Bp, Bp + n_row, 0);

    for (I n = 0; n < nnz; n++) {
        Bp[Ai[n]]++;
    }

    // cumsum the nnz per row to get Bp[]
    for (I i = 0, cumsum = 0; i < n_row; i++) {
        I temp = Bp[i];
        Bp[i]  = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz;

    // write Aj,Ax into Bj,Bx
    for (I n = 0; n < nnz; n++) {
        I row  = Ai[n];
        I dest = Bp[row];

        Bj[dest] = Aj[n];

        Bp[row]++;
    }

    for (I i = 0, last = 0; i <= n_row; i++) {
        I temp = Bp[i];
        Bp[i]  = last;
        last   = temp;
    }

    // now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/**
 * Source:
 * https://github.com/scipy/scipy/blob/3b36a57/scipy/sparse/sparsetools/coo.h#L34
 */
template <class I>
void cooToCsc(const I n_row,
              const I n_col,
              const I nnz,
              const I Ai[],
              const I Aj[],
              I Bp[],
              I Bi[])
{
    cooToCsr<I>(n_col, n_row, nnz, Aj, Ai, Bp, Bi);
}

/**
 * Source:
 * https://math.nist.gov/MatrixMarket/mmio-c.html
 */
template <class I>
void mm2coo(int argc,
            char* argv[],
            std::vector<I>& rows,
            std::vector<I>& cols,
            I& nnz,
            I& n)
{
    MM_typecode matcode;
    FILE* f;
    I r, c; // MxN dimensions (square matrix M=N)
    // double *val; // dont need this. Our matrices are binary 1 or zero

    // expecting a filename to read (./main <filename>)
    if (argc < 2) {
        printf("Missed command line arguements\n");
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else {
        if ((f = fopen(argv[1], "r")) == NULL) {
            printf("Can't open file\n");
            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    // what MM formats do you support?
    if (!(mm_is_matrix(matcode) && mm_is_coordinate(matcode) &&
          mm_is_pattern(matcode))) {
        printf("Sorry, this application does not support ");
        printf("Matrix Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((mm_read_mtx_crd_size(f, &r, &c, &nnz)) != 0)
        exit(1);
    n = r;

    rows.resize(nnz);
    cols.resize(nnz);

    I x, y = 0;
    for (I i = 0; i < nnz; i++) {
        fscanf(f, "%u %u\n", &x, &y);
        rows[i] = x;
        cols[i] = y;
        /* adjust from 1-based to 0-based */
        rows[i]--; 
        cols[i]--;
    }

    printf("Success, MM format is converted to COO\n");

    if (f != stdin) {
        fclose(f);
    }
}