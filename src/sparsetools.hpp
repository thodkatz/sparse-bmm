#pragma once

#include <algorithm>
#include "mmio.hpp"
#include "utils.hpp"

typedef struct MatrixInfo {
    uint32_t nRow = 0;
    uint32_t nCol = 0;
    uint32_t nnz  = 0;
} MatrixInfo;

class CSX {
  public:
    uint32_t nnz;
    uint32_t n;
    std::vector<uint32_t> pointer;
    std::vector<uint32_t> indices;

    CSX(uint32_t a, uint32_t b) : nnz{a}, n{b} {}
    CSX() : nnz{0}, n{0} {}
};

/**
 * Source:
 * https://math.nist.gov/MatrixMarket/mmio-c.html
 */
void mm2coo(int argc,
            char* argv,
            std::vector<uint32_t>& rows,
            std::vector<uint32_t>& cols,
            uint32_t& n_row,
            uint32_t& n_col,
            uint32_t& nnz)
{
    MM_typecode matcode;
    FILE* f;

    // expecting a filename to read (./main <filename>)
    if (argc < 2) {
        printf("Missed command line arguements\n");
        fprintf(stderr, "Usage: ./bin [martix-market-filename]\n");
        exit(1);
    }
    else {
        if ((f = fopen(argv, "r")) == NULL) {
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
    if ((mm_read_mtx_crd_size(f, &n_row, &n_col, &nnz)) != 0)
        exit(1);

    rows.resize(nnz);
    cols.resize(nnz);

    uint32_t x, y = 0;
    for (uint32_t i = 0; i < nnz; i++) {
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
void coo2csr(const uint32_t n_row,
             const uint32_t n_col,
             const uint32_t nnz,
             const std::vector<uint32_t>& Ai,
             const std::vector<uint32_t>& Aj,
             std::vector<uint32_t>& Bp,
             std::vector<uint32_t>& Bj)
{
    // compute number of non-zero entries per row of A
    Bp.resize(n_row + 1);
    Bj.resize(nnz);

    for (uint32_t n = 0; n < nnz; n++) {
        Bp[Ai[n]]++;
    }

    // cumsum the nnz per row to get Bp[]
    for (uint32_t i = 0, cumsum = 0; i < n_row; i++) {
        uint32_t temp = Bp[i];
        Bp[i]         = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz;

    // write Aj,Ax into Bj,Bx
    for (uint32_t n = 0; n < nnz; n++) {
        uint32_t row  = Ai[n];
        uint32_t dest = Bp[row];

        Bj[dest] = Aj[n];
        Bp[row]++;
    }

    for (uint32_t i = 0, last = 0; i <= n_row; i++) {
        uint32_t temp = Bp[i];
        Bp[i]         = last;
        last          = temp;
    }

    // now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/**
 * Source:
 * https://github.com/scipy/scipy/blob/3b36a57/scipy/sparse/sparsetools/coo.h#L34
 */
template <class I>
void coo2csc(const I n_row,
             const I n_col,
             const I nnz,
             const std::vector<uint32_t>& Ai,
             const std::vector<uint32_t>& Aj,
             std::vector<uint32_t>& Bi,
             std::vector<uint32_t>& Bp)
{
    coo2csr(n_col, n_row, nnz, Aj, Ai, Bp, Bi);
}

void mm2csr(char argc,
            char* argv,
            std::vector<uint32_t>& csrRow,
            std::vector<uint32_t>& csrCol,
            uint32_t& n_row,
            uint32_t& n_col,
            uint32_t& nnz)
{
    std::vector<uint32_t> rows;
    std::vector<uint32_t> cols;
    mm2coo(argc, argv, rows, cols, n_row, n_col, nnz);
    coo2csr(n_row, n_col, nnz, rows, cols, csrRow, csrCol);
}

void mm2csc(char argc,
            char* argv,
            std::vector<uint32_t>& cscRow,
            std::vector<uint32_t>& cscCol,
            uint32_t& n_row,
            uint32_t& n_col,
            uint32_t& nnz)
{
    std::vector<uint32_t> rows;
    std::vector<uint32_t> cols;
    mm2coo(argc, argv, rows, cols, n_row, n_col, nnz);
    coo2csc(n_row, n_col, nnz, rows, cols, cscRow, cscCol);
}