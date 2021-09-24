#include <algorithm>
#include <vector>
#include "mmio.hpp"
#include "utils.hpp"

/**
 * Source:
 * https://math.nist.gov/MatrixMarket/mmio-c.html
 */
void mm2coo(int argc, char* argv, std::vector<uint32_t>& cooRows, std::vector<uint32_t>& cooCols, MatrixInfo& arr)
{
    uint32_t& nRow = arr.nRow;
    uint32_t& nCol = arr.nCol;
    uint32_t& nnz  = arr.nnz;

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
    if (!(mm_is_matrix(matcode) && mm_is_coordinate(matcode) && mm_is_pattern(matcode))) {
        printf("Sorry, this application does not support ");
        printf("Matrix Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((mm_read_mtx_crd_size(f, &nRow, &nCol, &nnz)) != 0)
        exit(1);

    cooRows.resize(nnz);
    cooCols.resize(nnz);

    uint32_t x, y = 0;
    for (uint32_t i = 0; i < nnz; i++) {
        fscanf(f, "%u %u\n", &x, &y);
        cooRows[i] = x;
        cooCols[i] = y;
        /* adjust from 1-based to 0-based */
        cooRows[i]--;
        cooCols[i]--;
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
void coo2csr(const MatrixInfo& arr,
             const std::vector<uint32_t>& cooRows,
             const std::vector<uint32_t>& cooCols,
             CSX& csr)
{
    const uint32_t& nRow = arr.nRow;
    const uint32_t& nnz  = arr.nnz;

    std::vector<uint32_t>& csrRow = csr.pointer;
    std::vector<uint32_t>& csrCol = csr.indices;

    // compute number of non-zero entries per row of A
    csrRow.resize(nRow + 1);
    csrCol.resize(nnz);

    for (uint32_t i = 0; i < nnz; i++) {
        csrRow[cooRows[i]]++;
    }

    // cumsum the nnz per row to get rowPointer[]
    for (uint32_t i = 0, cumsum = 0; i < nRow; i++) {
        uint32_t temp = csrRow[i];
        csrRow[i]     = cumsum;
        cumsum += temp;
    }
    csrRow[nRow] = nnz;

    // write cooCols,Ax into colIndices,Bx
    for (uint32_t i = 0; i < nnz; i++) {
        uint32_t row  = cooRows[i];
        uint32_t dest = csrRow[row];

        csrCol[dest] = cooCols[i];
        csrRow[row]++;
    }

    for (uint32_t i = 0, last = 0; i <= nRow; i++) {
        uint32_t temp = csrRow[i];
        csrRow[i]     = last;
        last          = temp;
    }

    // now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/**
 * Source:
 * https://github.com/scipy/scipy/blob/3b36a57/scipy/sparse/sparsetools/coo.h#L34
 */
void coo2csc(const MatrixInfo& arr,
             const std::vector<uint32_t>& cooRows,
             const std::vector<uint32_t>& cooCols,
             CSX& csc)
{
    MatrixInfo swapArr;
    swapArr.nRow = arr.nCol;
    swapArr.nRow = arr.nCol;
    swapArr.nnz  = arr.nnz;
    coo2csr(swapArr, cooCols, cooRows, csc);
}

/**
 * Convert csr to its blocked version
 * TODO: Explain the blocking
 */
void csr2bcsr(MatrixInfo& arr, const CSX& csr, CSX& bcsr)
{
    const uint32_t& nRow       = arr.nRow;
    const uint32_t& nCol       = arr.nCol;
    const uint32_t& nnz        = arr.nnz;
    const uint32_t& blockSizeX = arr.blockSizeX;
    const uint32_t& blockSizeY = arr.blockSizeY;
    uint32_t& numBlockX        = arr.numBlockX;
    uint32_t& numBlockY        = arr.numBlockY;

    const std::vector<uint32_t>& csrRow  = csr.pointer;
    const std::vector<uint32_t>& csrCol  = csr.indices;
    std::vector<uint32_t>& bcsrRow = bcsr.pointer;
    std::vector<uint32_t>& bcsrCol = bcsr.indices;

    if (nCol / blockSizeX == 0 || nRow / blockSizeY == 0) {
        std::cout << "Block dimensions exceed the matrix dimensions" << std::endl;
        exit(-1);
    }
    numBlockX            = (nCol % blockSizeX != 0) ? nCol / blockSizeX + 1 : nCol / blockSizeX;
    numBlockY            = (nRow % blockSizeY != 0) ? nRow / blockSizeY + 1 : nRow / blockSizeY;
    uint32_t totalBlocks = numBlockX * numBlockY;

    bcsrRow.resize(totalBlocks * blockSizeY + 1);
    bcsrCol.resize(nnz);

    bcsrRow[0]            = 0;
    uint32_t nnzBlocks    = 0;
    uint32_t currentBlock = 0;
    for (uint32_t blockY = 0; blockY < numBlockY; blockY++) {
        for (uint32_t blockX = 0; blockX < numBlockX; blockX++) {
            for (uint32_t i = blockY * blockSizeY, blockRow = 0; i < (blockY + 1) * blockSizeY; i++, blockRow++) {
                currentBlock = blockY * numBlockX + blockX;

                // padding when out of bounds
                if (i >= nRow) {
                    bcsrRow[currentBlock * blockSizeY + blockRow + 1] = nnzBlocks;
                    continue;
                }

                // find the column elements of (blockx,blocky) for blockRow
                for (uint32_t j = csrRow[i]; j < csrRow[i + 1]; j++) {
                    if (csrCol[j] >= (blockX + 1) * blockSizeX)
                        break;
                    if (blockX * blockSizeX <= csrCol[j]) {
                        bcsrCol[nnzBlocks++] = csrCol[j];
                    }
                }

                bcsrRow[currentBlock * blockSizeY + blockRow + 1] = nnzBlocks;
            }
        }
    }
}

/**
 * Convert csr to its blocked version
 * TODO: Explain the blocking
 */
void csc2bcsc(const MatrixInfo& arr, const CSX& csc, CSX& bcsc)
{
    MatrixInfo swapArr;
    swapArr.nRow       = arr.nCol;
    swapArr.nCol       = arr.nRow;
    swapArr.nnz        = arr.nnz;
    swapArr.blockSizeX = arr.blockSizeY;
    swapArr.blockSizeY = arr.blockSizeX;
    swapArr.numBlockX  = arr.numBlockY;
    swapArr.numBlockY  = arr.numBlockX;

    csr2bcsr(swapArr, csc, bcsc);
}

/**
 * Convert bcsr or bcsc to its non-blocked version
 *
 * Note:
 * - Assumed that the bloking versions was created by \p csx2bcsx()
 */

void mm2csr(char argc, char* argv, CSX& csr, MatrixInfo& arr)
{
    std::vector<uint32_t> cooRows;
    std::vector<uint32_t> cooCols;
    mm2coo(argc, argv, cooRows, cooCols, arr);
    coo2csr(arr, cooRows, cooCols, csr);
}

void mm2csc(char argc, char* argv, CSX& csc, MatrixInfo& arr)
{
    std::vector<uint32_t> cooRows;
    std::vector<uint32_t> cooCols;
    mm2coo(argc, argv, cooRows, cooCols, arr);
    coo2csc(arr, cooRows, cooCols, csc);
}