#include <algorithm>
#include <vector>
#include <cmath>
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
void coo2csr(const MatrixInfo& arr, const std::vector<uint32_t>& cooRows, const std::vector<uint32_t>& cooCols, CSX& csr)
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
void coo2csc(const MatrixInfo& arr, const std::vector<uint32_t>& cooRows, const std::vector<uint32_t>& cooCols, CSX& csc)
{
    MatrixInfo swapArr;
    swapArr.nRow = arr.nCol;
    swapArr.nRow = arr.nCol;
    swapArr.nnz  = arr.nnz;
    coo2csr(swapArr, cooCols, cooRows, csc);
}

/**
 * Convert csr to its blocked version (padding version)
 *
 * In order to know the order of blocks we keep all of them including the zero ones padding the pointer array with the
 * count of nnz elements
 */
void csr2bcsrPad(MatrixInfo& arr, const CSX& csr, BSXPad& bcsr)
{
    if (arr.nCol / arr.blockSizeX == 0 || arr.nRow / arr.blockSizeY == 0) {
        std::cout << "Block dimensions exceed the matrix dimensions" << std::endl;
        exit(-1);
    }

    arr.numBlockX = (arr.nCol % arr.blockSizeX != 0) ? arr.nCol / arr.blockSizeX + 1 : arr.nCol / arr.blockSizeX;
    arr.numBlockY = (arr.nRow % arr.blockSizeY != 0) ? arr.nRow / arr.blockSizeY + 1 : arr.nRow / arr.blockSizeY;

    // this will overflow if blockSize too small TODO prevent happening
    uint32_t totalBlocks = arr.numBlockX * arr.numBlockY;

    bcsr.pointer.resize(totalBlocks * arr.blockSizeY + 1);
    bcsr.indices.resize(arr.nnz);

    bcsr.pointer[0]       = 0;
    uint32_t nnz          = 0;
    uint32_t currentBlock = 0;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockX = 0; blockX < arr.numBlockX; blockX++) {
            for (uint32_t i = blockY * arr.blockSizeY, blockRow = 0; i < (blockY + 1) * arr.blockSizeY; i++, blockRow++) {
                currentBlock = blockY * arr.numBlockX + blockX;

                // padding when out of bounds
                if (i >= arr.nRow) {
                    bcsr.pointer[currentBlock * arr.blockSizeY + blockRow + 1] = nnz;
                    continue;
                }

                // find the column elements of (blockx,blocky) for blockRow
                for (uint32_t j = csr.pointer[i]; j < csr.pointer[i + 1]; j++) {
                    if (csr.indices[j] >= (blockX + 1) * arr.blockSizeX)
                        break;
                    if (blockX * arr.blockSizeX <= csr.indices[j]) {
                        bcsr.indices[nnz++] = csr.indices[j];
                    }
                }

                bcsr.pointer[currentBlock * arr.blockSizeY + blockRow + 1] = nnz;
            }
        }
    }
}

/**
 * Convert csc to its blocked version
 */
void csc2bcscPad(MatrixInfo& arr, const CSX& csc, BSXPad& bcsc)
{
    MatrixInfo& swapArr = arr;
    swapArr.nRow        = arr.nCol;
    swapArr.nCol        = arr.nRow;
    swapArr.nnz         = arr.nnz;
    swapArr.blockSizeX  = arr.blockSizeY;
    swapArr.blockSizeY  = arr.blockSizeX;
    swapArr.numBlockX   = arr.numBlockY;
    swapArr.numBlockY   = arr.numBlockX;

    csr2bcsrPad(swapArr, csc, bcsc);
}

void bcsr2csrPad(const MatrixInfo& arr, const BSXPad& bcsr, CSX& csr)
{
    csr.pointer.resize(arr.nRow + 1);
    csr.indices.resize(arr.nnz);
    csr.pointer[0] = 0;

    uint32_t nnzCount   = 0;
    uint32_t countRow   = 0;
    uint32_t idBlock    = 0;
    uint32_t startBlock = 0;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockRow = 0; blockRow < arr.blockSizeY && countRow < arr.nRow; blockRow++) {
            for (uint32_t blockX = 0; blockX < arr.numBlockX; blockX++) {
                idBlock    = blockX + blockY * arr.numBlockX;
                startBlock = idBlock * arr.blockSizeY;
                for (uint32_t idCol = bcsr.pointer[startBlock + blockRow]; idCol < bcsr.pointer[startBlock + blockRow + 1]; idCol++) {
                    csr.indices[nnzCount++] = bcsr.indices[idCol];
                }
            }
            csr.pointer[++countRow] = nnzCount;
        }
    }
}

void bcsc2cscPad(const MatrixInfo& arr, const BSXPad& bcsc, CSX& csc)
{
    MatrixInfo swapArr;
    swapArr.nRow       = arr.nCol;
    swapArr.nCol       = arr.nRow;
    swapArr.nnz        = arr.nnz;
    swapArr.blockSizeX = arr.blockSizeY;
    swapArr.blockSizeY = arr.blockSizeX;
    swapArr.numBlockX  = arr.numBlockY;
    swapArr.numBlockY  = arr.numBlockX;
    bcsr2csrPad(swapArr, bcsc, csc);
}

/**
 * Convert csr to blocking csr (no padding version)
 *
 * In order to not include the zero blocks, we need to know the order of the blocks and how many of the non zero ones
 * exist per blockX and blockY. So we will add 2 new arrays to the csr structure. One with the size of numBlocksY that
 * will count the number of blockX per blockY and another one with the ids of the blocks. The last 2 arrays will be the
 * pointer and the nnz of each block similar to a simple csr structure.
 */
void csr2bcsrNoPad(MatrixInfo& arr, const CSX& csr, BSXNoPad& bcsr)
{
    if (arr.nCol / arr.blockSizeX == 0 || arr.nRow / arr.blockSizeY == 0) {
        std::cout << "Block dimensions exceed the matrix dimensions" << std::endl;
        exit(-1);
    }

    arr.numBlockX = (arr.nCol % arr.blockSizeX != 0) ? arr.nCol / arr.blockSizeX + 1 : arr.nCol / arr.blockSizeX;
    arr.numBlockY = (arr.nRow % arr.blockSizeY != 0) ? arr.nRow / arr.blockSizeY + 1 : arr.nRow / arr.blockSizeY;

    bcsr.indices.resize(arr.nnz);
    bcsr.pointer.push_back(0);
    bcsr.blockPointer.resize(arr.numBlockY + 1);

    bcsr.blockPointer[0]   = 0;
    uint32_t nnz           = 0;
    uint32_t currentBlock  = 0;
    uint32_t emptyRow      = 0;
    uint32_t nonzeroBlocks = 0;
    bool isFound           = false;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockX = 0; blockX < arr.numBlockX; blockX++) {
            emptyRow = 0;
            for (uint32_t i = blockY * arr.blockSizeY, blockRow = 0; i < (blockY + 1) * arr.blockSizeY; i++, blockRow++) {
                currentBlock = blockY * arr.numBlockX + blockX;

                // padding when out of bounds
                if (i >= arr.nRow) {
                    bcsr.pointer.push_back(nnz);
                    emptyRow++;
                    std::cout << "\nEmpty " << emptyRow;
                    std::cout << "\nOut of bounds, Block Row" << blockRow;
                    printVector(bcsr.pointer, " ");
                    continue;
                }

                // find the column elements of (blockX,blockY) for blockRow
                isFound = false;
                for (uint32_t j = csr.pointer[i]; j < csr.pointer[i + 1]; j++) {
                    if (csr.indices[j] >= (blockX + 1) * arr.blockSizeX) {
                        break;
                    }
                    if (blockX * arr.blockSizeX <= csr.indices[j]) {
                        bcsr.indices[nnz++] = csr.indices[j];
                        isFound             = true;
                    }
                }
                if (!isFound)
                    emptyRow++;

                std::cout << "\nEmpty " << emptyRow << std::endl;

                bcsr.pointer.push_back(nnz);
                std::cout << "Block Row" << blockRow;
                printVector(bcsr.pointer, " ");
            } // filled (blockX,blockY)

            // empty block detected, remove the padding
            if (emptyRow == arr.blockSizeY) {
                std::cout << "Before";
                printVector(bcsr.pointer, " ");
                bcsr.pointer.erase(bcsr.pointer.end() - arr.blockSizeY, bcsr.pointer.end());
                std::cout << "After";
                printVector(bcsr.pointer, " ");
            }
            else {
                bcsr.idBlock.push_back(currentBlock);
                nonzeroBlocks++;
            }
        }
        bcsr.blockPointer[blockY + 1] = nonzeroBlocks;
    }
}

void csc2bcscNoPad(MatrixInfo& arr, const CSX& csc, BSXNoPad& bcsc)
{
    MatrixInfo& swapArr = arr;
    swapArr.nRow        = arr.nCol;
    swapArr.nCol        = arr.nRow;
    swapArr.nnz         = arr.nnz;
    swapArr.blockSizeX  = arr.blockSizeY;
    swapArr.blockSizeY  = arr.blockSizeX;
    swapArr.numBlockX   = arr.numBlockY;
    swapArr.numBlockY   = arr.numBlockX;
    csr2bcsrNoPad(swapArr, csc, bcsc);
}

void bcsr2csrNoPad(const MatrixInfo& arr, const BSXNoPad& bcsrNoPad, CSX& csr)
{
    csr.pointer.resize(arr.nRow + 1);
    csr.indices.resize(arr.nnz);
    csr.pointer[0] = 0;

    uint32_t nnzCount   = 0;
    uint32_t countRow   = 0;
    uint32_t startBlock = 0;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockRow = 0; blockRow < arr.blockSizeY && countRow < arr.nRow; blockRow++) {
            for (uint32_t blockX = bcsrNoPad.blockPointer[blockY]; blockX < bcsrNoPad.blockPointer[blockY+1]; blockX++) {
                startBlock = blockX * arr.blockSizeY;
                for (uint32_t idCol = bcsrNoPad.pointer[startBlock + blockRow]; idCol < bcsrNoPad.pointer[startBlock + blockRow + 1]; idCol++) {
                    csr.indices[nnzCount++] = bcsrNoPad.indices[idCol];
                }
            }
            csr.pointer[++countRow] = nnzCount;
        }
    }
}

void bcsc2cscNoPad(const MatrixInfo& arr, const BSXNoPad& bcscNoPad, CSX& csc)
{
    MatrixInfo swapArr;
    swapArr.nRow       = arr.nCol;
    swapArr.nCol       = arr.nRow;
    swapArr.nnz        = arr.nnz;
    swapArr.blockSizeX = arr.blockSizeY;
    swapArr.blockSizeY = arr.blockSizeX;
    swapArr.numBlockX  = arr.numBlockY;
    swapArr.numBlockY  = arr.numBlockX;
    bcsr2csrNoPad(swapArr, bcscNoPad, csc);
}

/**
 * Convert bcsr or bcsc to its non-blocked version
 *
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