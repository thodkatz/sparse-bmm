#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include "utils.hpp"

/**
 * Source:
 * https://math.nist.gov/MatrixMarket/mmio-c.html
 */
void mm2coo(char* argv, std::vector<uint32_t>& cooRows, std::vector<uint32_t>& cooCols, MatrixInfo& arr)
{
    std::ifstream f(argv);
    std::string line;
    if (f.is_open()) {
        while (f.peek() == '%') {
            getline(f, line);
        }
        f >> arr.nRow >> arr.nCol >> arr.nnz;

        cooRows.resize(arr.nnz);
        cooCols.resize(arr.nnz);

        for (uint32_t i = 0; i < arr.nnz; i++) {
            f >> cooRows[i] >> cooCols[i];
            cooRows[i]--;
            cooCols[i]--;
        }
        f.close();
    }
    else {
        std::cout << "Unable to open file" << std::endl;
        exit(-1);
    }
}

/*
 * Source:
 * https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/coo.h
 */
void coo2csr(const MatrixInfo& arr, const std::vector<uint32_t>& cooRows, const std::vector<uint32_t>& cooCols, CSX& csr)
{
    // compute number of non-zero entries per row of A
    csr.pointer.resize(arr.nRow + 1);
    csr.indices.resize(arr.nnz);

    for (uint32_t i = 0; i < arr.nnz; i++) {
        csr.pointer[cooRows[i]]++;
    }

    // cumsum the arr.nnz per row to get rowPointer[]
    for (uint32_t i = 0, cumsum = 0; i < arr.nRow; i++) {
        uint32_t temp  = csr.pointer[i];
        csr.pointer[i] = cumsum;
        cumsum += temp;
    }
    csr.pointer[arr.nRow] = arr.nnz;

    // write cooCols,Ax into colIndices,Bx
    for (uint32_t i = 0; i < arr.nnz; i++) {
        uint32_t row  = cooRows[i];
        uint32_t dest = csr.pointer[row];

        csr.indices[dest] = cooCols[i];
        csr.pointer[row]++;
    }

    for (uint32_t i = 0, last = 0; i <= arr.nRow; i++) {
        uint32_t temp  = csr.pointer[i];
        csr.pointer[i] = last;
        last           = temp;
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
 * Convert csr to blocking csr
 *
 * In order to not include the zero blocks, we need to know the order of the blocks and how many of the non zero ones
 * exist per blockX and blockY. So we will add 2 new arrays to the csr structure. One with the size of numBlocksY that
 * will count the number of blockX per blockY and another one with the ids of the blocks. The last 2 arrays will be the
 * pointer and the nnz of each block similar to a simple csr structure.
 */
void csr2bcsr(MatrixInfo& arr, const CSX& csr, BSX& bcsr)
{
    if (arr.blockSizeX == 0 || arr.blockSizeY == 0) {
        std::cout << "Block size is set to 0" << std::endl;
        std::cout << "Pick default block size\n";
        arr.blockSizeX = arr.nCol / 2;
        arr.blockSizeY = arr.nRow / 2;
        exit(-1);
    }

    if (arr.nCol / arr.blockSizeX == 0 || arr.nRow / arr.blockSizeY == 0) {
        std::cout << "Block dimensions exceed the matrix dimensions" << std::endl;
        exit(-1);
    }

    arr.numBlockX = (arr.nCol % arr.blockSizeX != 0) ? arr.nCol / arr.blockSizeX + 1 : arr.nCol / arr.blockSizeX;
    arr.numBlockY = (arr.nRow % arr.blockSizeY != 0) ? arr.nRow / arr.blockSizeY + 1 : arr.nRow / arr.blockSizeY;

    bcsr.indices.resize(arr.nnz);
    bcsr.pointer.push_back(0);
    bcsr.blockPointer.resize(arr.numBlockY + 1);

    bcsr.blockPointer[0]    = 0;
    uint32_t nnz            = 0;
    uint32_t emptyRow       = 0;
    uint32_t nonzeroBlocks  = 0;
    uint32_t blockColOffset = 0;
    bool isFound            = false;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockX = 0; blockX < arr.numBlockX; blockX++) {
            emptyRow = 0;
            for (uint32_t i = blockY * arr.blockSizeY, blockRow = 0; i < (blockY + 1) * arr.blockSizeY; i++, blockRow++) {
                blockColOffset = blockX * arr.blockSizeX;

                // padding when out of bounds
                if (i >= arr.nRow) {
                    bcsr.pointer.push_back(nnz);
                    emptyRow++;
                    continue;
                }

                // find the column elements of (blockX,blockY) for blockRow
                isFound = false;
                for (uint32_t j = csr.pointer[i]; j < csr.pointer[i + 1]; j++) {
                    if (csr.indices[j] >= (blockX + 1) * arr.blockSizeX) {
                        break;
                    }
                    if (blockX * arr.blockSizeX <= csr.indices[j]) {
                        bcsr.indices[nnz++] = csr.indices[j] - blockColOffset;
                        isFound             = true;
                    }
                }
                if (!isFound)
                    emptyRow++;

                bcsr.pointer.push_back(nnz);
            } // filled (blockX,blockY)

            // empty block detected, remove the padding
            if (emptyRow == arr.blockSizeY) {
                bcsr.pointer.erase(bcsr.pointer.end() - arr.blockSizeY, bcsr.pointer.end());
            }
            else {
                bcsr.idBlock.push_back(blockX);
                nonzeroBlocks++;
            }
        }
        bcsr.blockPointer[blockY + 1] = nonzeroBlocks;
    }
}

void csc2bcsc(MatrixInfo& arr, const CSX& csc, BSX& bcsc)
{
    MatrixInfo& swapArr = arr;
    std::swap(swapArr.nCol, swapArr.nRow);
    csr2bcsr(swapArr, csc, bcsc);

    // re-adjust dimensions
    std::swap(swapArr.nCol, swapArr.nRow);
    std::swap(swapArr.blockSizeX, swapArr.blockSizeY);
    std::swap(swapArr.numBlockX, swapArr.numBlockY);
}

void bcsr2csr(const MatrixInfo& arr, const BSX& bcsr, CSX& csr)
{
    csr.pointer.resize(arr.nRow + 1);
    csr.indices.resize(arr.nnz);
    csr.pointer[0] = 0;

    uint32_t nnzCount       = 0;
    uint32_t countRow       = 0;
    uint32_t startBlock     = 0;
    uint32_t blockColOffset = 0;
    for (uint32_t blockY = 0; blockY < arr.numBlockY; blockY++) {
        for (uint32_t blockRow = 0; blockRow < arr.blockSizeY && countRow < arr.nRow; blockRow++) {
            for (uint32_t blockX = bcsr.blockPointer[blockY]; blockX < bcsr.blockPointer[blockY + 1]; blockX++) {
                startBlock     = blockX * arr.blockSizeY;
                blockColOffset = bcsr.idBlock[blockX] * arr.blockSizeX;
                for (uint32_t idCol = bcsr.pointer[startBlock + blockRow]; idCol < bcsr.pointer[startBlock + blockRow + 1]; idCol++) {
                    csr.indices[nnzCount++] = bcsr.indices[idCol] + blockColOffset;
                }
            }
            csr.pointer[++countRow] = nnzCount;
        }
    }
}

void bcsc2csc(const MatrixInfo& arr, const BSX& bcsc, CSX& csc)
{
    MatrixInfo swapArray = arr;
    std::swap(swapArray.nCol, swapArray.nRow);
    std::swap(swapArray.blockSizeX, swapArray.blockSizeY);
    std::swap(swapArray.numBlockX, swapArray.numBlockY);
    bcsr2csr(swapArray, bcsc, csc);
}

/**
 * Convert bcsr or bcsc to its non-blocked version
 *
 */
void mm2csr(char* argv, CSX& csr, MatrixInfo& arr)
{
    std::vector<uint32_t> cooRows;
    std::vector<uint32_t> cooCols;
    mm2coo(argv, cooRows, cooCols, arr);
    // printVector(cooRows, " ");
    // printVector(cooCols, " ");
    coo2csr(arr, cooRows, cooCols, csr);
}

void mm2csc(char* argv, CSX& csc, MatrixInfo& arr)
{
    std::vector<uint32_t> cooRows;
    std::vector<uint32_t> cooCols;
    mm2coo(argv, cooRows, cooCols, arr);
    coo2csc(arr, cooRows, cooCols, csc);
}

void readInput(char* argv[], MatrixInfo& A, MatrixInfo& B, MatrixInfo& F, CSX& csrA, CSX& cscB, CSX& csrF)
{
    {
        Timer time("Time converting input to CSX: \n");
        mm2csr(argv[1], csrA, A);
        mm2csc(argv[2], cscB, B);
        mm2csr(argv[3], csrF, F);
    }

    std::cout << "A: Rows: " << A.nRow << ", Cols: " << A.nCol << ", nnz: " << A.nnz << std::endl
              << "B: Rows: " << B.nRow << ", Cols: " << B.nCol << ", nnz: " << B.nnz << std::endl
              << "F: Rows: " << F.nRow << ", Cols: " << F.nCol << ", nnz: " << F.nnz << std::endl;

    if (A.nCol != B.nRow) {
        std::cout << "Multiplication dimensions mismatch" << std::endl;
        exit(-1);
    }
    if (F.nRow != A.nRow || F.nCol != B.nCol) {
        std::cout << "Masked matrix dimensions mismatch" << std::endl;
        exit(-1);
    }
}