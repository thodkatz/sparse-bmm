#ifndef SPARSETOOLS_HPP
#define SPARSETOOLS_HPP

#include <cstdint>
#include <vector>

class MatrixInfo {
  public:
    uint32_t nRow;
    uint32_t nCol;
    uint32_t nnz;
    uint32_t blockSizeX;
    uint32_t blockSizeY;
    uint32_t numBlockX;
    uint32_t numBlockY;

    MatrixInfo() : nRow(0), nCol(0), nnz(0), blockSizeX(0), blockSizeY(0), numBlockX(0), numBlockY(0) {}
    MatrixInfo(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e, uint32_t f, uint32_t g)
        : nRow(a), nCol(b), nnz(c), blockSizeX(d), blockSizeY(e), numBlockX(f), numBlockY(g)
    {
    }
};

/**
 */
class CSX {
  public:
    std::vector<uint32_t> pointer;
    std::vector<uint32_t> indices;
};

/**
 */
class BSXPad : public CSX {
};

/**
 */
class BSXNoPad : public CSX {
  public:
    std::vector<uint32_t> idBlock;
    std::vector<uint32_t> blockPointer;
};

void mm2coo(int argc, char* argv, std::vector<uint32_t>& cooRows, std::vector<uint32_t>& cooCols, MatrixInfo& arr);

void coo2csr(const MatrixInfo& arr,
             const std::vector<uint32_t>& cooRows,
             const std::vector<uint32_t>& cooCols,
             CSX& csr);

void coo2csc(const MatrixInfo& arr,
             const std::vector<uint32_t>& cooRows,
             const std::vector<uint32_t>& cooCols,
             CSX& csc);

void csr2bcsrPad(MatrixInfo& arr, const CSX& csr, BSXPad& bcsr);

void csc2bcscPad(MatrixInfo& arr, const CSX& csc, BSXPad& bcsc);

void bcsr2csrPad(const MatrixInfo& arr, const BSXPad& bcsr, CSX& csr);

void bcsc2cscPad(const MatrixInfo& arr, const BSXPad& bcsc, CSX& csc);

void csr2bcsrNoPad(MatrixInfo& arr, const CSX& csr, BSXNoPad& bcsr);

void csc2bcscNoPad(MatrixInfo& arr, const CSX& csc, BSXNoPad& bcsc);

void bcsr2csrNoPad(const MatrixInfo& arr, const BSXNoPad& bcsr, CSX& csr);

void bcsc2cscNoPad(const MatrixInfo& arr, const BSXNoPad& bcsc, CSX& csc);

void mm2csr(char argc, char* argv, CSX& csr, MatrixInfo& arr);

void mm2csc(char argc, char* argv, CSX& csc, MatrixInfo& arr);

#endif