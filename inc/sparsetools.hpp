#ifndef SPARSETOOLS_HPP
#define SPARSETOOLS_HPP

#include <cstdint>
#include <vector>

typedef struct MatrixInfo {
    uint32_t nRow       = 0;
    uint32_t nCol       = 0;
    uint32_t nnz        = 0;
    uint32_t blockSizeX = 0;
    uint32_t blockSizeY = 0;
    uint32_t numBlockX  = 0;
    uint32_t numBlockY  = 0;
} MatrixInfo;

/**
 * \param n Length of the pointer matrix (rows+1 if csr or cols+1 if csc)
 */
class CSX {
  public:
    uint32_t nnz;
    uint32_t n;
    std::vector<uint32_t> pointer;
    std::vector<uint32_t> indices;

    CSX(uint32_t a, uint32_t b) : nnz{a}, n{b} {}
    CSX() : nnz{0}, n{0} {}
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

void csr2bcsr(MatrixInfo& arr, const CSX& csr, CSX& bcsr);

void csc2bcsc(const MatrixInfo& arr, const CSX& csc, CSX& bcsc);

void mm2csr(char argc, char* argv, CSX& csr, MatrixInfo& arr);

void mm2csc(char argc, char* argv, CSX& csc, MatrixInfo& arr);

#endif