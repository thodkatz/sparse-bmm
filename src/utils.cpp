#include <fstream>
#include <iostream>
#include <algorithm>
#include "utils.hpp"
#include "sparsetools.hpp"

void toDense(CSX csx, uint32_t rows, uint32_t cols, sparseType type, uint32_t pointerOffset, uint32_t indicesOffset)
{
    for(auto& i : csx.pointer) {
        i-=pointerOffset;
    }

    for(auto& i : csx.indices) {
        i-=indicesOffset;
    }

    if (type == sparseType::CSC) {
        std::swap(rows, cols);
    }

    std::vector<std::vector<uint32_t>> dense(rows, std::vector<uint32_t>(cols, 0));
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = csx.pointer[i]; j < csx.pointer[i + 1]; j++) {
            dense[i][csx.indices[j]] = 1;
        }
    }

    std::cout << std::endl;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            uint32_t value = (type == sparseType::CSR) ? dense[i][j] : dense[j][i];
            std::cout << value;
            std::string delimiter = (j != cols - 1) ? " " : "\n";
            std::cout << delimiter;
        }
    }
}

void csxWriteFile(CSX& csx, std::string filename)
{
    uint32_t n = csx.pointer.size()-1;
    uint32_t nnz = csx.indices.size();
    std::ofstream csxMulFile;
    csxMulFile.open(filename);
    for (uint32_t i = 0; i < n; i++) {
        csxMulFile << csx.pointer[i];
        if (i != n-1) {
            csxMulFile << ",";
        }
    }
    csxMulFile << "\n";
    for (uint32_t i = 0; i < nnz; i++) {
        csxMulFile << csx.indices[i];
        if (i != nnz - 1) {
            csxMulFile << ",";
        }
    }
    csxMulFile.close();
}