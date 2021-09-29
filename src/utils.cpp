#include <fstream>
#include "sparsetools.hpp"

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