#include "sparsetools.hpp"
#include <array>

CSX csxMul(const CSX& csr, const CSX& csc)
{
    uint32_t rows = csr.pointer.size()-1;

    CSX mulResult;
    mulResult.pointer.resize(rows+1);
    mulResult.pointer.push_back(0);

    uint32_t nnz = 0;
    for (uint32_t i = 0, idxCol = 0, idxRow = 0; i < rows; i++) {
        for (uint32_t j = 0; j < rows; j++) {
            idxCol = csr.pointer[i];
            idxRow = csc.pointer[j];
            if (idxCol == csr.pointer[i + 1]) {
                break;
            }

            // find common elements row ith and column ith (sorted)
            // two pointers version
            while (idxCol < csr.pointer[i + 1] && idxRow < csc.pointer[j + 1]) {
                if (csr.indices[idxCol] == csc.indices[idxRow]) {
                    mulResult.indices.push_back(j);
                    nnz++;
                    break;
                }
                if (csr.indices[idxCol] > csc.indices[idxRow]) {
                    idxRow++;
                }
                if (csr.indices[idxCol] < csc.indices[idxRow]) {
                    idxCol++;
                }
            }
        }
        mulResult.pointer[i+1] = nnz;
    }

    return mulResult;
}
