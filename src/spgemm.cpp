#include "sparsetools.hpp"

CSX csxMul(CSX& csr, CSX& csc)
{
    CSX mulResult;
    mulResult.n = csr.n;
    uint32_t n = csr.n-1;

    mulResult.pointer.push_back(0);
    for (uint32_t i = 0, idxCol = 0, idxRow = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            idxCol = csr.pointer[i];
            idxRow = csc.pointer[j];
            if (idxCol == csr.pointer[i + 1]) {
                break;
            }
            // find common elements row ith and column ith (sorted)
            while (idxCol < csr.pointer[i + 1] && idxRow < csc.pointer[j + 1]) {
                if (csr.indices[idxCol] == csc.indices[idxRow]) {
                    mulResult.indices.push_back(j);
                    mulResult.nnz++;
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
        mulResult.pointer.push_back(mulResult.nnz);
    }

    return mulResult;
}
