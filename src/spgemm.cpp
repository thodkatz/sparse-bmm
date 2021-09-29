#include <array>
#include <algorithm>
#include "sparsetools.hpp"
#include "utils.hpp"

CSX csxMul(const CSX& csr, const CSX& csc)
{
    uint32_t rows = csr.pointer.size() - 1;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
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
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

CSX csxMulSTL(const CSX& csr, const CSX& csc)
{
    uint32_t rows = csr.pointer.size() - 1;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
    mulResult.pointer.push_back(0);

    uint32_t nnz        = 0;
    uint32_t idStartCol = 0;
    uint32_t idEndCol   = 0;
    uint32_t idStartRow = 0;
    uint32_t idEndRow   = 0;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < rows; j++) {
            idStartCol = csr.pointer[i];
            idEndCol   = csr.pointer[i + 1];
            if (idStartCol == idEndCol)
                break;
            idStartRow = csc.pointer[j];
            idEndRow   = csc.pointer[j + 1];

            // auto result = std::find_first_of(csr.indices.begin() + idStartCol, csr.indices.begin() + idEndCol, csc.indices.begin() + idStartRow, csc.indices.begin() + idEndRow);
            // if (result != csr.indices.begin() + idEndCol) {
            //     mulResult.indices.push_back(j);
            //     nnz++;
            // }

            bool isCommon = hasCommon(csr.indices.begin() + idStartCol, csr.indices.begin() + idEndCol, csc.indices.begin() + idStartRow, csc.indices.begin() + idEndRow);
            if (isCommon) {
                mulResult.indices.push_back(j);
                nnz++;
            }

            // std::vector<uint32_t> result;
            // std::set_intersection(csr.indices.begin() + idStartCol, csr.indices.begin() + idEndCol, csc.indices.begin() + idStartRow, csc.indices.begin() + idEndRow, std::back_inserter(result));
            // if (result.size() != 0) {
            //     mulResult.indices.push_back(j);
            //     nnz++;
            // }

            // bool result = has_common_elements(csr.indices.begin() + idStartCol, csr.indices.begin() + idEndCol, csc.indices.begin() + idStartRow, csc.indices.begin() + idEndRow);
            // if (result) {
            //      mulResult.indices.push_back(j);
            //      nnz++;
            // }
        }
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

BSXPad csxMulPad(const BSXPad& bcsr, const BSXPad& bcsc)
{
    CSX subCSR;
    CSX subCSC;
}
