#ifndef SPGEMM_HPP
#define SPGEMM_HPP

#include "sparsetools.hpp"

template <class InputIterator1, class InputIterator2, class OutputIterator>
OutputIterator indicesIntersection(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, OutputIterator result)
{
    InputIterator1 first1Start = first1;
    InputIterator2 first2Start = first2;

    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2)
            ++first1;
        else if (*first2 < *first1)
            ++first2;
        else {
            *result = first1 - first1Start;
            ++result;
            *result = first2 - first2Start;
            ++result;
            ++first1;
            ++first2;
        }
    }
    return result;
}

/* -------------------- Non Blocking -------------------- */
CSX csxMul(const CSX& csr, const CSX& csc);

CSX csxMulSTL(const CSX& csr, const CSX& csc);

CSX csxMul(const CSX& csrA, const CSX& cscB, const CSX& csrF);

/* ---------------------- Blocking ---------------------- */
CSX bmmPerBlock(const BSXNoPad& csrA, const BSXNoPad& cscB, const CSX& csrF, uint32_t pointerOffsetA, uint32_t pointerOffsetB, uint32_t blockSize);

BSXNoPad bmmBlock(const MatrixInfo& F, const BSXNoPad& bcsrA, const BSXNoPad& bcscB, const BSXNoPad& bcsrF);

void appendResult(BSXNoPad& result, const CSX& csrResultBlock, uint32_t blockSizeY);

CSX subBlockMul(const BSXNoPad& bcsrA,
                const BSXNoPad& bcscB,
                CSX& csrMask,
                uint32_t blockSizeY,
                uint32_t blockSizeX,
                uint32_t indexBlockStartCol,
                uint32_t indexBlockEndCol,
                uint32_t indexBlockStartRow,
                uint32_t indexBlockEndRow);

CSX updateMask(const CSX& current, const CSX& old);

CSX updateBlockC(const CSX& current, const CSX& old);

CSX fillMaskOnes(uint32_t blockSizeY, uint32_t blockSizeX);

void getBlock(const BSXNoPad& bcsx, CSX& block, uint32_t nnzBlocksPassed, uint32_t blockSizeY);

#endif