#include <array>
#include <algorithm>
#include <utility>
#include "spgemm.hpp"
#include "sparsetools.hpp"
#include "utils.hpp"

#ifdef HYBRID
#define OPENMP
#endif

#ifdef OPENMP
#include <cstdlib>
#include <unistd.h>
#include <omp.h>
#endif

CSX bmm(const CSX& csrA, const CSX& cscB)
{
    uint32_t rows = csrA.pointer.size() - 1;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
    mulResult.pointer[0] = 0;

    uint32_t nnz        = 0;
    uint32_t idStartCol = 0;
    uint32_t idEndCol   = 0;
    uint32_t idStartRow = 0;
    uint32_t idEndRow   = 0;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < rows; j++) {
            idStartCol = csrA.pointer[i];
            idEndCol   = csrA.pointer[i + 1];
            if (idStartCol == idEndCol) {
                break;
            }

            idStartRow = cscB.pointer[j];
            idEndRow   = cscB.pointer[j + 1];

            bool isCommon = hasCommon(
                csrA.indices.begin() + idStartCol, csrA.indices.begin() + idEndCol, cscB.indices.begin() + idStartRow, cscB.indices.begin() + idEndRow);
            if (isCommon) {
                mulResult.indices.push_back(j);
                nnz++;
            }
        }
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

CSX bmm(const CSX& csrA, const CSX& cscB, const CSX& csrF)
{
    uint32_t rows = csrA.pointer.size() - 1;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
    mulResult.pointer[0];

    uint32_t nnz        = 0;
    uint32_t idStartCol = 0;
    uint32_t idEndCol   = 0;
    uint32_t idStartRow = 0;
    uint32_t idEndRow   = 0;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t idCol = csrF.pointer[i], j = 0; idCol < csrF.pointer[i + 1]; idCol++) {
            j = csrF.indices[idCol];

            idStartCol = csrA.pointer[i];
            idEndCol   = csrA.pointer[i + 1];
            if (idStartCol == idEndCol) {
                break;
            }

            idStartRow = cscB.pointer[j];
            idEndRow   = cscB.pointer[j + 1];

            bool isCommon = hasCommon(
                csrA.indices.begin() + idStartCol, csrA.indices.begin() + idEndCol, cscB.indices.begin() + idStartRow, cscB.indices.begin() + idEndRow);
            if (isCommon) {
                mulResult.indices.push_back(j);
                nnz++;
            }
        }
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

CSX bmmPerBlock(const BSX& csrA, const BSX& cscB, const CSX& csrF, uint32_t pointerOffsetA, uint32_t pointerOffsetB, uint32_t blockSize)
{
    uint32_t rows = blockSize;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
    mulResult.pointer[0] = 0;

    uint32_t nnz        = 0;
    uint32_t idStartCol = 0;
    uint32_t idEndCol   = 0;
    uint32_t idStartRow = 0;
    uint32_t idEndRow   = 0;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t idCol = csrF.pointer[i], j = 0; idCol < csrF.pointer[i + 1]; idCol++) {
            j = csrF.indices[idCol];

            idStartCol = csrA.pointer[i + pointerOffsetA];
            idEndCol   = csrA.pointer[i + pointerOffsetA + 1];
            if (idStartCol == idEndCol) {
                break;
            }

            idStartRow = cscB.pointer[j + pointerOffsetB];
            idEndRow   = cscB.pointer[j + pointerOffsetB + 1];

            bool isCommon = hasCommon(
                csrA.indices.begin() + idStartCol, csrA.indices.begin() + idEndCol, cscB.indices.begin() + idStartRow, cscB.indices.begin() + idEndRow);
            if (isCommon) {
                mulResult.indices.push_back(j);
                nnz++;
            }
        }
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

/**
 *  Blocking Matrix-Matrix Multiplication C = F .* (A * B), where C-BCSR, A-BCSR, B-BCSC, F-CSR
 *  Note:
 *  The blocking techniques used is the no padding, keep tracking the non zero blocks.
 */
BSX bmmBlock(const MatrixInfo& A, const BSX& bcsrA, const BSX& bcscB, const BSX& bcsrF)
{
#ifdef OPENMP
    if (std::getenv("OMP_NUM_THREADS") == nullptr) {
        std::cout << "Environment variable not found\n";
        std::cout << "Usage: export OMP_NUM_THREADS=<number>\n";
        exit(-1);
    }
    uint32_t numOfThreads = std::stoi(std::getenv("OMP_NUM_THREADS"));
    std::cout << "Number of threads: " << numOfThreads << std::endl;

    std::vector<BSX> results(numOfThreads);
    for_each(results.begin(), results.end(), [](BSX& i) { i.pointer.push_back(0); });
    for_each(results.begin(), results.end(), [](BSX& i) { i.blockPointer.push_back(0); });

    std::vector<uint32_t> nnzBlocks(numOfThreads);
#else
    BSX result;
    result.pointer.push_back(0);
    result.blockPointer.push_back(0);
    uint32_t nnzBlocks = 0;
#endif

#ifdef OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (uint32_t blockY = 0; blockY < A.numBlockY; blockY++) {
        uint32_t indexBlockStartCol = 0;
        uint32_t indexBlockEndCol   = 0;
        uint32_t indexBlockStartRow = 0;
        uint32_t indexBlockEndRow   = 0;

#ifdef OPENMP
        uint32_t tid = omp_get_thread_num();
#endif

        uint32_t idBlockCol = 0;
        for (uint32_t indexBlockX = bcsrF.blockPointer[blockY]; indexBlockX < bcsrF.blockPointer[blockY + 1]; indexBlockX++) {

            idBlockCol = bcsrF.idBlock[indexBlockX];

            indexBlockStartCol = bcsrA.blockPointer[blockY];
            indexBlockEndCol   = bcsrA.blockPointer[blockY + 1];
            if (indexBlockStartCol == indexBlockEndCol) {
                break;
            }

            indexBlockStartRow = bcscB.blockPointer[idBlockCol];
            indexBlockEndRow   = bcscB.blockPointer[idBlockCol + 1];

            CSX csrBlockMask;
            getBlock(bcsrF, csrBlockMask, indexBlockX, A.blockSizeY);

            CSX csrResultBlock =
                subBlockMul(bcsrA, bcscB, csrBlockMask, A.blockSizeY, A.blockSizeX, indexBlockStartCol, indexBlockEndCol, indexBlockStartRow, indexBlockEndRow);

            // check if empty block
            if (csrResultBlock.pointer[A.blockSizeY] != 0) {
#ifdef OPENMP
                results[tid].idBlock.push_back(idBlockCol);
                appendResult(results[tid], csrResultBlock, A.blockSizeY);
                nnzBlocks[tid]++;
#else
                result.idBlock.push_back(idBlockCol);
                appendResult(result, csrResultBlock, A.blockSizeY);
                nnzBlocks++;
#endif
            }
        }
#ifdef OPENMP
        results[tid].blockPointer.push_back(nnzBlocks[tid]);
#else
        result.blockPointer.push_back(nnzBlocks);
#endif
    }

#ifdef OPENMP
    // concatenate BSX per thread to a single one
    BSX result = concatBSX(results);
#endif

    return result;
}

/**
 *  Take a vector of BSX and return one unified. Assume that vector is ordered with respect blockY to handle offset requirements for pointers arrays.
 */
BSX concatBSX(std::vector<BSX>& result)
{
    BSX ret;
    ret.blockPointer.push_back(0);
    ret.pointer.push_back(0);

    uint32_t offsetPointer      = 0;
    uint32_t offsetBlockPointer = 0;
    uint32_t lenPointer         = 0;
    uint32_t lenBlockPointer    = 0;

    for (uint32_t i = 0; i < result.size(); i++) {
        if (i > 0) {
            lenBlockPointer    = result[i - 1].blockPointer.size();
            offsetBlockPointer = result[i - 1].blockPointer[lenBlockPointer - 1];
            for_each(result[i].blockPointer.begin(), result[i].blockPointer.end(), [offsetBlockPointer](uint32_t& i) { i += offsetBlockPointer; });

            lenPointer    = result[i - 1].pointer.size();
            offsetPointer = result[i - 1].pointer[lenPointer - 1];
            for_each(result[i].pointer.begin(), result[i].pointer.end(), [offsetPointer](uint32_t& i) { i += offsetPointer; });
        }
        ret.blockPointer.insert(
            ret.blockPointer.end(), make_move_iterator(result[i].blockPointer.begin() + 1), make_move_iterator(result[i].blockPointer.end()));
        ret.idBlock.insert(ret.idBlock.end(), make_move_iterator(result[i].idBlock.begin()), make_move_iterator(result[i].idBlock.end()));
        ret.pointer.insert(ret.pointer.end(), make_move_iterator(result[i].pointer.begin() + 1), make_move_iterator(result[i].pointer.end()));
        ret.indices.insert(ret.indices.end(), make_move_iterator(result[i].indices.begin()), make_move_iterator(result[i].indices.end()));
    }

    return ret;
}

/**
 * Given a BSX format, return an isolated block in CSX format
 */
void getBlock(const BSX& bcsx, CSX& block, uint32_t nnzBlocksPassed, uint32_t blockSizeY)
{
    uint32_t startPointer = nnzBlocksPassed * blockSizeY;
    uint32_t endPointer   = startPointer + blockSizeY;
    uint32_t startIndices = bcsx.pointer[startPointer];
    uint32_t endIndices   = bcsx.pointer[endPointer];

    std::copy(bcsx.pointer.begin() + startPointer, bcsx.pointer.begin() + endPointer + 1, std::back_inserter(block.pointer));
    std::copy(bcsx.indices.begin() + startIndices, bcsx.indices.begin() + endIndices, std::back_inserter(block.indices));

    uint32_t pointerOffsetF = bcsx.pointer[nnzBlocksPassed * blockSizeY];
    for_each(block.pointer.begin(), block.pointer.end(), [pointerOffsetF](uint32_t& i) { i -= pointerOffsetF; });
}

/**
 * Append a CSX format to BSX
 */
void appendResult(BSX& result, const CSX& csrResultBlock, uint32_t blockSizeY)
{
    uint32_t offset = result.pointer[result.pointer.size() - 1];

    // append pointer and indices content
    result.pointer.insert(result.pointer.end(), csrResultBlock.pointer.begin() + 1, csrResultBlock.pointer.end());
    result.indices.insert(result.indices.end(), csrResultBlock.indices.begin(), csrResultBlock.indices.end());

    std::for_each(result.pointer.end() - blockSizeY, result.pointer.end(), [offset](uint32_t& i) { i += offset; });
}

/**
 * Take the blocks of a row-blocked and a column-blocked and perform the multiplication
 */
CSX subBlockMul(const BSX& bcsrA,
                const BSX& bcscB,
                CSX& csrMask,
                uint32_t blockSizeY,
                uint32_t blockSizeX,
                uint32_t indexBlockStartCol,
                uint32_t indexBlockEndCol,
                uint32_t indexBlockStartRow,
                uint32_t indexBlockEndRow)
{
    // find the indices of the common elements
    // pairs are stored in a contiguous manner
    // e.g. (1,2),(3,3)
    std::vector<uint32_t> indicesOfCommon;
    indicesIntersection(bcsrA.idBlock.begin() + indexBlockStartCol,
                        bcsrA.idBlock.begin() + indexBlockEndCol,
                        bcscB.idBlock.begin() + indexBlockStartRow,
                        bcscB.idBlock.begin() + indexBlockEndRow,
                        std::back_inserter(indicesOfCommon));

    CSX csrBlockCold;
    CSX csrBlockCnew;
    std::fill_n(std::back_inserter(csrBlockCnew.pointer), blockSizeY + 1, 0);

    uint32_t pointerOffsetA = 0;
    uint32_t pointerOffsetB = 0;

    // for every pair of aligned blocks perform the bmm
    for (uint32_t i = 0, first = 0, second = 0; i < indicesOfCommon.size(); i += 2) {
        first          = indicesOfCommon[i];
        second         = indicesOfCommon[i + 1];
        pointerOffsetA = (first + indexBlockStartCol) * blockSizeY;
        pointerOffsetB = (second + indexBlockStartRow) * blockSizeY;

        if (i == 0) {
            csrBlockCnew = bmmPerBlock(bcsrA, bcscB, csrMask, pointerOffsetA, pointerOffsetB, blockSizeY);
            csrBlockCold = csrBlockCnew;
            csrMask      = updateMask(csrMask, csrBlockCnew);
            continue;
        }
        csrBlockCnew = bmmPerBlock(bcsrA, bcscB, csrMask, pointerOffsetA, pointerOffsetB, blockSizeY);

        // no need update mask last iteration
        if (i != indicesOfCommon.size() - 2) {
            csrMask = updateMask(csrMask, csrBlockCnew);
        }

        csrBlockCnew = updateBlockC(csrBlockCnew, csrBlockCold);
        csrBlockCold = csrBlockCnew;
    }

    return csrBlockCnew;
}

/**
 *  Update the mask based on the previous result and the new one
 */
CSX updateMask(const CSX& csrNewMask, const CSX& csrOldMask)
{
    uint32_t rows = csrNewMask.pointer.size() - 1;
    uint32_t nnz  = 0;

    CSX ret;
    ret.pointer.resize(rows + 1);
    ret.pointer[0] = 0;

    uint32_t prevSize = 0;
    for (uint32_t i = 0; i < rows; i++) {
        prevSize = ret.indices.size();
        std::set_symmetric_difference(csrNewMask.indices.begin() + csrNewMask.pointer[i],
                                      csrNewMask.indices.begin() + csrNewMask.pointer[i + 1],
                                      csrOldMask.indices.begin() + csrOldMask.pointer[i],
                                      csrOldMask.indices.begin() + csrOldMask.pointer[i + 1],
                                      std::back_inserter(ret.indices));
        nnz += ret.indices.size() - prevSize;
        ret.pointer[i + 1] = nnz;
    }

    return ret;
}

/**
 *  Concatenate two CSR blocks
 *
 *  We expect the new CSR block to have different values of indices from the old CSR due to masking (check only the false (0) coordinates)
 *  Requirement: The values of the indices array should not have values in common, in order for the update to be correct, otherwise we will have duplicates
 */
CSX updateBlockC(const CSX& csrNew, const CSX& csrOld)
{
    uint32_t rows = csrNew.pointer.size() - 1;
    uint32_t nnz  = 0;

    CSX ret;
    ret.pointer.resize(rows + 1);
    ret.pointer[0] = 0;

    uint32_t prevSize = 0;
    for (uint32_t i = 0; i < rows; i++) {

        // empty lines
        // check if there is any performance difference
        // Update: there is a slight improvement including this check
        if (csrNew.pointer[i] == csrNew.pointer[i + 1] && csrOld.pointer[i] == csrOld.pointer[i + 1]) {
            ret.pointer[i + 1] = nnz;
            continue;
        }

        prevSize = ret.indices.size();
        std::merge(csrNew.indices.begin() + csrNew.pointer[i],
                   csrNew.indices.begin() + csrNew.pointer[i + 1],
                   csrOld.indices.begin() + csrOld.pointer[i],
                   csrOld.indices.begin() + csrOld.pointer[i + 1],
                   std::back_inserter(ret.indices));

        nnz += ret.indices.size() - prevSize;
        ret.pointer[i + 1] = nnz;
    }

    return ret;
}

/**
 * Fill a CSX format with ones
 */
CSX fillMaskOnes(uint32_t blockSizeY, uint32_t blockSizeX)
{
    CSX csrMask;
    csrMask.pointer.resize(blockSizeY + 1);
    csrMask.indices.resize(blockSizeY * blockSizeX);
    csrMask.pointer[0] = 0;

    uint32_t nnz = 0;
    for (uint32_t i = 0; i < blockSizeY; i++) {
        for (uint32_t j = 0; j < blockSizeX; j++) {
            csrMask.indices[nnz++] = j;
        }
        csrMask.pointer[i + 1] = nnz;
    }

    return csrMask;
}
