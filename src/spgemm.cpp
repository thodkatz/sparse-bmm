#include <array>
#include <algorithm>
#include <utility>
#include "spgemm.hpp"
#include "sparsetools.hpp"
#include "utils.hpp"

/**
 *  C = A * B, where A CSR, B CSC and C CSR format
 */
CSX csxMul(const CSX& csrA, const CSX& cscB)
{
    uint32_t rows = csrA.pointer.size() - 1;

    CSX mulResult;
    mulResult.pointer.resize(rows + 1);
    mulResult.pointer[0];

    uint32_t nnz = 0;
    for (uint32_t i = 0, idxCol = 0, idxRow = 0; i < rows; i++) {
        for (uint32_t j = 0; j < rows; j++) {
            idxCol = csrA.pointer[i];
            if (idxCol == csrA.pointer[i + 1]) {
                break;
            }
            idxRow = cscB.pointer[j];

            // find common elements row ith and column ith (sorted)
            // two pointers version
            while (idxCol < csrA.pointer[i + 1] && idxRow < cscB.pointer[j + 1]) {
                if (csrA.indices[idxCol] == cscB.indices[idxRow]) {
                    mulResult.indices.push_back(j);
                    nnz++;
                    break;
                }
                if (csrA.indices[idxCol] > cscB.indices[idxRow]) {
                    idxRow++;
                }
                if (csrA.indices[idxCol] < cscB.indices[idxRow]) {
                    idxCol++;
                }
            }
        }
        mulResult.pointer[i + 1] = nnz;
    }

    return mulResult;
}

/**
 * An optimized compressed SpGEMM using STL variant of std::set_intersection()
 */
CSX csxMulSTL(const CSX& csrA, const CSX& cscB)
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

CSX csxMul(const CSX& csrA, const CSX& cscB, const CSX& csrF)
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

CSX bmmPerBlock(const BSXNoPad& csrA, const BSXNoPad& cscB, const CSX& csrF, uint32_t pointerOffsetA, uint32_t pointerOffsetB, uint32_t blockSize)
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
BSXNoPad bmmBlock(const MatrixInfo& F, const BSXNoPad& bcsrA, const BSXNoPad& bcscB, const BSXNoPad& bcsrF)
{
    // CSX csrResultBlock;
    BSXNoPad result;
    result.pointer.push_back(0);
    result.blockPointer.push_back(0);

    uint32_t indexBlockStartCol = 0;
    uint32_t indexBlockEndCol   = 0;
    uint32_t indexBlockStartRow = 0;
    uint32_t indexBlockEndRow   = 0;
    uint32_t nnzBlocks          = 0;
    for (uint32_t blockY = 0; blockY < F.numBlockY; blockY++) {
        for (uint32_t indexBlockX = bcsrF.blockPointer[blockY], idBlockCol = 0; indexBlockX < bcsrF.blockPointer[blockY + 1]; indexBlockX++) {
            idBlockCol = bcsrF.idBlock[indexBlockX];

            indexBlockStartCol = bcsrA.blockPointer[blockY];
            indexBlockEndCol   = bcsrA.blockPointer[blockY + 1];
            if (indexBlockStartCol == indexBlockEndCol) {
                break;
            }

            indexBlockStartRow = bcscB.blockPointer[idBlockCol];
            indexBlockEndRow   = bcscB.blockPointer[idBlockCol + 1];

            CSX csrBlockMask;
            getBlock(bcsrF, csrBlockMask, indexBlockX, F.blockSizeY);

            CSX csrResultBlock =
                subBlockMul(bcsrA, bcscB, csrBlockMask, F.blockSizeY, F.blockSizeX, indexBlockStartCol, indexBlockEndCol, indexBlockStartRow, indexBlockEndRow);

            // check if empty block
            if (csrResultBlock.pointer[F.blockSizeY] != 0) {
                result.idBlock.push_back(idBlockCol);
                appendResult(result, csrResultBlock, F.blockSizeY);
                nnzBlocks++;
            }
        }
        result.blockPointer.push_back(nnzBlocks);
    }

    return result;
}

void getBlock(const BSXNoPad& bcsx, CSX& block, uint32_t nnzBlocksPassed, uint32_t blockSizeY)
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
 *  Outer block level
 */
void appendResult(BSXNoPad& result, const CSX& csrResultBlock, uint32_t blockSizeY)
{
    uint32_t offset = result.pointer[result.pointer.size() - 1];

    // append pointer and indices content
    result.pointer.insert(result.pointer.end(), csrResultBlock.pointer.begin() + 1, csrResultBlock.pointer.end());
    result.indices.insert(result.indices.end(), csrResultBlock.indices.begin(), csrResultBlock.indices.end());

    std::for_each(result.pointer.end() - blockSizeY, result.pointer.end(), [offset](uint32_t& i) { i += offset; });
}

/**
 *  Outer block level
 */
CSX subBlockMul(const BSXNoPad& bcsrA,
                const BSXNoPad& bcscB,
                CSX& csrMask,
                uint32_t blockSizeY,
                uint32_t blockSizeX,
                uint32_t indexBlockStartCol,
                uint32_t indexBlockEndCol,
                uint32_t indexBlockStartRow,
                uint32_t indexBlockEndRow)
{
    // find the indices of the common elements
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

        if (i != indicesOfCommon.size() - 2) { // change this to 2
            csrMask = updateMask(csrMask, csrBlockCnew); // no need to calculate last iteration
        }

        csrBlockCnew = updateBlockC(csrBlockCnew, csrBlockCold);
        csrBlockCold = csrBlockCnew;
    }

    return csrBlockCnew;
}

/**
 *  Inner block level
 *
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
 *  Inner block level
 *
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
        // if (csrNew.pointer[i] == csrNew.pointer[i + 1] && csrOld.pointer[i] == csrOld.pointer[i + 1]) {
        //     ret.pointer[i+1] = nnz;
        //     continue;
        // }

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