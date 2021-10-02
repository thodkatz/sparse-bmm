#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>
#include <random>
#include <algorithm>
#include "sparsetools.hpp"
#include "utils.hpp"
#include "spgemm.hpp"

int main(int argc, char* argv[])
{
    if (argc < 4) {
        printf("Missed command line arguements\n");
        fprintf(stderr, "Usage: ./bin [martix-market-filename A] [matrix-market-filename B] [matrix-market-filename F]\n");
        exit(-1);
    }

    /* -------------------- Prepare data -------------------- */
    MatrixInfo A;
    CSX csrA;

    MatrixInfo B;
    CSX cscB;

    MatrixInfo F;
    CSX csrF;

    {
        Timer time("Time converting input to CSX: \n");
        mm2csr(argv[1], csrA, A);
        mm2csc(argv[2], cscB, B);
        mm2csr(argv[3], csrF, F);
    }

    std::cout << "A: Rows: " << A.nRow << ", " << A.nCol << ", " << A.nnz << std::endl;
    std::cout << "B: Rows: " << B.nRow << ", " << B.nCol << ", " << B.nnz << std::endl;
    std::cout << "F: Rows: " << F.nRow << ", " << F.nCol << ", " << F.nnz << std::endl;

    if (A.nCol != B.nRow) {
        std::cout << "Multiplication dimensions mismatch" << std::endl;
        exit(-1);
    }
    if (F.nRow != A.nRow || F.nCol != B.nCol) {
        std::cout << "Masked matrix dimensions mismatch" << std::endl;
        exit(-1);
    }

    /* ------- Masked CSR Matrix-Matrix Multiplication ------ */
    CSX csrCmask;
    {
        Timer time("Time masked serial no blocking SpGEMM: \n");
        csrCmask = csxMul(csrA, cscB, csrF);
    }
    std::cout << "C.nnz: " << csrCmask.indices.size() << std::endl;

    /* --- Write results for validation with python script -- */
    csxWriteFile(csrCmask, "test/csrMul.txt");

    // ====================================================== //
    // ====================== BLOCKING ====================== //
    // ====================================================== //

    // A.blockSizeX = std::min((uint32_t)(A.nRow / log10(A.nRow)), A.nRow / 2);
    // A.blockSizeY = std::min((uint32_t)(A.nCol / log10(A.nCol)), A.nCol / 2);
    A.blockSizeX = std::min((uint32_t)(A.nRow / log2(A.nRow)), A.nRow / 2);
    A.blockSizeY = std::min((uint32_t)(A.nCol / log2(A.nCol)), A.nCol / 2);
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;
    F.blockSizeX = A.blockSizeX;
    F.blockSizeY = A.blockSizeY;

    BSXNoPad bcsrA;
    BSXNoPad bcscB;
    BSXNoPad bcsrF;
    {
        Timer time("Time creating BSX\n");
        csr2bcsrNoPad(A, csrA, bcsrA);
        csc2bcscNoPad(B, cscB, bcscB);
        csr2bcsrNoPad(F, csrF, bcsrF);
    }

    std::cout << "A: "
              << "blockX: " << A.blockSizeX << ", numBlocksX: " << A.numBlockX << ", blockY: " << A.blockSizeY << ", numBlocksY: " << A.numBlockY << std::endl;
    std::cout << "B: "
              << "blockX: " << B.blockSizeX << ", numBlocksX: " << B.numBlockX << ", blockY: " << B.blockSizeY << ", numBlocksY: " << B.numBlockY << std::endl;
    std::cout << "F: "
              << "blockX: " << F.blockSizeX << ", numBlocksX: " << F.numBlockX << ", blockY: " << F.blockSizeY << ", numBlocksY: " << F.numBlockY << std::endl;

    CSX csrArevert;
    CSX cscBrevert;
    CSX csrFrevert;
    {
        Timer time("Time Revert Blocking\n");
        bcsr2csrNoPad(A, bcsrA, csrArevert);
        bcsc2cscNoPad(B, bcscB, cscBrevert);
        bcsr2csrNoPad(F, bcsrF, csrFrevert);
    }

    std::cout << "\nValidation" << std::endl;

    /* ---------- Validation Blocked data structure --------- */
    std::cout << "Check conversion from BSX -> CSX\n";
    isEqualCSX(csrA, csrArevert);
    isEqualCSX(cscB, cscBrevert);
    isEqualCSX(csrF, csrFrevert);

    std::cout << "\nBLOCK BMM" << std::endl;

    BSXNoPad ret;
    {
        Timer time("Time BLOCK-BMM\n");
        ret = bmmBlock(F, bcsrA, bcscB, bcsrF);
    }

    MatrixInfo C;
    C     = F;
    C.nnz = ret.indices.size();

    CSX csrRet;
    {
        Timer time("Time BLOCK-BMM result convert BSX->CSX\n");
        bcsr2csrNoPad(C, ret, csrRet);
    }

    std::cout << "\nBlock BMM validation\n";
    isEqualCSX(csrCmask, csrRet);

    return 0;
}