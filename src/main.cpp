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

using namespace std::chrono;

void isEqualCSX(const CSX& reference, const CSX& got)
{
    uint32_t n = reference.pointer.size() - 1;
    /* ------------------ Check dimensions ------------------ */
    if (reference.pointer.size() != got.pointer.size()) {
        std::cout << "Test Failed\n";
        std::cout << "Pointer dimension mismatch\n";
    }
    if (reference.indices.size() != got.indices.size()) {
        std::cout << "Test Failed\n";
        std::cout << "Indices dimension mismatch\n";
    }

    /* ------------------- Check pointers ------------------- */
    for (uint32_t i = 0; i <= n; i++) {
        if (reference.pointer[i] != got.pointer[i]) {
            std::cout << "Test Failed" << std::endl;
            std::cout << "Pointers mismatch" << std::endl;
            std::cout << "Expected";
            printVector(reference.pointer, " ");
            std::cout << "Got";
            printVector(got.pointer, " ");
            exit(-1);
        }
    }

    /* -------------------- Check indices ------------------- */
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = reference.pointer[i]; j < reference.pointer[i + 1]; j++) {
            if (reference.indices[j] != got.indices[j]) {
                std::cout << "Test Failed" << std::endl;
                std::cout << "Indices mismatch" << std::endl;
                std::cout << "Expected";
                printVector(reference.indices, " ");
                std::cout << "Got";
                printVector(got.indices, " ");

                exit(-1);
            }
        }
    }

    std::cout << "Test Passed!" << std::endl;
}

int main(int argc, char* argv[])
{
    // expecting a filename to read (./main <filename>)
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
        Timer time("Time converting input to compressed formats: \n");
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

    // std::cout << "\nCSR" << std::endl;
    // printVector<uint32_t>(csrA.pointer, " ");
    // printVector<uint32_t>(csrA.indices, " ");
    // std::cout << "\nCSC" << std::endl;
    // printVector<uint32_t>(cscB.pointer, " ");
    // printVector<uint32_t>(cscB.indices, " ");

    /* ---------- CSR Matrix-Matrix Multiplication ---------- */
    // CSX csrC;
    // {
    //     Timer time("Time serial no blocking SpGEMM: \n");
    //     csrC = csxMulSTL(csrA, cscB);
    // }
    // std::cout << "C.nnz: " << csrC.indices.size() << std::endl;

    /* ------- Masked CSR Matrix-Matrix Multiplication ------ */
    CSX csrCmask;
    {
        Timer time("Time masked serial no blocking SpGEMM: \n");
        csrCmask = csxMul(csrA, cscB, csrF);
    }
    std::cout << "C.nnz: " << csrCmask.indices.size() << std::endl;

;
    /* --- Masked Blocked CSR Matrix-Matrix Multiplication -- */

    /* --- Write results for validation with python script -- */
    csxWriteFile(csrCmask, "test/csrMul.txt");

#ifdef BLOCK
    // std::cout << "\nBlocking Pad" << std::endl;
    BSXPad bcsrA;
    BSXPad bcscB;
    A.blockSizeX = std::min((uint32_t)(A.nRow / log10(A.nRow)), A.nRow / 2);
    A.blockSizeY = std::min((uint32_t)(A.nCol / log10(A.nCol)), A.nCol / 2);
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;

    MatrixInfo copyA = A;
    MatrixInfo copyB = B;

    {
        steady_clock::time_point tic = steady_clock::now();

        std::cout << "\nBlocks size: (x,y) " << A.blockSizeX << "," << A.blockSizeY << std::endl;
        csrA2bcsrAPad(A, csrA, bcsrA);
        // printVector<uint32_t>(bcsrA.pointer,",");
        // printVector<uint32_t>(bcsrA.indices,",");

        cscB2bcscBPad(B, cscB, bcscB);
        // printVector<uint32_t>(bcscB.pointer,",");
        // printVector<uint32_t>(bcscB.indices,",");
        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed bcsx: " << tspan.count() << " (s) " << std::endl;
    }

    CSX revertBlocked;
    bcsr2csrAPad(A, bcsrA, revertBlocked);
    // std::cout << "\nRevert Blocking Pad" << std::endl;
    // printVector<uint32_t>(revertBlocked.pointer," ");
    // printVector<uint32_t>(revertBlocked.indices," ");

    CSX cscBRevert;
    bcsc2cscPad(B, bcscB, cscBRevert);
    // printVector<uint32_t>(cscBRevert.pointer," ");
    // printVector<uint32_t>(cscBRevert.indices," ");

#endif

    A.blockSizeX = std::min((uint32_t)(A.nRow / log10(A.nRow)), A.nRow / 2);
    A.blockSizeY = std::min((uint32_t)(A.nCol / log10(A.nCol)), A.nCol / 2);
    // A.blockSizeX = 3;
    // A.blockSizeY = 2;
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;
    F.blockSizeX = A.blockSizeX;
    F.blockSizeY = A.blockSizeY;

    BSXNoPad bcsrA;
    BSXNoPad bcscB;
    BSXNoPad bcsrF;
    {
        Timer time("Time Blocking\n");
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
        // printVector<uint32_t>(revertBlockedNoPad.pointer," ");
        // printVector<uint32_t>(revertBlockedNoPad.indices," ");

        bcsc2cscNoPad(B, bcscB, cscBrevert);
        // printVector<uint32_t>(cscBRevertNoPad.pointer," ");
        // printVector<uint32_t>(cscBRevertNoPad.indices," ");

        bcsr2csrNoPad(F, bcsrF, csrFrevert);
    }

    std::cout << "\nValidation" << std::endl;

    /* ---------- Validation Blocked data structure --------- */
    isEqualCSX(csrA, csrArevert);
    isEqualCSX(cscB, cscBrevert);
    isEqualCSX(csrF, csrFrevert);

    std::cout << "\nA dense:\n";
    toDense(csrA, A.nRow, A.nCol, sparseType::CSR, 0, 0);
    // printVector(bcsrA.blockPointer, " ");
    // printVector(bcsrA.idBlock, " ");
    // printVector(bcsrA.pointer, " ");
    // printVector(bcsrA.indices, " ");

    std::cout << "\nB dense:\n";
    toDense(cscB, B.nRow, B.nCol, sparseType::CSC, 0, 0);
    // printVector(bcscB.blockPointer, " ");
    // printVector(bcscB.idBlock, " ");
    // printVector(bcscB.pointer, " ");
    // printVector(bcscB.indices, " ");

    std::cout << "\nF dense:\n";
    toDense(csrF, F.nRow, F.nCol, sparseType::CSR, 0, 0);
    // printVector(bcsrF.blockPointer, " ");
    // printVector(bcsrF.idBlock, " ");
    // printVector(bcsrF.pointer, " ");
    // printVector(bcsrF.indices, " ");

#ifndef BMM_BLOCK
    BSXNoPad ret = bmmBlock(F, bcsrA, bcscB, bcsrF);

    MatrixInfo C;
    C = F;
    C.nnz = ret.indices.size();

    std::cout << "\nBlocked BMM result";
    CSX csrRet;
    bcsr2csrNoPad(C, ret, csrRet);
    toDense(csrRet, F.nRow, F.nCol, sparseType::CSR, 0, 0);
    printCSX(csrRet);

    std::cout << "\nCorrect result";
    toDense(csrCmask, F.nRow, F.nCol, sparseType::CSR, 0, 0);
    printCSX(csrCmask);

    isEqualCSX(csrCmask,csrRet);

#endif

#ifdef CHECK_INTERSECTION
    /* ------------- Check indices intersection ------------- */
    std::cout << std::endl;
    std::srand(unsigned(std::time(nullptr)));
    std::vector<int> vec1(10);
    std::generate(begin(vec1), end(vec1), std::rand);
    std::for_each(vec1.begin(), vec1.end(), [](int& a) { a %= 10; });
    std::for_each(vec1.begin(), vec1.end(), [](int& a) { std::cout << a << " "; });
    std::vector<int> vec2 = vec1;
    unsigned seed         = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(vec2.begin(), vec2.end(), std::default_random_engine(seed));
    vec2.resize(5);
    std::cout << std::endl;
    std::for_each(vec2.begin(), vec2.end(), [](int& a) { std::cout << a << " "; });
    std::cout << std::endl;

    std::vector<uint32_t> indicesOfCommon;
    indicesIntersection(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), std::back_inserter(indicesOfCommon));

    printVector(indicesOfCommon, " ");

#endif

#ifdef CHECK
    /* ----------------- Check updateMask() ----------------- */
    std::cout << "\nUpdate Mask, Symmetric Difference\n";

    CSX fillOnes = fillMaskOnes(A.nRow, B.nCol);
    // toDense(fillOnes,A.nRow,B.nCol,sparseType=CSR);

    CSX csrUpdateMask = updateMask(csrA, fillOnes);
    toDense(csrUpdateMask, A.nRow, B.nCol, sparseType = CSR);

    /* ------------ Check updateBlockC() function ----------- */
    std::cout << "\nUpdate Block, Merge\n";
    CSX csrUpdate = updateBlockC(csrA, csrUpdateMask);
    toDense(csrUpdate, A.nRow, A.nCol, sparseType = CSR);
#endif

    return 0;
}