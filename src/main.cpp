#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include "sparsetools.hpp"
#include "utils.hpp"
#include "spgemm.hpp"

using namespace std::chrono;

int main(int argc, char* argv[])
{
    /* -------------------- Prepare data -------------------- */
    MatrixInfo A;
    CSX csr;

    MatrixInfo B;
    CSX csc;

    {
        steady_clock::time_point tic = steady_clock::now();

        mm2csr(argc, argv[1], csr, A);

        mm2csc(argc, argv[1], csc, B);

        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed csx: " << tspan.count() << " (s) " << std::endl;
    }

    std::cout << "Rows: " << A.nRow << std::endl;
    std::cout << "Cols: " << A.nCol << std::endl;
    std::cout << "nnz:  " << A.nnz << std::endl;

    // std::cout << "\nCSR" << std::endl;
    // printVector<uint32_t>(csr.pointer, " ");
    // printVector<uint32_t>(csr.indices, " ");
    // std::cout << "\nCSC" << std::endl;
    // printVector<uint32_t>(csc.pointer, " ");
    // printVector<uint32_t>(csc.indices, " ");

    /* ---------- CSR Matrix-Matrix Multiplication ---------- */
    CSX csrMulResult;
    {
        steady_clock::time_point tic = steady_clock::now();

        csrMulResult                 = csxMulSTL(csr, csc);

        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed SpGEMM: " << tspan.count() << " (s) " << std::endl;
        std::cout << "Mul NNZ: " << csrMulResult.indices.size() << std::endl;
    }

    /* --- Write results for validation with python script -- */
    csxWriteFile(csrMulResult, "test/csrMul.txt");

#ifdef BLOCK
    // std::cout << "\nBlocking Pad" << std::endl;
    BSXPad bcsr;
    BSXPad bcsc;
    A.blockSizeX = std::min((uint32_t)(A.nRow / log10(A.nRow)), A.nRow/2);
    A.blockSizeY = std::min((uint32_t)(A.nCol / log10(A.nCol)), A.nCol/2);
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;

    MatrixInfo copyA = A;
    MatrixInfo copyB = B;

    {
        steady_clock::time_point tic = steady_clock::now();

        std::cout << "\nBlocks size: (x,y) " << A.blockSizeX << "," << A.blockSizeY << std::endl;
        csr2bcsrPad(A, csr, bcsr);
        // printVector<uint32_t>(bcsr.pointer,",");
        // printVector<uint32_t>(bcsr.indices,",");

        csc2bcscPad(B, csc, bcsc);
        //printVector<uint32_t>(bcsc.pointer,",");
        //printVector<uint32_t>(bcsc.indices,",");
        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed bcsx: " << tspan.count() << " (s) " << std::endl;
    }
 
    CSX revertBlocked;
    bcsr2csrPad(A,bcsr,revertBlocked);
    // std::cout << "\nRevert Blocking Pad" << std::endl;
    // printVector<uint32_t>(revertBlocked.pointer," ");
    // printVector<uint32_t>(revertBlocked.indices," ");

    CSX cscRevert;
    bcsc2cscPad(B,bcsc,cscRevert);
    // printVector<uint32_t>(cscRevert.pointer," ");
    // printVector<uint32_t>(cscRevert.indices," ");

    A = copyA;
    BSXNoPad bcsrNoPad;
    csr2bcsrNoPad(A,csr,bcsrNoPad);
    // std::cout << "\nBlocking No Pad" << std::endl;
    // printVector<uint32_t>(bcsrNoPad.indices," ");
    // printVector<uint32_t>(bcsrNoPad.pointer," ");
    // printVector<uint32_t>(bcsrNoPad.idBlock," ");
    // printVector<uint32_t>(bcsrNoPad.blockPointer," ");

    B = copyB;
    BSXNoPad bcscNoPad;
    csc2bcscNoPad(B,csc,bcscNoPad);
    // printVector<uint32_t>(bcscNoPad.indices," ");
    // printVector<uint32_t>(bcscNoPad.pointer," ");
    // printVector<uint32_t>(bcscNoPad.idBlock," ");
    // printVector<uint32_t>(bcscNoPad.blockPointer," ");


    CSX revertBlockedNoPad;
    bcsr2csrNoPad(A,bcsrNoPad,revertBlockedNoPad);
    // std::cout << "\nRevert Blocking No Pad" << std::endl;
    // printVector<uint32_t>(revertBlockedNoPad.pointer," ");
    // printVector<uint32_t>(revertBlockedNoPad.indices," ");

    CSX cscRevertNoPad;
    bcsc2cscNoPad(B,bcscNoPad,cscRevertNoPad);
    // printVector<uint32_t>(cscRevertNoPad.pointer," ");
    // printVector<uint32_t>(cscRevertNoPad.indices," ");
#endif

    return 0;
}