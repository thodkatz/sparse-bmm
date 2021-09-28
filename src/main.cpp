#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>
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
        csr.nnz = A.nnz;
        csr.n   = A.nRow + 1;

        mm2csc(argc, argv[1], csc, B);

        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed csx: " << tspan.count() << " (s) " << std::endl;
    }

    csc.nnz = B.nnz;
    csc.n   = B.nCol + 1;
    std::cout << "Rows: " << A.nRow << std::endl;
    std::cout << "Cols: " << A.nCol << std::endl;
    std::cout << "nnz:  " << A.nnz << std::endl;

    // printVector<uint32_t>(csr.pointer, " ");
    // printVector<uint32_t>(csr.indices, " ");
    // printVector<uint32_t>(csc.pointer, " ");
    // printVector<uint32_t>(csc.indices, " ");

    /* ---------- CSR Matrix-Matrix Multiplication ---------- */
    CSX csrMulResult;
    {
        steady_clock::time_point tic = steady_clock::now();
        // csrMulResult                 = csxMul(csr, csc);
        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed SpGEMM: " << tspan.count() << " (s) " << std::endl;
        std::cout << "Mul NNZ: " << csrMulResult.nnz << std::endl;
    }

    /* --- Write results for validation with python script -- */
    // csxWriteFile(csrMulResult, "test/csrMul.txt");

    CSX bcsr;
    CSX bcsc;
    A.blockSizeX = A.nRow / log10(A.nRow);
    A.blockSizeY = A.nCol / log10(A.nCol);

    {
        steady_clock::time_point tic = steady_clock::now();

        std::cout << "\nBlocks: (x,y) " << A.blockSizeX << "," << A.blockSizeY << std::endl;
        csr2bcsr(A, csr, bcsr);
        // printVector<uint32_t>(bcsr.pointer,",");
        // printVector<uint32_t>(bcsr.indices,",");

        csr2bcsr(A, csc, bcsc);
        // printVector<uint32_t>(bcsc.pointer,",");
        // printVector<uint32_t>(bcsc.indices,",");
        steady_clock::time_point toc = steady_clock::now();
        duration<double> tspan       = duration_cast<duration<double>>(toc - tic);
        std::cout << "Time elapsed bcsx: " << tspan.count() << " (s) " << std::endl;
    }

    return 0;
}