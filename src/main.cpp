#include <iostream>
#include <vector>
#include <fstream>
#include "sparsetools.hpp"
#include "utils.hpp"
#include "spgemm.hpp"

int main(int argc, char* argv[])
{
    /* -------------------- Prepare data -------------------- */
    MatrixInfo A;
    CSX csr;

    mm2csr(argc, argv[1], csr, A);
    csr.nnz = A.nnz;
    csr.n   = A.nRow+1;

    MatrixInfo B;
    CSX csc;
    mm2csc(argc, argv[1], csc, B);
    csc.nnz = B.nnz;
    csc.n   = B.nCol+1;
    std::cout << "Rows: "<< A.nRow << std::endl;
    std::cout << "Cols: "<< A.nCol << std::endl;
    std::cout << "nnz:  "<< A.nnz << std::endl;

    // printVector<uint32_t>(csr.pointer, " ");
    // printVector<uint32_t>(csr.indices, " ");
    // printVector<uint32_t>(csc.pointer, " ");
    // printVector<uint32_t>(csc.indices, " ");

    /* ---------- CSR Matrix-Matrix Multiplication ---------- */
    CSX csrMulResult = csxMul(csr, csc);

    /* --- Write results for validation with python script -- */
    csxWriteFile(csrMulResult,"test/csrMul.txt");

    CSX bcsr;
    A.blockSizeX = 4;
    A.blockSizeY = 3;
    std::cout << "Blocks: (x,y) " << A.blockSizeX << "," << A.blockSizeY << std::endl;
    csr2bcsr(A,csr,bcsr);
    printVector<uint32_t>(bcsr.pointer,",");
    printVector<uint32_t>(bcsr.indices,",");

    CSX bcsc;
    csr2bcsr(A,csc,bcsc);
    //printVector<uint32_t>(bcsc.pointer,",");
    //printVector<uint32_t>(bcsc.indices,",");

    return 0;
}