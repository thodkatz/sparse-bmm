#include <iostream>
#include <vector>
#include <fstream>
#include "sparsetools.hpp"
#include "utils.hpp"

CSX csxMul(CSX& csr, CSX& csc)
{
    uint32_t n = csr.n;
    CSX mulResult;
    mulResult.n = n;

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

int main(int argc, char* argv[])
{
    /* -------------------- Prepare data -------------------- */
    MatrixInfo A;
    CSX csr;
    CSX csc;

    // std::vector<uint32_t> csrRow,csrCol,cscRow,cscCol;
    mm2csr(argc, argv[1], csr.pointer, csr.indices, A.nRow, A.nCol, A.nnz);
    csr.nnz = A.nnz;
    csr.n   = A.nRow;
    // mm2csc<uint32_t>(argc,argv,cscRow,cscCol,n_row,n_col,nnz);

    MatrixInfo B;
    mm2csc(argc, argv[1], csc.indices, csc.pointer, B.nRow, B.nCol, B.nnz);
    csc.nnz = B.nnz;
    csc.n   = B.nCol;

    std::cout << "Matrix A: " << A.nRow << "," << A.nCol << "," << A.nnz
              << std::endl;
    std::cout << "Matrix B: " << B.nRow << "," << B.nCol << "," << B.nnz
              << std::endl;
    printVector<uint32_t>(csr.pointer, " ");
    printVector<uint32_t>(csr.indices, " ");
    printVector<uint32_t>(csc.pointer, " ");
    printVector<uint32_t>(csc.indices, " ");

    CSX csrMulResult = csxMul(csr, csc);
    std::cout << "\nResult" << std::endl;
    std::cout << csrMulResult.nnz << " " << csrMulResult.n << std::endl;
    //printVector<uint32_t>(mulResult.pointer, " ");
    //printVector<uint32_t>(mulResult.indices, " ");

    // write to file the results and run python script to compate if matrix mul done right
    std::ofstream csrMulFile;
    csrMulFile.open("test/csrMul.txt");
    for (uint32_t i = 0; i <= csrMulResult.n; i++)
    {
        csrMulFile << csrMulResult.pointer[i];
        if(i!=csrMulResult.n) {
            csrMulFile << ",";
        }
    }
    csrMulFile << "\n";
    for (uint32_t i = 0; i < csrMulResult.nnz; i++)
    {
        csrMulFile << csrMulResult.indices[i];
        if(i!=csrMulResult.nnz-1) {
            csrMulFile << ",";
        }
    }
    csrMulFile.close();

    return 0;
}