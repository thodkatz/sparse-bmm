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

#ifdef HYBRID
#define OPENMPI
#define OPENMP
#include <mpi.h>
// define tags to group types of communication with respect to MASTER (rank 0) and the rest WORKERS
#define FROM_MASTER  1
#define FROM_WORKERS 2
#endif

void print(std::vector<int>& arr) {
    for(auto i : arr) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

void readInput(char* argv[], MatrixInfo& A, MatrixInfo& B, MatrixInfo& F, CSX& csrA, CSX& cscB, CSX& csrF)
{
    {
        Timer time("Time converting input to CSX: \n");
        mm2csr(argv[1], csrA, A);
        mm2csc(argv[2], cscB, B);
        mm2csr(argv[3], csrF, F);
    }

    std::cout << "A: Rows: " << A.nRow << ", "
              << ", Cols: " << A.nCol << ", "
              << ", nnz: " << A.nnz << std::endl;
    std::cout << "B: Rows: " << B.nRow << ", "
              << ", Cols: " << B.nCol << ", "
              << ", nnz: " << B.nnz << std::endl;
    std::cout << "F: Rows: " << F.nRow << ", "
              << ", Cols: " << F.nCol << ", "
              << ", nnz: " << F.nnz << std::endl;

    if (A.nCol != B.nRow) {
        std::cout << "Multiplication dimensions mismatch" << std::endl;
        exit(-1);
    }
    if (F.nRow != A.nRow || F.nCol != B.nCol) {
        std::cout << "Masked matrix dimensions mismatch" << std::endl;
        exit(-1);
    }
}

int main(int argc, char* argv[])
{
#ifdef OPENMPI
    int numtasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    std::cout << "Number of processors: " << numtasks << std::endl;
    uint32_t numWorkers = numtasks - 1;
    MPI_Status status;

    if (numtasks < 2) {
        std::cout << "Run the application with at least 2 processes\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }
#endif

    if (argc < 4) {
        printf("Missed command line arguements\n");
        fprintf(stderr, "Usage: ./bin [martix-market-filename A] [matrix-market-filename B] [matrix-market-filename F]\n");
        exit(-1);
    }

#ifdef OPENMPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MatrixInfo A, B, F, C;
    CSX csrA, cscB, csrF;
    uint32_t offsetRows = 0;

    std::vector<int> testPointer;
    std::vector<int> testIndices;

    if (rank == 0) {
        readInput(argv, A, B, F, csrA, cscB, csrF);
        uint32_t averageRows = A.nRow / numWorkers;
        uint32_t extra       = A.nRow % numWorkers;

        C.nRow = A.nRow;
        C.nCol = B.nCol;

        // select and send a range of rows of matrices A and F to workers
        uint32_t rows       = 0;
        uint32_t nnzA       = 0;
        uint32_t nnzF       = 0;
        uint32_t nnzOffsetA = 0;
        uint32_t nnzOffsetF = 0;
        for (uint32_t dest = 1; dest <= numWorkers; dest++) {
            rows = (dest <= extra) ? averageRows + 1 : averageRows;
            nnzA = csrA.pointer[rows] - csrA.pointer[offsetRows];
            nnzF = csrF.pointer[rows] - csrF.pointer[offsetRows];

            MPI_Send(&offsetRows, 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&nnzA, 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&nnzF, 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);

            MPI_Send(&csrA.pointer + offsetRows, rows + 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&csrA.indices + nnzOffsetA, nnzA, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);

            MPI_Send(&csrF.pointer + offsetRows, rows + 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&csrF.indices + nnzOffsetF, nnzF, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);

            offsetRows += rows;
            nnzOffsetA += nnzA;
            nnzOffsetF += nnzF;

            for(int i = 0; i < cscB.pointer.size(); i++) {
                testPointer.push_back(cscB.pointer[i]);
            }
            for(int i = 0; i < cscB.indices.size(); i++) {
                testIndices.push_back(cscB.indices[i]);
            }
        }
    }

    // broadcast matrix B to all workers
    MPI_Bcast(&B.nRow, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B.nCol, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B.nnz, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    if(rank > 0) {
        cscB.pointer.resize(B.nCol+1);
        cscB.indices.resize(B.nnz);
        //testPointer.resize(B.nCol+1);
        //testIndices.resize(B.nnz);
    }
    MPI_Bcast(&cscB.pointer.front(), B.nCol + 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cscB.indices.front(), B.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    CSX csrRet;
    if (rank > 0) {
        MPI_Recv(&offsetRows, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&A.nRow, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&A.nnz, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&F.nnz, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        F.nRow = A.nRow;
        F.nCol = B.nCol;
        A.nCol = B.nRow;

        MPI_Recv(&csrA.pointer, A.nRow + 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&csrA.indices, A.nnz, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);

        MPI_Recv(&csrF.pointer, F.nRow, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&csrF.indices, F.nnz, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);

        /* ---------------------- Blocking ---------------------- */
        BSXNoPad bcsrA;
        BSXNoPad bcscB;
        BSXNoPad bcsrF;

        // A.blockSizeX = (uint32_t)(A.nRow / 20);
        // A.blockSizeY = (uint32_t)(A.nCol / 20);
        A.blockSizeX = 2;
        A.blockSizeY = 2;
        B.blockSizeX = A.blockSizeX;
        B.blockSizeY = A.blockSizeY;
        F.blockSizeX = B.blockSizeY;
        F.blockSizeY = A.blockSizeY;

        csr2bcsrNoPad(A, csrA, bcsrA);
        csc2bcscNoPad(B, cscB, bcscB);
        csr2bcsrNoPad(F, csrF, bcsrF);

        BSXNoPad ret;
        {
            Timer time("Time BLOCK-BMM\n");
            ret = bmmBlock(A, bcsrA, bcscB, bcsrF);
        }

        C     = F;
        C.nnz = ret.indices.size();

        {
            Timer time("Time BLOCK-BMM result convert BSX->CSX\n");
            bcsr2csrNoPad(C, ret, csrRet);
        }

        MPI_Reduce(&C.nnz, &C.nnz, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Gather results
    CSX result;
    if (rank == 0) {
        result.pointer.resize(A.nRow + 1);
        result.indices.resize(C.nnz);
    }

    int nnz = C.nnz;

    MPI_Gather(&csrRet.pointer, nnz, MPI_UINT32_T, &csrRet.pointer, nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Gather(&csrRet.indices, nnz, MPI_UINT32_T, &csrRet.indices, nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        csxWriteFile(result, "test/csrMul.txt");
    }

    MPI_Finalize();
#endif

#ifdef SERIAL
    /* -------------------- Prepare data -------------------- */
    MatrixInfo A, B, F;
    CSX csrA, cscB, csrF;
    readInput(argv, A, B, F, csrA, cscB, csrF);

#ifdef DEBUG
    std::cout << "\nA";
    printCSX(csrA);

    std::cout << "\nB";
    printCSX(cscB);

    std::cout << "\nF";
    printCSX(csrF);
#endif

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
    // A.blockSizeX = std::min((uint32_t)(A.nRow / log2(A.nRow)), A.nRow / 2);
    // A.blockSizeY = std::min((uint32_t)(A.nCol / log2(A.nCol)), A.nCol / 2);
    // A.blockSizeX = std::min((uint32_t)(A.nRow / sqrt(A.nRow)), A.nRow / 2);
    // A.blockSizeY = std::min((uint32_t)(A.nCol / sqrt(A.nCol)), A.nCol / 2);
    // A.blockSizeX = (uint32_t)(A.nRow / pow(A.nRow,1.0/3.0));
    // A.blockSizeY = (uint32_t)(A.nCol / pow(A.nCol,1.0/3.0));
    A.blockSizeX = std::min((uint32_t)(A.nRow / 20), A.nRow / 2);
    A.blockSizeY = std::min((uint32_t)(A.nCol / 20), A.nCol / 2);
    // A.blockSizeX = 2;
    // A.blockSizeY = 1;
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;
    F.blockSizeX = B.blockSizeY;
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

#ifdef DEBUG
    printBSX(bcsrA);
    printBSX(bcscB);
    printBSX(bcsrF);
    std::cout << "\nA dense:\n";
    toDense(csrA, A.nRow, A.nCol, sparseType::CSR, 0, 0);

    std::cout << "\nB dense:\n";
    toDense(cscB, B.nRow, B.nCol, sparseType::CSC, 0, 0);

    std::cout << "\nF dense:\n";
    toDense(csrF, F.nRow, F.nCol, sparseType::CSR, 0, 0);
#endif

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
        ret = bmmBlock(A, bcsrA, bcscB, bcsrF);
    }

    MatrixInfo C;
    C     = F;
    C.nnz = ret.indices.size();

    CSX csrRet;
    {
        Timer time("Time BLOCK-BMM result convert BSX->CSX\n");
        bcsr2csrNoPad(C, ret, csrRet);
    }

#ifdef DEBUG

    std::cout << "\nBlocked BMM result";
    toDense(csrRet, F.nRow, F.nCol, sparseType::CSR, 0, 0);
    printCSX(csrRet);

    std::cout << "\nCorrect result";
    toDense(csrCmask, F.nRow, F.nCol, sparseType::CSR, 0, 0);
    printCSX(csrCmask);
#endif

    std::cout << "\nBlock BMM validation\n";
    isEqualCSX(csrCmask, csrRet);

#endif // SERIAL
    return 0;
}