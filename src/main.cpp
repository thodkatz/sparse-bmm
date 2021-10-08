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

void print(std::vector<uint32_t>& arr)
{
    for (auto i : arr) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

void printInt(std::vector<int>& arr)
{
    for (auto i : arr) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
#ifdef OPENMPI
    int numtasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    std::cout << "Number of processors: " << numtasks << std::endl;
    int numWorkers = numtasks - 1;
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
    CSX csrAmaster, csrFmaster;
    CSX csrA, cscB, csrF;

    std::vector<int> nnzSizeA;
    std::vector<int> nnzOffsetA;

    std::vector<int> nnzSizeF;
    std::vector<int> nnzOffsetF;

    std::vector<int> rowOffset;
    std::vector<int> rowSize;

    // read and distribute A and F
    if (rank == 0) {
        readInput(argv, A, B, F, csrAmaster, cscB, csrFmaster);
        int averageRows = A.nRow / numWorkers;
        int extra       = A.nRow % numWorkers;

        C.nRow = A.nRow;
        C.nCol = B.nCol;

        nnzSizeA.resize(numtasks);
        nnzOffsetA.resize(numtasks);

        nnzSizeF.resize(numtasks);
        nnzOffsetF.resize(numtasks);

        rowOffset.resize(numtasks);
        rowSize.resize(numtasks);

        for (int dest = 1, rows = 0, countRows = 0; dest <= numWorkers; dest++) {
            rowSize[dest]    = (dest <= extra) ? averageRows + 1 : averageRows;
            rowOffset[dest]  = rowOffset[dest - 1] + rowSize[dest - 1];
            rows             = rowSize[dest];
            countRows        = rowOffset[dest];
            nnzSizeA[dest]   = csrAmaster.pointer[countRows + rows] - csrAmaster.pointer[countRows];
            nnzSizeF[dest]   = csrFmaster.pointer[countRows + rows] - csrFmaster.pointer[countRows];
            nnzOffsetA[dest] = nnzOffsetA[dest - 1] + nnzSizeA[dest - 1];
            nnzOffsetF[dest] = nnzOffsetF[dest - 1] + nnzSizeF[dest - 1];

            // send size to each process to allocate memory
            MPI_Send(&nnzSizeA[dest], 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&nnzSizeF[dest], 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&rowSize[dest], 1, MPI_UINT32_T, dest, FROM_MASTER, MPI_COMM_WORLD);

            // remove pointer offsets
            std::for_each(
                csrAmaster.pointer.begin() + countRows, csrAmaster.pointer.begin() + countRows + rows, [start = csrAmaster.pointer[countRows]](uint32_t& i) { i -= start; });
            std::for_each(
                csrFmaster.pointer.begin() + countRows, csrFmaster.pointer.begin() + countRows + rows, [start = csrFmaster.pointer[countRows]](uint32_t& i) { i -= start; });
        }
    }

    // broadcast matrix B to all workers
    // just read it actually, but keep it for now
    {
        if (rank == 0)
            Timer time("Broadcasting B\n");
        MPI_Bcast(&B.nRow, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(&B.nCol, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(&B.nnz, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        if (rank > 0) {
            cscB.pointer.resize(B.nCol + 1);
            cscB.indices.resize(B.nnz);
        }
        MPI_Bcast(&cscB.pointer.front(), B.nCol + 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(&cscB.indices.front(), B.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }

    // receive matrix info and allocate memory
    if (rank > 0) {
        MPI_Recv(&A.nnz, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&F.nnz, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        MPI_Recv(&A.nRow, 1, MPI_UINT32_T, 0, FROM_MASTER, MPI_COMM_WORLD, &status);
        F.nRow = A.nRow;
        A.nCol = B.nRow;
        F.nCol = B.nCol;

        csrA.pointer.resize(A.nRow);
        csrA.indices.resize(A.nnz);
        csrF.pointer.resize(F.nRow);
        csrF.indices.resize(F.nnz);
    }

    {
        if (rank == 0)
            Timer time("Master scattering data\n");
        MPI_Scatterv(&csrAmaster.pointer.front(), &rowSize.front(), &rowOffset.front(), MPI_UINT32_T, &csrA.pointer.front(), A.nRow, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrAmaster.indices.front(), &nnzSizeA.front(), &nnzOffsetA.front(), MPI_UINT32_T, &csrA.indices.front(), A.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrFmaster.pointer.front(), &rowSize.front(), &rowOffset.front(), MPI_UINT32_T, &csrF.pointer.front(), F.nRow, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrFmaster.indices.front(), &nnzSizeF.front(), &nnzOffsetF.front(), MPI_UINT32_T, &csrF.indices.front(), F.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }

    // add the count of nnz to the end of the pointer array
    if (rank > 0) {
        csrA.pointer.push_back(A.nnz);
        csrF.pointer.push_back(F.nnz);
    }

    CSX csrRet;
    std::vector<int> nnzRetSize;
    std::vector<int> nnzRetOffset;
    if (rank == 0) {
        nnzRetSize.resize(numtasks);
        nnzRetOffset.resize(numtasks);
    }

    if (rank > 0) {
        /* ---------------------- Blocking ---------------------- */
        BSXNoPad bcsrA;
        BSXNoPad bcscB;
        BSXNoPad bcsrF;

        A.blockSizeY = (uint32_t)(A.nRow / 20); // the slices are mapped to the potential threads that can be used for multithreading
        A.blockSizeX = (uint32_t)(A.nCol / 1); // small block size of course reduce the block overhead and have better results
        B.blockSizeX = A.blockSizeX;
        B.blockSizeY = A.blockSizeY;
        F.blockSizeX = B.blockSizeY;
        F.blockSizeY = A.blockSizeY;

        csr2bcsrNoPad(A, csrA, bcsrA);
        csc2bcscNoPad(B, cscB, bcscB);
        csr2bcsrNoPad(F, csrF, bcsrF);

        std::cout << "A: "
                  << "blockX: " << A.blockSizeX << ", numBlocksX: " << A.numBlockX << ", blockY: " << A.blockSizeY << ", numBlocksY: " << A.numBlockY << std::endl;
        std::cout << "B: "
                  << "blockX: " << B.blockSizeX << ", numBlocksX: " << B.numBlockX << ", blockY: " << B.blockSizeY << ", numBlocksY: " << B.numBlockY << std::endl;
        std::cout << "F: "
                  << "blockX: " << F.blockSizeX << ", numBlocksX: " << F.numBlockX << ", blockY: " << F.blockSizeY << ", numBlocksY: " << F.numBlockY << std::endl;

        std::cout << "A: Rows: " << A.nRow << ", "
                  << "Cols: " << A.nCol << ", "
                  << "nnz: " << A.nnz << std::endl;
        std::cout << "B: Rows: " << B.nRow << ", "
                  << "Cols: " << B.nCol << ", "
                  << "nnz: " << B.nnz << std::endl;
        std::cout << "F: Rows: " << F.nRow << ", "
                  << "Cols: " << F.nCol << ", "
                  << "nnz: " << F.nnz << std::endl;

        std::cout << "\nWorker:" << rank << " BLOCK BMM" << std::endl;
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
    }

    // gather nnz from all processes
    MPI_Gather(&C.nnz, 1, MPI_UINT32_T, &nnzRetSize.front(), 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    // Gather results
    CSX result;
    if (rank == 0) {
        // accumulate nnzRet
        C.nnz = std::accumulate(nnzRetSize.begin(), nnzRetSize.end(), 0);
        for (int i = 1; i <= numWorkers; i++) {
            nnzRetOffset[i] += nnzRetOffset[i - 1] + nnzRetSize[i - 1];
        }
        result.pointer.resize(A.nRow);
        result.indices.resize(C.nnz);
    }

    MPI_Gatherv(&csrRet.pointer.front(), C.nRow, MPI_UINT32_T, &result.pointer.front(), &rowSize.front(), &rowOffset.front(), MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&csrRet.indices.front(), C.nnz, MPI_UINT32_T, &result.indices.front(), &nnzRetSize.front(), &nnzRetOffset.front(), MPI_UINT32_T, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // fix pointer offset
        for (int i = 1; i <= numWorkers; i++) {
            std::for_each(result.pointer.begin() + rowOffset[i], result.pointer.begin() + rowOffset[i] + rowSize[i], [start = nnzRetOffset[i]](uint32_t& i) { i += start; });
        }
        result.pointer.push_back(C.nnz);
        csxWriteFile(result, "test/csrMul.txt");
    }

    MPI_Finalize();
#endif // end MPI

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
    // A.blockSizeX = 3;
    // A.blockSizeY = 3;
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