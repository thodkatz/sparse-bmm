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
#endif

#ifdef OPENMPI
#include <mpi.h>
#endif

#ifndef HYBRID
#ifdef OPENMP
#define SERIAL
#endif
#endif

int main(int argc, char* argv[])
{
#ifdef OPENMPI
    int numtasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
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
    if (rank == 0)
        std::cout << "Number of processors: " << numtasks << std::endl;
    MatrixInfo A, B, F, C;
    CSX csrA, cscB, csrF;
    CSX csrAmaster, csrFmaster;

    std::vector<int> nnzSizeA, nnzOffsetA;

    std::vector<int> nnzSizeF, nnzOffsetF;

    std::vector<int> rowOffset, rowSize;

    /* --------- Configure CSR-A and CSR-F slices for work distribution --------- */
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
            MPI_Send(&nnzSizeA[dest], 1, MPI_UINT32_T, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&nnzSizeF[dest], 1, MPI_UINT32_T, dest, 1, MPI_COMM_WORLD);
            MPI_Send(&rowSize[dest], 1, MPI_UINT32_T, dest, 1, MPI_COMM_WORLD);

            // remove pointer offsets
            std::for_each(
                csrAmaster.pointer.begin() + countRows, csrAmaster.pointer.begin() + countRows + rows, [start = csrAmaster.pointer[countRows]](uint32_t& i) { i -= start; });
            std::for_each(
                csrFmaster.pointer.begin() + countRows, csrFmaster.pointer.begin() + countRows + rows, [start = csrFmaster.pointer[countRows]](uint32_t& i) { i -= start; });
        }
    }

    /* ----------------- Broadcast matrix B ----------------- */
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

    /* ---------------- Receive matrix info for memory allocation --------------- */
    if (rank > 0) {
        MPI_Recv(&A.nnz, 1, MPI_UINT32_T, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&F.nnz, 1, MPI_UINT32_T, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&A.nRow, 1, MPI_UINT32_T, 0, 1, MPI_COMM_WORLD, &status);
        F.nRow = A.nRow;
        A.nCol = B.nRow;
        F.nCol = B.nCol;

        csrA.pointer.resize(A.nRow);
        csrA.indices.resize(A.nnz);
        csrF.pointer.resize(F.nRow);
        csrF.indices.resize(F.nnz);
    }

    /* ----------------------- Distribute CSR-A and CSR-F ----------------------- */
    {
        if (rank == 0)
            Timer time("Master scattering data\n");
        MPI_Scatterv(&csrAmaster.pointer.front(), &rowSize.front(), &rowOffset.front(), MPI_UINT32_T, &csrA.pointer.front(), A.nRow, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrAmaster.indices.front(), &nnzSizeA.front(), &nnzOffsetA.front(), MPI_UINT32_T, &csrA.indices.front(), A.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrFmaster.pointer.front(), &rowSize.front(), &rowOffset.front(), MPI_UINT32_T, &csrF.pointer.front(), F.nRow, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&csrFmaster.indices.front(), &nnzSizeF.front(), &nnzOffsetF.front(), MPI_UINT32_T, &csrF.indices.front(), F.nnz, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }

    // no overlapping is allowed when sending chuncks of data, so
    // add the count of nnz to the end of the pointer array
    if (rank > 0) {
        csrA.pointer.push_back(A.nnz);
        csrF.pointer.push_back(F.nnz);
    }

    CSX csrRet;
    // keep track of the bmm's nnz each process will return
    std::vector<int> nnzRetSize;
    std::vector<int> nnzRetOffset;
    if (rank == 0) {
        nnzRetSize.resize(numtasks);
        nnzRetOffset.resize(numtasks);
    }

    // start tracking time
    Timer time("Master BMM total, rank: " + std::to_string(rank) + "\n");

    /* ---------------- Per Process block-bmm --------------- */
    if (rank > 0) {
        BSX bcsrA, bcscB, bcsrF;

#ifdef HYBRID
        // for the hybrid mode we pick 49x1 blockSize, because the maximum pair of tasks will be 7 (processors) x 7 (threads per processor)
        uint32_t numBlockGoal    = 36;
        uint32_t numBlockPerProc = numBlockGoal / numWorkers;
        A.blockSizeY             = (uint32_t)(A.nRow / numBlockPerProc); // the slices are mapped to the potential threads that can be used for multithreading
        A.blockSizeX             = (uint32_t)(A.nCol / 1);               // small block size of course reduce the block overhead and have better results
#else
        uint32_t numBlockGoal    = 20;
        uint32_t numBlockPerProc = numBlockGoal / numWorkers;
        A.blockSizeY             = (uint32_t)(A.nRow / numBlockPerProc);
        A.blockSizeX             = (uint32_t)(A.nCol / 20);
#endif
        B.blockSizeX = A.blockSizeX;
        B.blockSizeY = A.blockSizeY;
        F.blockSizeX = B.blockSizeY;
        F.blockSizeY = A.blockSizeY;

        csr2bcsr(A, csrA, bcsrA);
        csc2bcsc(B, cscB, bcscB);
        csr2bcsr(F, csrF, bcsrF);

        std::cout << "Rank: " << rank << std::endl
                  << "A: blockX: " << A.blockSizeX << ", numBlocksX: " << A.numBlockX << ", blockY: " << A.blockSizeY << ", numBlocksY: " << A.numBlockY << std::endl
                  << "B: blockX: " << B.blockSizeX << ", numBlocksX: " << B.numBlockX << ", blockY: " << B.blockSizeY << ", numBlocksY: " << B.numBlockY << std::endl
                  << "F: blockX: " << F.blockSizeX << ", numBlocksX: " << F.numBlockX << ", blockY: " << F.blockSizeY << ", numBlocksY: " << F.numBlockY << std::endl
                  << "A: Rows: " << A.nRow << ", Cols: " << A.nCol << ", nnz: " << A.nnz << std::endl
                  << "B: Rows: " << B.nRow << ", Cols: " << B.nCol << ", nnz: " << B.nnz << std::endl
                  << "F: Rows: " << F.nRow << ", Cols: " << F.nCol << ", nnz: " << F.nnz << std::endl
                  << std::endl;

        BSX ret;
        {
            Timer time("Rank: " + std::to_string(rank) + " Time BLOCK-BMM\n");
            ret = bmmBlock(A, bcsrA, bcscB, bcsrF);
        }

        C     = F;
        C.nnz = ret.indices.size();

        {
            Timer time("Rank: " + std::to_string(rank) + " Time BLOCK-BMM result convert BSX->CSX\n");
            bcsr2csr(C, ret, csrRet);
        }
    }

    /* ------------------------- Gathering from workers and unifying result------------------------- */
    {
        Timer time("Gathering nnz, rank: " + std::to_string(rank) + "\n");
        MPI_Gather(&C.nnz, 1, MPI_UINT32_T, &nnzRetSize.front(), 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }

    // log time
    if (rank == 0) {
        std::ofstream logFile;
        logFile.open("log.csv", std::ios_base::app);
        std::string type;
#ifdef HYBRID
        type = "hybrid";
#else
        type                     = "openmpi";
#endif
        logFile << numWorkers << "," << type << "," << argv[1] << "," << time.elapsed() << "\n";
        logFile.close();
    }

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
    /* ------------------------------ Prepate Data ------------------------------ */
    MatrixInfo A, B, F;
    CSX csrA, cscB, csrF;
    readInput(argv, A, B, F, csrA, cscB, csrF);

    /* --------------------------- Masked CSR-CSC BMM --------------------------- */
    CSX csrCmask;
    {
        Timer time("Time masked serial no blocking SpGEMM: \n");
        csrCmask = bmm(csrA, cscB, csrF);
    }
    // std::cout << "C.nnz: " << csrCmask.indices.size() << std::endl;

    /* -------------------------------- Blocking -------------------------------- */
    // blockSizeX ~= blockSizeY is supported
    //A.blockSizeX = A.nRow / 20;
    //A.blockSizeY = A.nCol / 20;

    A.blockSizeX = A.nRow / 2;
    A.blockSizeY = A.nCol / 2;
    B.blockSizeX = A.blockSizeX;
    B.blockSizeY = A.blockSizeY;
    F.blockSizeX = B.blockSizeY;
    F.blockSizeY = A.blockSizeY;

    BSX bcsrA, bcscB, bcsrF;
    {
        Timer time("Time creating BSX\n");
        csr2bcsr(A, csrA, bcsrA);
        csc2bcsc(B, cscB, bcscB);
        csr2bcsr(F, csrF, bcsrF);
    }

    printBSX(bcsrA);
    toDense(csrA,4,4,sparseType::CSR,0,0);

    std::cout << "A: blockX: " << A.blockSizeX << ", numBlocksX: " << A.numBlockX << ", blockY: " << A.blockSizeY << ", numBlocksY: " << A.numBlockY << std::endl
              << "B: blockX: " << B.blockSizeX << ", numBlocksX: " << B.numBlockX << ", blockY: " << B.blockSizeY << ", numBlocksY: " << B.numBlockY << std::endl
              << "F: blockX: " << F.blockSizeX << ", numBlocksX: " << F.numBlockX << ", blockY: " << F.blockSizeY << ", numBlocksY: " << F.numBlockY << std::endl;

    std::cout << "\nBLOCK BMM" << std::endl;

    BSX ret;
    {
        Timer time("Time BLOCK-BMM\n");
        ret = bmmBlock(A, bcsrA, bcscB, bcsrF);

        // log time
        std::ofstream bmmFile;
        bmmFile.open("log.csv", std::ios_base::app);
        std::string type;
        uint32_t numOfThreads;
#ifdef OPENMP
        type         = "openmp";
        numOfThreads = std::stoi(std::getenv("OMP_NUM_THREADS"));
#else
        type                     = "serial";
        numOfThreads             = 1;
#endif
        bmmFile << numOfThreads << "," << type << "," << argv[1] << "," << time.elapsed() << "\n";
        bmmFile.close();
    }

    MatrixInfo C;
    C     = F;
    C.nnz = ret.indices.size();

    CSX csrRet;
    {
        Timer time("Time BLOCK-BMM result convert BSX->CSX\n");
        bcsr2csr(C, ret, csrRet);
    }

    // std::cout << "\nBlock BMM validation with non block BMM\n";
    // isEqualCSX(csrCmask, csrRet);

    /* ------------ Export results for validation with python script ------------ */
    csxWriteFile(csrRet, "test/csrMul.txt");

#endif // SERIAL
    return 0;
}