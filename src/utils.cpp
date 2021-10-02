#include <fstream>
#include <iostream>
#include <algorithm>
#include "utils.hpp"
#include "sparsetools.hpp"

/**
 * Print CSX format to Dense Matrix
 */
void toDense(CSX csx, uint32_t rows, uint32_t cols, sparseType type, uint32_t pointerOffset, uint32_t indicesOffset)
{
    for (auto& i : csx.pointer) {
        i -= pointerOffset;
    }

    for (auto& i : csx.indices) {
        i -= indicesOffset;
    }

    if (type == sparseType::CSC) {
        std::swap(rows, cols);
    }

    std::vector<std::vector<uint32_t>> dense(rows, std::vector<uint32_t>(cols, 0));
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = csx.pointer[i]; j < csx.pointer[i + 1]; j++) {
            dense[i][csx.indices[j]] = 1;
        }
    }

    std::cout << std::endl;
    for (uint32_t i = 0; i < rows; i++) {
        for (uint32_t j = 0; j < cols; j++) {
            uint32_t value = (type == sparseType::CSR) ? dense[i][j] : dense[j][i];
            std::cout << value;
            std::string delimiter = (j != cols - 1) ? " " : "\n";
            std::cout << delimiter;
        }
    }
}

/**
 *  Check if two CSX objects are the same
 */
void isEqualCSX(const CSX& reference, const CSX& got)
{
    uint32_t n = reference.pointer.size() - 1;
    /* ------------------ Check dimensions ------------------ */
    if (reference.pointer.size() != got.pointer.size()) {
        std::cout << "Test Failed\n";
        std::cout << "Pointer dimension mismatch\n";
        exit(-1);
    }
    if (reference.indices.size() != got.indices.size()) {
        std::cout << "Test Failed\n";
        std::cout << "Indices dimension mismatch\n";
        exit(-1);
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

void csxWriteFile(CSX& csx, std::string filename)
{
    uint32_t n   = csx.pointer.size();
    uint32_t nnz = csx.indices.size();
    std::ofstream csxMulFile;
    csxMulFile.open(filename);
    for (uint32_t i = 0; i < n; i++) {
        csxMulFile << csx.pointer[i];
        if (i != n - 1) {
            csxMulFile << ",";
        }
    }
    csxMulFile << "\n";
    for (uint32_t i = 0; i < nnz; i++) {
        csxMulFile << csx.indices[i];
        if (i != nnz - 1) {
            csxMulFile << ",";
        }
    }
    csxMulFile.close();
}

void printVector(const std::vector<uint32_t> arr, std::string formatter)
{
    std::cout << std::endl;
    for (const auto& i : arr) {
        std::cout << i << formatter;
    }
    std::cout << std::endl;
};

void printCSX(const CSX& csx)
{
    std::cout << std::endl;
    for (const auto& i : csx.pointer) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for (const auto& i : csx.indices) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
};

void printBSX(const BSXNoPad& bcsx)
{
    std::cout << std::endl;
    for (const auto& i : bcsx.blockPointer) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for (const auto& i : bcsx.idBlock) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for (const auto& i : bcsx.pointer) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    for (const auto& i : bcsx.indices) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

void printCoo(std::vector<uint32_t>& rows, std::vector<uint32_t>& cols)
{
    if (rows.size() != cols.size())
        return;

    std::cout << std::endl;
    for (uint32_t i = 0; i < rows.size(); i++) {
        std::cout << rows[i] << " " << cols[i] << std::endl;
    }
    std::cout << std::endl;
};