#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <chrono>
#include "sparsetools.hpp"

using namespace std::chrono;

class Timer {
    private:
    std::string text;
    steady_clock::time_point tic;

    public:
    Timer(std::string text) {
        this->text = text;
        tic = steady_clock::now();
    }

    ~Timer() {
        auto toc = steady_clock::now();
        duration<double> tspan = duration_cast<duration<double>>(toc-tic);
        std::cout << "\n" << text << tspan.count() << " (s) \n" << std::endl;
    }
};

/**
 * Two pointer technique to find the first common element
 */
template <class InputIterator1, class InputIterator2>
bool hasCommon(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2)
{
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2)
            ++first1;
        else if (*first2 < *first1)
            ++first2;
        else {
            return true;
        }
    }
    return false;
}


/* ------------------- Print utilities ------------------ */
void printVector(const std::vector<uint32_t>& arr, std::string formatter);

void printBSX(const BSXNoPad& bcsx);

void printCSX(const CSX& csx);

void printCoo(std::vector<uint32_t>& rows, std::vector<uint32_t>& cols);

enum class sparseType { CSR, CSC };
void toDense(CSX csx, uint32_t rows, uint32_t cols, sparseType type, uint32_t pointerOffset, uint32_t indicesOffset);

/**
 *  Write \p csr data structure to file with name \p filename
 */
void csxWriteFile(CSX& csr, std::string filename);

void isEqualCSX(const CSX& reference, const CSX& got);

#endif