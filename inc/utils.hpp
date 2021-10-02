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
 * source: https://stackoverflow.com/questions/27131628/check-whether-two-elements-have-a-common-element-in-c/27132690
 */
template <class ForwardIt1, class ForwardIt2>
bool has_common_elements(ForwardIt1 first, ForwardIt1 last, ForwardIt2 s_first, ForwardIt2 s_last)
{
    auto it   = first;
    auto s_it = s_first;
    while (it < last && s_it < s_last) {
        if (*it == *s_it) {
            return true;
        }

        *it < *s_it ? ++it : ++s_it; // increase the smaller of both
    }
    return false;
}

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
void printVector(const std::vector<uint32_t> arr, std::string formatter);
void printCSX(const CSX& csx);
void printCoo(std::vector<uint32_t>& rows, std::vector<uint32_t>& cols);

/**
 *  Write \p csr data structure to file with name \p filename
 */
void csxWriteFile(CSX& csr, std::string filename);

enum class sparseType { CSR, CSC };
void toDense(CSX csx, uint32_t rows, uint32_t cols, sparseType type, uint32_t pointerOffset, uint32_t indicesOffset);

#endif