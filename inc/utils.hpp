#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include "sparsetools.hpp"

template <class T>
void printVector(std::vector<T>& arr, std::string formatter)
{
    std::cout << std::endl;
    for (const auto& i : arr) {
        std::cout << i << formatter;
    }
    std::cout << std::endl;
};

template <class T>
void printCoo(std::vector<T>& rows, std::vector<T>& cols)
{
    if (rows.size() != cols.size())
        return;

    std::cout << std::endl;
    for (T i = 0; i < rows.size(); i++) {
        std::cout << rows[i] << " " << cols[i] << std::endl;
    }
    std::cout << std::endl;
};

void csxWriteFile(CSX& csr, std::string filename);

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

#endif