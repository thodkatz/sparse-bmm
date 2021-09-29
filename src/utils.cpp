#include <fstream>
#include "sparsetools.hpp"

void csxWriteFile(CSX& csx, std::string filename)
{
    uint32_t n = csx.pointer.size()-1;
    uint32_t nnz = csx.indices.size();
    std::ofstream csxMulFile;
    csxMulFile.open(filename);
    for (uint32_t i = 0; i < n; i++) {
        csxMulFile << csx.pointer[i];
        if (i != n-1) {
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

/**
 * source: https://stackoverflow.com/questions/27131628/check-whether-two-elements-have-a-common-element-in-c/27132690
 */
template< class ForwardIt1, class ForwardIt2 >
bool has_common_elements( ForwardIt1 first, ForwardIt1 last, ForwardIt2 s_first, ForwardIt2 s_last )
{
    auto it=first;
    auto s_it=s_first;
    while(it<last && s_it<s_last)
    {
        if(*it==*s_it)
        {
            return true;
        }

        *it<*s_it ? ++it : ++s_it;  //increase the smaller of both
    }
    return false;
}