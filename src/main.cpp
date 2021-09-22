#include <iostream>
#include <vector>
#include "sparsetools.hpp"
#include "utils.hpp"

int main(int argc, char* argv[])
{
    /* -------------------- Prepare data -------------------- */
    uint32_t n   = 0;
    uint32_t nnz = 0;
    std::vector<uint32_t> rows,cols;
    mm2coo<uint32_t>(argc, argv, rows, cols, nnz, n);
    printCoo<uint32_t>(rows,cols);

    return 0;
}