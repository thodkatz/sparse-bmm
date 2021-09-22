template <class T>
void printVector(std::vector<T>& arr)
{
    std::cout << std::endl;
    for(const auto& i : arr) {
        std::cout << i << ",";
    }
    std::cout << std::endl;
}

template <class T>
void printCoo(std::vector<T>& rows, std::vector<T>& cols)
{
    if(rows.size() != cols.size()) return;

    std::cout << std::endl;
    for(T i = 0; i < rows.size(); i++) {
        std::cout<< rows[i] <<  " " << cols[i] << std::endl;
    }
    std::cout << std::endl;
}