#include "sparse_matrix.hpp"

namespace algebra {

    template<typename T, StorageOrdering Order>
    std::array<std::size_t, 2> SparseMatrix<T, Order>::convertCoords(std::size_t i, std::size_t j) const {
        return (Order == StorageOrdering::RowOrdering) ? std::array{i, j} : std::array{j, i};
    }

    // Constructor that takes the number of rows and columns
    template<typename T, StorageOrdering Order>
    SparseMatrix<T, Order>::SparseMatrix(std::size_t rows, std::size_t cols)
        : numrows(rows), numcols(cols), compressed(false) {}

    //another constructor that takes the values?
    

    // non const version: can add elements in uncompressed format, can only modify existing elements in compressed format
    template<typename T, StorageOrdering Order>
    T& SparseMatrix<T, Order>::operator()(std::size_t i, std::size_t j) {
        auto coords = convertCoords(i, j); // Get the coords according to the ordering
        if (!isCompressed()) {
            // Uncompressed format, I can add new elements 
            if(i >=numrows || j >=numcols){
                resize(i,j); //if I need to add an element that increases the size of the matrix
            }
            return uncompressed_data[coords]; // Creates new data if element not found 
        }
        else{
            // Compressed format, I CANNOT add new elements
            if(i >=numrows || j >=numcols){
                throw std::out_of_range("Matrix in compressed form, cannot add new elements!");
            }
            else{
            // Look for the element
            //if j is in the interval [ outer[inner[i]], outer[inner[i+1]] ), then A[i,j] exists
            auto i = coords[0];
            auto j = coords[1];
            std::size_t row_start = compressed_inner[i];
            std::size_t row_end = compressed_inner[i + 1];
            auto it = std::find (compressed_outer.begin() + row_start,
                                 compressedouter.begin() + row_end, j); 
        if (it != compressed_outer.begin() + row_end && *it == j) {
            std::size_t index = std::distance(compressed_outer.begin(), it);
            return compressed_data[index];
            //if it doesn't exist, raise an error
            else{
                throw std::out_of_range("Matrix in compressed form, cannot add new elements!");
                }
            }
        }       
    }
}


    // const version 
    template<typename T, StorageOrdering Order>
    const T& SparseMatrix<T, Order>::operator()(std::size_t i, std::size_t j) const {
        if (i >= numrows || j >= numcols) {
              throw std::out_of_range("Index out of boundary"); //if indexes out of range
        }
        auto coords = convertCoords(i, j);
        // Uncompressed matrix
        if !(compressed){
            auto it = uncompressed_data.find(coords);
            if (it != uncompressed_data.end()) {
                return it->second;
            }
            else{
                 throw std::out_of_range("Element not found");
            }
        }
        // Compressed format
        else{
            //if j is in the interval [ outer[inner[i]], outer[inner[i+1]] ), A[i,j] exists
            auto i = coords[0];
            auto j = coords[1];
            std::size_t row_start = compressed_inner[i];
            std::size_t row_end = compressed_inner[i + 1];
            auto it = std::find (compressed_outer.begin() + row_start,
                                 compressed_outer.begin() + row_end, j); 
        if (it != compressed_outer.begin() + row_end && *it == j) {
            std::size_t index = std::distance(compressed_outer.begin(), it);
            return compressed_data[index];
         }
        }
        static T default_value{};
        return default_value; //if element not found: default value??????????
    }


    // Resize method 
    template<typename T, StorageOrder Order>
        void SparseMatrix<T, Order>::resize(std::size_t rows, std::size_t cols) {
            if(!isCompressed){
                numrows = rows;
                numcols = cols;
            }
        }


    template<typename T, StorageOrdering Order>
    void SparseMatrix<T, Order>::compress() {
       if (compressed) {
            return; // Matrix is already compressed, no need to compress again
        }
        // Clear existing compressed data if any
        compressed_inner.clear();
        compressed_outer.clear();
        compressed_data.clear();
        if constexpr (Order == StorageOrdering::RowOrdering) {
            compressRowOrder();
        } else if constexpr (Order == StorageOrdering::ColumnOrdering) {
            compressColumnOrder();
        // Clear uncompressed data after compression
        uncompressed_data.clear();
        compressed = true; // Set compressed flag
    }

   template<typename T, StorageOrdering Order>
    void SparseMatrix<T, Order>::compressRowOrder() {
        // Reserve space for the indexes
        compressed_inner.reserve(numrows + 1); 
        compressed_outer.reserve(uncompressed_data.size()); 
        compressed_data.reserve(uncompressed_data.size()); 
        auto it= uncompressed_data.begin();
        const auto& coords = it->first;
        const auto& value = it->second;
        compressed_inner.push_back(coords[0]); 
        compressed_outer.push_back(coords[1]);
        compressed_data.pushback(value);

        std::size_t current_row = 0;
        while(it!= uncompressed_data.end() &&current_row < numrows){
            ++it;
            coords = it-> first;
            value = it->second;
            while (coords[0] > current_row) {
                // Move to the next row in compressed format
                current_row++;
                compressed_inner.push_back(coords[0]);
                compressed_outer.push_back(coords[1]);
                compressed_data.push_back(value); 
            }
            compressed_outer.push_back(coords[1]);
            compressed_data.push_back(value); 
        }
    }


    template<typename T, StorageOrdering Order>
    void SparseMatrix<T, Order>::compressColumnOrder() {
        // Reserve space for the indexes
        compressed_inner.reserve(numcols + 1); 
        compressed_outer.reserve(uncompressed_data.size()); 
        compressed_data.reserve(uncompressed_data.size()); 
        auto it= uncompressed_data.begin();
        const auto& coords = it->first;
        const auto& value = it->second;
        compressed_inner.push_back(coords[1]); 
        compressed_outer.push_back(coords[0]);
        compressed_data.pushback(value);

        std::size_t current_col = 0;
        while(it!= uncompressed_data.end() &&current_col < numcols){
            ++it;
            coords = it -> first;
            value = it -> second;
            while (coords[1] > current_col) {
                // Move to the next column in compressed format
                current_col++;
                compressed_inner.push_back(coords[1]);
                compressed_outer.push_back(coords[0]);
                compressed_data.push_back(value); 
            }
            // if the column is the same as previous
            compressed_outer.push_back(coords[0]);
            compressed_data.push_back(value); 
            }
    }


    template<typename T, StorageOrdering Order>
    bool SparseMatrix<T, Order>::isCompressed() const {
        return compressed;
    }

    template<typename T, StorageOrdering Order>
    void SparseMatrix<T, Order>::print() const {
        if (compressed) {
            // Print compressed format (CSR or CSC)
            std::cout << "Compressed Sparse Matrix:" << std::endl;
            // Print compressed_data, inner_index, outer_index as needed
        } else {
            std::cout << "Uncompressed Sparse Matrix:" << std::endl;
            for (const auto& [coords, value] : uncompressed_data) {
                std::cout << "(" << coords[0] << "," << coords[1] << "): " << value << std::endl;
            }
        }
    }

    
    
// Friend function for matrix-vector multiplication
template<typename T, StorageOrdering Order>
friend std::vector<T> operator*(const SparseMatrix<T, Order>& matrix, const std::vector<T>& vec) {
    std::vector<T> result(matrix.numrows, T{}); // Initialize result vector

    if (!matrix.isCompressed()) {
        // Uncompressed matrix multiplication
        for (const auto& [coords, value] : matrix.uncompressed_data) {
            std::size_t i = coords[0];
            std::size_t j = coords[1];
            result[i] += value * vec[j];
        }
    } else {
        // Compressed matrix multiplication (CSR or CSC)
        for (std::size_t i = 0; i < matrix.numrows; ++i) {
            std::size_t row_start = matrix.compressed_inner[i];
            std::size_t row_end = matrix.compressed_inner[i + 1];
            for (std::size_t k = row_start; k < row_end; ++k) {
                std::size_t j = matrix.compressed_outer[k];
                result[i] += matrix.compressed_data[k] * vec[j];
            }
        }
    }

    return result;
}
    












    // Explicit instantiation for the types we use
    template class SparseMatrix<int, StorageOrdering::RowOrdering>;
    template class SparseMatrix<int, StorageOrdering::ColumnOrdering>;

} // namespace algebra
