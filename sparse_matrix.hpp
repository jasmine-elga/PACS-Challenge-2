/**
 * @file sparse_matrix.hpp
 * @brief Contains the definition of Matrix class and its member functions.
 */

#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include "sparse_matrix_traits.hpp"
#include <iostream>
#include <iomanip> 
#include <map>
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include<utility>
#include<cstddef>
#include<stdexcept>
#include<numeric>
#include<cmath>
#include <fstream>
#include <sstream>
#include <random>
#include<complex>

namespace algebra {

    /**
     * @brief Enumerator indicating the Storage ordering
     * @param RowOrdering indicates ordering by rows
     * @param ColumnORdering indicates ordering by columns
     */
    enum class StorageOrder {RowOrdering, ColumnOrdering };


    /**
     * @brief Enumerator indicating the type of norm to be computed
     * @param One indicates the one-norm
     * @param Infinity indicates the infinity-norm
     * @param Frobenius indicates the Frobenius norm
    */
    enum class NormType{One, Infinity, Frobenius};

     
    /**
     * @brief Compares two complex numbers by their magnitudes.
     * 
     * @tparam T The type of elements in the complex numbers.
     * @param lhs The left-hand side complex number.
     * @param rhs The right-hand side complex number.
     * @return true if the magnitude of lhs is less than that of rhs, false otherwise.
     */
    template<typename T>
    bool complexLess(const std::complex<T>& lhs, const std::complex<T>& rhs) {
         return std::abs(lhs)<std::abs(rhs);
    }

    /**
     * @brief Helper struct for comparison based on storage ordering
     * 
     * @tparam Order The storage order
     */
    template <StorageOrder Order>
    struct CompareHelper {
        /**
         * @brief Compares two indices based on storage ordering
         * 
         * @param lhs Left-hand side index
         * @param rhs Right-hand side index
         * @return true if lhs should be ordered before rhs, false otherwise
         */
        bool operator()(const std::array<std::size_t, 2>& lhs, const std::array<std::size_t, 2>& rhs) const {
            if constexpr (Order == StorageOrder::RowOrdering) {
                return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]); 
            } else {
                return lhs[1] < rhs[1] || (lhs[1] == rhs[1] && lhs[0] < rhs[0]);
            }
        }
    };

    // Declaration of class Matrix (needed for the functions generateRandomVector and operator*)
    template<RealOrComplex T, StorageOrder Order >
    class Matrix;

    // Declaration of generateRandomVector (definition in the source file)
    template<RealOrComplex T, StorageOrder Order >
    std::vector<T> generateRandomVector(const Matrix<T, Order>& matrix);

    // Declaration and definition of operator* (matrix-vector multiplication)
    template<RealOrComplex T, StorageOrder Order >
    std::vector<T> operator*(const Matrix<T, Order>& matrix, const std::vector<T>& vec){
            std::vector<T> result(matrix.numrows, T{0}); // Initialize result vector with zeros
            if (!matrix.is_compressed()) {
                 // Uncompressed format
                std::size_t i{0};
                std::size_t j{0};
                
                // Traverse the matrix and compute dot product
                for (const auto& [coords, value] : matrix.uncompressed_data) {
                    i = coords[0];
                    j = coords[1];
                    result[i] += value * vec[j];
                }
            } else {
                // Compressed format
                if constexpr(Order == StorageOrder::RowOrdering){
                // Row ordering (CSR): traverse the matrix and perform classical row-times-vector algorithm
                    for (std::size_t i = 0; i < matrix.numrows; ++i) {
                        std::size_t row_start = matrix.compressed_inner[i];
                        std::size_t row_end = matrix.compressed_inner[i + 1];    

                        for (std::size_t k = row_start; k < row_end; ++k) {
                            std::size_t col_index = matrix.compressed_outer[k];
                            T value = matrix.compressed_data[k];
                            result[i] += value * vec[col_index]; // Compute dot product
                        }
                    }
                }
                else{
                    // Column ordering (CSC)
                    for (std::size_t j = 0; j < matrix.numcols; ++j) {
                        std::size_t col_start = matrix.compressed_inner[j];
                        std::size_t col_end = matrix.compressed_inner[j + 1];
                        
                        for (std::size_t k = col_start; k < col_end; ++k) {
                            std::size_t row_index = matrix.compressed_outer[k];
                            T value = matrix.compressed_data[k];
                            result[row_index] += value * vec[j]; // Compute linear combination
                        }
                    }  
                }
            }
            return result;
        }

    // Declaration and definition of operator* (matrix-matrix multiplication)
    template<RealOrComplex T, StorageOrder Order >
    std::vector<T> operator*(const Matrix<T, Order>& matrix, const Matrix<T, Order>& vec) {
            if (vec.numcols != 1) {
                throw std::invalid_argument("The second matrix must have one column.");
            }
            // Transform the matrix with only one column in a std::vector
            std::vector<T> vecColumn; // vector to store the column
            if(vec.is_compressed()){
                 for (size_t i = 0; i < vec.numrows; ++i){
                    vecColumn.push_back(vec.compressed_data[i]);
                 }
            }
            else{
                for (const auto& [coords, value] : vec.uncompressed_data){
                    vecColumn.push_back(value);
                    }
            }
            // Perform multiplication using the operator* overloaded for matrix-vector multiplication
            return matrix * vecColumn;
        }

    /**
     * @brief Sparse matrix class supporting compressed and uncompressed storage
     * 
     * @tparam T The type of elements in the matrix
     * @tparam Order The storage order of the matrix
     */
    template<RealOrComplex T, StorageOrder Order >
    class Matrix {
    private:
        bool compressed;  //!< indicates if the matrix is in compressed format.
        std::size_t numrows;  /*!< number of rows of the matrix*/
        std::size_t numcols; /*!< number of columns of the matrix*/
        
        //!< type of data in uncompressed state
        using MapData = std::map<std::array<std::size_t, 2>, T, CompareHelper<Order>>; 
        MapData uncompressed_data;  //!< stores data in uncompressed state
       
        // stores data in compressed state
        std::vector<std::size_t> compressed_inner; //!< stores inner index of compressed state
        std::vector<std::size_t> compressed_outer; //!< stores outer index of compressed state
        std::vector<T> compressed_data; //!< stores data of compressed state

    public:
        /**
         * @brief Constructor: constructs a new Sparse Matrix object
         * 
         * @param rows Number of rows in the matrix
         * @param cols Number of columns in the matrix
         */
        Matrix(std::size_t rows, std::size_t cols): compressed(false), numrows(rows), numcols(cols){};

        /**
         * @brief Provides non-const access to matrix elements, can add element only for uncompressed matrix
         * 
         * @param i Row index
         * @param j Column index
         * @return T& Reference to the element at (i, j)
         */
        T& operator()(std::size_t i, std::size_t j);

        /**
         * @brief Provides const access to matrix elements; returns 0 if the element is inside matrix bounds but is not present
         * 
         * @param i Row index
         * @param j Column index
         * @return const T& Const reference to the element at (i, j)
         */
        const T& operator()(std::size_t i, std::size_t j) const; //const version

    
        /**
         * @brief Compresses the matrix data
         */
        void compress();

        /**
         * @brief Uncompresses the matrix data
         */
        void uncompress(); 
        

        /**
         * @brief Utility: checks if the matrix is in compressed format
         * 
         * @return true if compressed, false otherwise
         */
        bool is_compressed() const{ return compressed;};

        /**
         * @brief Utility: Prints the matrix
         */
        void print() const;

        /**
         * @brief Utility: Provides access to matrix elements in compressed format
         * 
         * @param i Row index
         * @param j Column index
         * @return auto Reference to the element at (i, j) in compressed format
         */
        auto compressed_access(std::size_t i, std::size_t j);

        /**
         * @brief Utility: Provides access to matrix elements in compressed format
         * 
         * @param i Row index
         * @param j Column index
         * @return const auto Const reference to the element at (i, j) in compressed format
         */
        const auto compressed_access(std::size_t i, std::size_t j) const; // const version

        /**
         * @brief Utility: Resizes the matrix
         * 
         * @param rows New number of rows
         * @param cols New number of columns
         */
        void resize(std::size_t rows, std::size_t cols);
        
        /**
         * @brief Utility: Reads matrix data in from a file written in matrix market format
         * 
         * @param file_name Name of the file to read from
         */
        void  read(const std::string& file_name);
        

         /**
         * @brief Computes the norm of the matrix
         * 
         * @tparam N Type of norm to compute (One, Infinity, Frobenius)
         * @return T& Reference to the norm value
         */
        template<NormType N>
        T norm() const;
        


        // ##### FRIEND FUNCTIONS ####

        /**
         * @brief Overloaded operator for matrix-vector multiplication
         * 
         * @param matrix Matrix object
         * @param vec Vector to multiply with
         * @return std::vector<T> Resulting vector
         */
        friend std::vector<T> operator*<T,Order>(const Matrix<T, Order>& matrix, const std::vector<T>& vec);

        /**
         * @brief Overloaded operator for multiplication between a matrix and another matrix with just one column. 
         * Note: it works only if the two matrices are stored with the same order.
         * 
         * @param matrix Matrix object
         * @param vec Matrix with only one column to multiply with
         * @return std::vector<T> Resulting vector
         */
        friend std::vector<T> operator*<T,Order>(const Matrix<T, Order>& matrix, const Matrix<T, Order>& vec);



        /**
         * @brief Utility function to generate a random vector of suitable length to perform matrix-vector multiplication
         * 
         * @param matrix Matrix object
         * @return std::vector<T> Random vector
         */
        friend std::vector<T> generateRandomVector<T,Order>(const Matrix<T, Order>& matrix);




    };
}// namespace algebra



#endif // SPARSE_MATRIX_HPP















