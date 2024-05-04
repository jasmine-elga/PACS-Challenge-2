/**
 * @file sparse_matrix.cpp
 * @brief Contains the implementation of the Matrix memeber functions
 */

#include "sparse_matrix.hpp"


namespace algebra {

    // Utility function to access elements in compressed format
    template<RealOrComplex T, StorageOrder Order>
    auto Matrix<T, Order>::compressed_access(std::size_t i, std::size_t j) {
            // Look for the element
            //if j is in the interval [ outer[inner[i]], outer[inner[i+1]] ), then A[i,j] exists
            std::size_t row_start = compressed_inner[i];
            std::size_t row_end = compressed_inner[i + 1];
            auto it = std::find (compressed_outer.begin() + row_start,
                                 compressed_outer.begin() + row_end, j); 
            return it;
    }



     // Utility function to access elements in compressed format (const version)
    template<RealOrComplex T, StorageOrder Order>
    const auto Matrix<T, Order>::compressed_access(std::size_t i, std::size_t j) const{
            // Look for the element
            //if j is in the interval [ outer[inner[i]], outer[inner[i+1]] ), then A[i,j] exists
            std::size_t row_start = compressed_inner[i];
            std::size_t row_end = compressed_inner[i + 1];
            auto it = std::find (compressed_outer.begin() + row_start,
                                 compressed_outer.begin() + row_end, j); 
            return it;
    }



    // Call operator, non const version: can add elements if matrix in uncompressed format,
    // can only modify existing elements if in compressed format
    template<RealOrComplex T, StorageOrder Order>
    T& Matrix<T, Order>::operator()(std::size_t i, std::size_t j) {
        if (!is_compressed()) {
            // Uncompressed format, I can add new elements 
            if(i >=numrows || j >=numcols){
                resize(i,j); // Adding the element increases the size of the matrix
            }
            return uncompressed_data[{i,j}]; // Creates new data if element not found 
        }
        else{
            // Compressed format, I CANNOT add new elements
            if(i >=numrows || j >=numcols){
                // If position out of matrix bounds
                throw std::out_of_range("Matrix in compressed form, cannot add new elements!");
            }
            else{
                // Look for the element
                auto it = compressed_outer.begin();
                std::size_t row_end{0};   
                if constexpr(Order == StorageOrder::RowOrdering){
                    // Row ordering
                    it = compressed_access(i,j); // access to element ij in compress format (row ordering)
                    row_end = compressed_inner[i + 1];
                }
                else{
                    //column ordering
                    it = compressed_access(j,i); // access to element ij in compress format (column ordering)
                    row_end = compressed_inner[j + 1];
                    }
                
                if (it != compressed_outer.begin() + row_end) {
                    // If element is present
                    std::size_t index = std::distance(compressed_outer.begin(), it);
                    return compressed_data.at(index);  
                }
                else{
                    // If element is not present
                    throw std::out_of_range("Matrix in compressed form, cannot add new elements!");
                }
            }
        }
    }



    // Call operator, const version; returns 0 if the element is in matrix range but not present
    template<RealOrComplex T, StorageOrder Order>
    const T& Matrix<T, Order>::operator()(std::size_t i, std::size_t j) const {
        if (i >= numrows || j >= numcols) {
              throw std::out_of_range("Index out of boundary"); //if indexes out of range
        }

        // Uncompressed format
        if (!is_compressed()){
            auto it = uncompressed_data.find({i,j});
            // If we found the element
            if (it != uncompressed_data.end()) {
                return it->second;
            }
            // If we didn't find the element, return 0
            else{ 
                 static T default_value{};
                 return default_value; 
            }
        }
        // Compressed format
        else{
            auto it = compressed_outer.begin();
            std::size_t row_end{0};

            if constexpr(Order == StorageOrder::RowOrdering){
                // Row ordering
                it =compressed_access(i,j);
                row_end = compressed_inner[i + 1];
                }
            else{
                // Column ordering
                it = compressed_access(j,i);
                row_end = compressed_inner[j + 1];
                }          

        // If we found the element
        if (it != compressed_outer.begin() + row_end) {
            std::size_t index = std::distance(compressed_outer.begin(), it);
            return compressed_data.at(index);
            }

        //If we didn't find the element, return 0
        else{
            static T default_value{};
            return default_value; 
            }
        }
    }



    // Compresses an uncompressed matrix
    template<RealOrComplex T, StorageOrder Order>
    void Matrix<T, Order>::compress() {
       if (is_compressed()) {
            return; // Matrix is already compressed, no need to compress again
        }

        // Clear existing compressed data if any
        compressed_inner.clear();
        compressed_outer.clear();
        compressed_data.clear();
        std::size_t sz{0};

        if constexpr (Order == StorageOrder::RowOrdering) {
            sz = numrows;
        } else if constexpr (Order == StorageOrder::ColumnOrdering) {
            sz = numcols;
        }

        // reserve space
        compressed_inner.reserve(sz + 1); 
        compressed_outer.reserve(uncompressed_data.size()); 
        compressed_data.reserve(uncompressed_data.size()); 

        // Consider one row/column at a time
        std::size_t idx{0}; //counter for current row/column 
        std::size_t count{0}; //counter for current element
        auto start = uncompressed_data.begin();
        auto end = uncompressed_data.end();
        auto coords = start->first;
         
        while(idx!=sz){
            compressed_inner.push_back(count); 
            if constexpr(Order == StorageOrder::RowOrdering){
                // Consider row idx
                start = uncompressed_data.lower_bound({idx, 0});
                end = uncompressed_data.upper_bound({idx, std::numeric_limits<std::size_t>::max()});
            }
            else{
                // Consider column idx
                start = uncompressed_data.lower_bound({0, idx});
                end = uncompressed_data.upper_bound({std::numeric_limits<std::size_t>::max(), idx});
            }

            for (auto it = start ; it!= end; ++it){
                // Traverse the row/column
                coords = it->first;
                if constexpr(Order == StorageOrder::RowOrdering){
                    // Store the column index in the outer vector
                    compressed_outer.push_back(coords[1]); 
                }
                else{
                    // Store the row index in the outer vector
                    compressed_outer.push_back(coords[0]); 
                }
                // Store the value 
                compressed_data.push_back(it->second);  
                count++;
                }
            idx++; // increase current row/column
        }
        // For the last row/column: add fictitious index after the last valid index
        compressed_inner.push_back(count);

        // Clear uncompressed data after compression
        uncompressed_data.clear();
        compressed = true; // Set compressed flag
    }



    // Uncompresses a compressed matrix
    template<RealOrComplex T, StorageOrder Order>
    void Matrix<T, Order>::uncompress() {
       if (!is_compressed()) {
            return; // Matrix is already uncompressed, no need to uncompress again
        }
        uncompressed_data.clear();
        if constexpr(Order == StorageOrder::RowOrdering){
            // Compressed format (CSR)
            for (std::size_t i = 0; i < numrows; ++i) {
                // Traverse the matrix by rows
                std::size_t row_start = compressed_inner[i];
                std::size_t row_end = compressed_inner[i + 1];
                // Consider elements of row i                 
                for (std::size_t k = row_start; k < row_end; ++k) {
                        std::size_t col_index = compressed_outer[k];
                        T value = compressed_data[k];
                        uncompressed_data[{i, col_index}] = value;
                }
            }
        }
        else{
             // Compressed format (CSC)
             for (std::size_t j = 0; j < numcols; ++j) {
                // Traverse matrix by columns
                std::size_t col_start = compressed_inner[j];
                std::size_t col_end = compressed_inner[j + 1];
                // Consider elements of column j 
                for (std::size_t k = col_start; k < col_end; ++k) {
                        std::size_t row_index = compressed_outer[k];
                        T value = compressed_data[k];
                        uncompressed_data[{row_index, j}] = value;
                }
            }  
        }

        // clear compressed data and set compressed flag as false
        compressed_inner.clear();
        compressed_outer.clear();
        compressed_data.clear();
        compressed = false;
    }



    // Prints matrix of not too big dimensions
    // If the matrix is in uncompressed format, it renders the view and prints also zeros
    // If the matrix is in compressed format, it prints the 3 vectors (inner, outer, data)
    template<RealOrComplex T, StorageOrder Order>
    void Matrix<T, Order>::print() const {
            if (!is_compressed()) {
                // Print elements in uncompressed form
                std::cout << "Matrix (" << numrows << "x" << numcols << ") in non-compressed form:\n";
                if(numrows>20||numcols>20){
                    std::cerr<<"Matrix too big to be printed."<<std::endl;
                    return;
                    }
                for (std::size_t i = 0; i < numrows; ++i) {
                    for (std::size_t j = 0; j < numcols; ++j) {
                        T v = operator()(i,j);
                        std::cout<<v<< " ";
                    }
                    std::cout << '\n'; // Move to the next line after each row
                }
            } else {
                // Print elements in compressed format
                std::cout << "Matrix (" << numrows << "x" << numcols << ") in compressed form:\n";
                if(numrows>20||numcols>20){
                    std::cerr<<"Matrix too big to be printed."<<std::endl;
                    return;
                }
                std::cout << "Inner Index: ";
                for (size_t i = 0; i < compressed_inner.size(); ++i) {
                    std::cout << compressed_inner[i] << " ";
                }
                std::cout << std::endl;

                std::cout << "Outer Index: ";
                for (size_t i = 0; i < compressed_outer.size(); ++i) {
                    std::cout << compressed_outer[i] << " ";
                }
                std::cout << std::endl;

                std::cout << "Compressed Data: ";
                for (size_t i = 0; i < compressed_data.size(); ++i) {
                    std::cout << compressed_data[i] << " ";
                }
                std::cout << std::endl;
            }
        }
    

    
    
    // Resize method 
    template<RealOrComplex T, StorageOrder Order>
    void Matrix<T, Order>::resize(std::size_t rows, std::size_t cols) {
            if(!is_compressed()){
                numrows = rows;
                numcols = cols;
            }
    }



   
    // Reads matrix in matrix market format from a file
    template<RealOrComplex T, StorageOrder Order>
    void  Matrix<T, Order>::read(const std::string& file_name){
        std::ifstream file(file_name);
        if (!file.is_open()) {
            std::cerr << "Error, to open the file: " << file_name << std::endl;
            return;
        }

        std::string line;
        std::getline(file, line); // read the first line 
        
        // Check on the format of the first line
        if (line.substr(0, 14) != "%%MatrixMarket") {
            std::cerr << "Eror the file is not in format Matrix Market" << std::endl;
        }

        while (std::getline(file, line) && line[0] == '%') {
            // Ignore comments
        }

        // Read numbers of rows,columns and non zero elements
        std::size_t numNonZero;
        std::string numRowsStr, numColsStr, numNonZeroStr;
        std::istringstream iss(line);
        if(!(iss >> numRowsStr >> numColsStr >> numNonZeroStr)) {
            throw std::runtime_error("Error during the reading");
        }

        auto numRows = std::stoul(numRowsStr);
        auto numCols = std::stoul(numColsStr);
        numNonZero = std::stoul(numNonZeroStr);
        std::cout<<"Matrix read from file, in uncompressed format!"<<std::endl;
        std::cout <<"rows: "<< numRows<<", columns: "<< numCols<< ", non zero elements: "<< numNonZero<<std::endl;
        resize(numRows,numCols);

        // Read the non zero element
        while(std::getline(file,line)){
            std::istringstream elementStream(line);
            std::string nrow,ncol,val;
            elementStream >> nrow >> ncol >> val;
            std::size_t row,col;
            double value;
            row = std::stoul(nrow);
            col = std::stoul(ncol);
            value = std::stod(val);
            uncompressed_data[{row-1,col-1}]= static_cast<T>(value);
        }
        // Close the file
        file.close();
    }



    // Computes the norm of the matrix: options are One norm, Infinity norm and Frobenius norm
    template<RealOrComplex T, StorageOrder Order>
    template<NormType N>
    T Matrix<T, Order>::norm() const {
        T norm_value = 0;

        if(is_compressed()){
            // COMPRESSED format
            if constexpr (N == NormType::Frobenius) {
                // Frobenius norm computation for compressed matrix (same for row or column ordering)
                for (const auto& value : compressed_data) {
                    norm_value += std::abs(value) * std::abs(value);
                }
                norm_value = std::sqrt(norm_value);
            }
        
            else if constexpr (Order == StorageOrder::RowOrdering) {    
                // Row ordering in compressed format
                if constexpr (N == NormType::One){
                    // One - norm computation for compressed matrix, row ordering
                    std::vector<T> norms(numcols, 0); // vector to store the sums by column
                    for (std::size_t i = 0; i < numrows; ++i) {
                        for (std::size_t k = compressed_inner[i]; k < compressed_inner[i + 1]; ++k){
                            norms[compressed_outer[k]] += std::abs(compressed_data[k]); 
                        }
                    }
                    // Take the max of the sums by column
                    norm_value = *std::max_element(norms.begin(), norms.end(), complexLess<double>);
                }else if constexpr(N == NormType::Infinity){
                    // Infinity - norm computation for compressed matrix, row ordering
                    T row_sum{0}; //variable to store the sum over the row i
                    for (std::size_t i = 0; i < numrows; ++i) {
                        row_sum=0;
                        // Sum the elements of row i
                        for (std::size_t k = compressed_inner[i]; k < compressed_inner[i + 1]; ++k) {
                            row_sum += std::abs(compressed_data[k]);
                        }
                        // Take the max
                        norm_value = std::max(norm_value, row_sum, complexLess<double>);
                    }
                }
            }

            else if constexpr(Order == StorageOrder::ColumnOrdering){
                // Column ordering in compressed format
                if constexpr (N == NormType::One){
                    // One - norm computation for compressed matrix, column ordering
                    T col_sum{0};  //variable to store the sum over column j
                    for (std::size_t j = 0; j < numcols; ++j) {
                        col_sum = 0;
                        // Sum the elements of column j
                        for (std::size_t k = compressed_inner[j]; k < compressed_inner[j + 1]; ++k) {
                            col_sum += std::abs(compressed_data[k]);
                        }
                        // Take the max
                        norm_value = std::max(norm_value, col_sum, complexLess<double>);
                    }
                }
                else if constexpr(N == NormType::Infinity){
                    // Infinity - norm computation for compressed matrix, column ordering
                    std::vector<T> norms(numrows, 0); // vector to store the sums by row
                    for (std::size_t j = 0; j < numcols; ++j) {
                        for (std::size_t k = compressed_inner[j]; k < compressed_inner[j + 1]; ++k){
                            norms[compressed_outer[k]] += std::abs(compressed_data[k]); 
                        }
                    }
                    // Take max of the sums by row
                    norm_value = *std::max_element(norms.begin(), norms.end(), complexLess<double>);
                }
            }
        }else{
        // UNCOMPRESSED format
        if constexpr (N == NormType::Frobenius) {
            // Frobenius norm computation for compressed matrix (same for row or column ordering)
            for (const auto& [coords, value] : uncompressed_data) {
                    norm_value += std::abs(value) * std::abs(value);
            }
            norm_value = std::sqrt(norm_value);
        }

        else if constexpr (Order == StorageOrder::RowOrdering) {
            // Row ordering in uncompressed format
            if constexpr (N == NormType::One) {
                // One - norm computation for uncompressed matrix, row ordering
                std::vector<T> norms(numcols, 0); // vector to store the sums by column
                for (const auto& [coords, value] : uncompressed_data) {
                    norms[coords[1]] += std::abs(value);
                }
                // Take maximum
                norm_value = *std::max_element(norms.begin(), norms.end(),complexLess<double> );

            } else if constexpr (N == NormType::Infinity) {
                // Infinity-norm computation for uncompressed matrix, row ordering
                T row_sum{0};
                for (std::size_t i = 0; i < numrows; ++i) {
                    row_sum = 0;
                    // Take row i
                    auto start = uncompressed_data.lower_bound({i, 0});
                    auto end = uncompressed_data.upper_bound({i, std::numeric_limits<std::size_t>::max()});
                    for (auto it = start ; it!= end; ++it){
                    auto coords = it->first;
                    T value = it->second;
                    row_sum += std::abs(value);
                    }
                    norm_value = std::max(norm_value, row_sum, complexLess<double>);
                }
            }

        } else if constexpr (Order == StorageOrder::ColumnOrdering) {
            // Column ordering in uncompressed format
            if constexpr (N == NormType::One) {
                // One - norm computation for uncompressed matrix, column ordering
                for (std::size_t j = 0; j < numcols; ++j) {
                        T col_sum = 0;
                        auto start = uncompressed_data.lower_bound({0, j});
                        auto end = uncompressed_data.upper_bound({std::numeric_limits<std::size_t>::max(), j});
                        for (auto it = start ; it!= end; ++it){
                        auto coords = it->first;
                        T value = it->second;
                        col_sum += std::abs(value);
                        }
                        norm_value = std::max(norm_value, col_sum, complexLess<double>);
                }

            } else if constexpr (N == NormType::Infinity) {
                // Infinity-norm computation for uncompressed matrix, row ordering
                std::vector<T> norms(numrows, 0);
                for (const auto& [coords, value] : uncompressed_data) {
                    norms[coords[0]] += std::abs(value);
                }
                norm_value = *std::max_element(norms.begin(), norms.end(), complexLess<double>);
                } 
            }
        }
        return norm_value;
    }




    // Explicit instantiation for the types we use
    template class Matrix<double, StorageOrder::RowOrdering>;
    template class Matrix<double, StorageOrder::ColumnOrdering>;
    // Explicit instantiation for std::complex<double>
    template class algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering>;
    template class algebra::Matrix<std::complex<double>, algebra::StorageOrder::ColumnOrdering>;
    
    // Explicit instantiation for One-norm
    template double algebra::Matrix<double, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::One>() const;
    template double algebra::Matrix<double, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::One>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::One>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::One>() const;
    

    // Explicit instantiation for Infinity-norm
    template double algebra::Matrix<double, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::Infinity>() const;
    template double algebra::Matrix<double, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::Infinity>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::Infinity>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::Infinity>() const;

    // Explicit instantiation for Frobenius-norm
    template double algebra::Matrix<double, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::Frobenius>() const;
    template double algebra::Matrix<double, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::Frobenius>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering>::norm<algebra::NormType::Frobenius>() const;
    template std::complex<double> algebra::Matrix<std::complex<double>, algebra::StorageOrder::ColumnOrdering>::norm<algebra::NormType::Frobenius>() const;


} // namespace algebra
