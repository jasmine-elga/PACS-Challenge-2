#include <iostream>
#include "sparse_matrix.hpp"
#include <chrono>



/*

int main() {
    
    std::cout<<"## Test with a sparse matrix stored in ROW ordering ##"<<std::endl;
    
    //initialization of the matrix
    algebra::Matrix<double,algebra::StorageOrder::RowOrdering> M1(0,0);
    
    //reading the matrix from a file
    M1.read("lnsp_131.mtx");
    std::cout<<"\nPrinting the matrix..."<<std::endl;
    M1.print();
    
    // creating a random vector to perform multiplication
    std::vector<double> randomV = generateRandomVector(M1);

    // perform multiplication
    auto t0_1 = std::chrono::high_resolution_clock::now();
    std::vector<double> res_uncompressed = M1*randomV;
    auto t1_1 = std::chrono::high_resolution_clock::now();
    auto delta_t_1 = std::chrono::duration_cast<std::chrono::microseconds>(t1_1-t0_1);  
    std::cout<<"\nM*v, UNCOMPRESSED format (row ordering) "<<std::endl;
    std::cout << ("The operation takes") << " : " << delta_t_1.count() << " ms" << std::endl;

    std::cout<< "\nLet's now compress the matrix!"<<std::endl;
    M1.compress();
    
    std::cout<<"\nPrinting the matrix..."<<std::endl;
    M1.print();
    
    // perform multiplication
    auto t0_2 = std::chrono::high_resolution_clock::now();
    std::vector<double> res_compressed = M1*randomV;
    auto t1_2 = std::chrono::high_resolution_clock::now();
    auto delta_t_2 = std::chrono::duration_cast<std::chrono::microseconds>(t1_2-t0_2);  
    std::cout<<"\nM*v, COMPRESSED format (row ordering)"<<std::endl;
    std::cout << ("The operation takes") << " : " << delta_t_2.count() << " ms" << std::endl;




    std::cout<<"\n\n\n## Test with a sparse matrix stored in COLUMN ordering ##"<<std::endl;
    
    //initialization of the matrix
    algebra::Matrix<double,algebra::StorageOrder::ColumnOrdering> M2(0,0);
    
    //reading the matrix from a file
    M2.read("lnsp_131.mtx");
    std::cout<<"\nPrinting the matrix..."<<std::endl;
    M2.print();
    
    // perform multiplication
    auto t0_3 = std::chrono::high_resolution_clock::now();
    res_uncompressed = M2*randomV;
    auto t1_3 = std::chrono::high_resolution_clock::now();
    auto delta_t_3 = std::chrono::duration_cast<std::chrono::microseconds>(t1_3-t0_3);  
    std::cout<<"\nM*v, UNCOMPRESSED format (column ordering)"<<std::endl;
    std::cout << ("The operation takes") << " : " << delta_t_3.count() << " ms" << std::endl;

    std::cout<< "\nLet's now compress the matrix!"<<std::endl;
    M2.compress();
    
    std::cout<<"\nPrinting the matrix..."<<std::endl;
    M2.print();
    
    // perform multiplication
    auto t0_4 = std::chrono::high_resolution_clock::now();
    res_compressed = M2*randomV;
    auto t1_4 = std::chrono::high_resolution_clock::now();
    auto delta_t_4 = std::chrono::duration_cast<std::chrono::microseconds>(t1_4-t0_4);  
    std::cout<<"\nM*v, COMPRESSED format (column ordering)"<<std::endl;
    std::cout << ("The operation takes") << " : " << delta_t_4.count() << " ms" << std::endl;





    std::cout<<"\n\n\n## Smaller test case ##"<<std::endl;
    algebra::Matrix<double,algebra::StorageOrder::ColumnOrdering> A(5,3);

    A(0,0)=1;
    A(0,2)=3;
    A(1,0)=4;
    A(1,1)=5+1*i;
    A(2,1)=8;
    A(2,2)=6;
    A(3,1)=1;
    A(4,0) =2; 

    A.print();

    // Compute norms
    std::cout << "\nOne-Norm of A (Uncompressed): " << A.norm<algebra::NormType::One>() << std::endl;
    std::cout << "Infinity-Norm of A (Uncompressed): " << A.norm<algebra::NormType::Infinity>() << std::endl;
    std::cout << "Frobenius-Norm of A (Uncompressed): " << A.norm<algebra::NormType::Frobenius>() << std::endl;

    std::vector<double> v = {1,2,3};

    std::cout<<"\nProduct matrix vector in uncompressed format:  "<<std::endl;
    std::cout<<"v = "<<"["<<v[0]<<" "<<v[1]<<" "<<v[2]<<"]"<<std::endl;
    auto t0_1 = std::chrono::high_resolution_clock::now();
    std::vector<double> res_uncompressed = A*v;
    auto t1_1 = std::chrono::high_resolution_clock::now();
    auto delta_t = std::chrono::duration_cast<std::chrono::microseconds>(t1_1-t0_1);
    std::cout<<"A*v, uncompressed format"<<std::endl;
    std::cout<<"Result: "<<"["<<res_uncompressed[0]<<" "<<res_uncompressed[1]<<" "<<res_uncompressed[2]<<"]"<<std::endl;
    std::cout << "The operation takes : " << delta_t.count() << " ms" << std::endl;


    std::cout<<"\n\nNow let us compress the matrix!\n\n"<<std::endl;
    A.compress();
    A.print();

    // Compute norms
    std::cout << "\nOne-Norm of A (Compressed): " << A.norm<algebra::NormType::One>() << std::endl;
    std::cout << "Infinity-Norm of A (Compressed): " << A.norm<algebra::NormType::Infinity>() << std::endl;
    std::cout << "Frobenius-Norm of A (Compressed): " << A.norm<algebra::NormType::Frobenius>() << std::endl;

    //A(1,1) =9;
    std::cout<<"\nProduct matrix vector in compressed form: "<<std::endl;
    auto t0_2 = std::chrono::high_resolution_clock::now();
    std::vector<double> res_compressed = A*v;
    auto t1_2 = std::chrono::high_resolution_clock::now();
    auto delta_t_2 = std::chrono::duration_cast<std::chrono::microseconds>(t1_2-t0_2);
  

    std::cout<<"A*v, compressed format"<<std::endl;
    std::cout<<"Result: "<<"["<<res_compressed[0]<<" "<<res_compressed[1]<<" "<<res_compressed[2]<<"]"<<std::endl;
    std::cout << "The operation takes : " << delta_t_2.count() << " ms" << std::endl;

    //A.uncompress();
    //A.print();

    algebra::Matrix<double,algebra::StorageOrder::ColumnOrdering> B(3,1);

    B(0,0) = 1;
    B(1,0) = 2;
    B(2,0) = 3;
    B.print();

    if(A.is_compressed()){
        B.compress();
    }

    std::cout<<"\nProduct matrix - matrix with only one row in compressed form: "<<std::endl;
    auto t0_3 = std::chrono::high_resolution_clock::now();
    std::vector<double> res_new = A*B;
    auto t1_3 = std::chrono::high_resolution_clock::now();
    auto delta_t_3 = std::chrono::duration_cast<std::chrono::microseconds>(t1_3-t0_3);
  

    std::cout<<"A*v, compressed format"<<std::endl;
    std::cout<<"Result: "<<"["<<res_new[0]<<" "<<res_new[1]<<" "<<res_new[2]<<"]"<<std::endl;
    std::cout << "The operation takes : " << delta_t_3.count() << " ms" << std::endl;


    return 0; 
} */


int main() {
    // Define a sparse matrix of complex numbers
    algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering> complexMatrix(3, 3);
    

    // Fill the matrix with complex values
    complexMatrix(0, 0) = {1.0, 2.0};
    complexMatrix(1, 1) = {3.0, 4.0};
    complexMatrix(2, 2) = {5.0, 6.0};

    complexMatrix.print();
    
    // Define a vector of complex numbers
    algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowOrdering> complexVec(3, 1);
    complexVec(0, 0) = {1.0, 1.0};
    complexVec(1, 0) ={2.0, 2.0};
    complexVec(2, 0) = {3.0, 3.0};

    
    // Perform matrix-vector multiplication
    auto result = complexMatrix * complexVec;
    
    // Print the result vector
    std::cout << "Resulting vector:" << std::endl;
    for (const auto& elem : result) {
        std::cout << elem << std::endl;
    }


    complexMatrix.compress();
    // Compute norms
    std::cout << "\nOne-Norm of complexMatrix (Uncompressed): " << complexMatrix.norm<algebra::NormType::One>() << std::endl;
    std::cout << "Infinity-Norm of complexMatrix (Uncompressed): " << complexMatrix.norm<algebra::NormType::Infinity>() << std::endl;
    std::cout << "Frobenius-Norm of complexMatrix (Uncompressed): " << complexMatrix.norm<algebra::NormType::Frobenius>() << std::endl;

    
    return 0;
}