## PACS Challenge 2 : A sparse matrix

This repositiory contains the code for the 2nd Challenge of the course Advanced Programming for Scientific Computing (PACS) at Politecnico di Milano.


# Description of the program
This program allows the user to work with sparse matrices stored in a compressed or an uncompressed format. The storage ordering can be row-wise or column-wise. The sparse matrix is implemented through a class.

In the **uncompressed** format, the matrix is stored in a `std::map`, where the couple (i,j) acts as the key. 

In the **compressed** format, the matrix is stored using three vectors: `compressed_inner`, `compressed_outer`, and `compressed_data`, representing the CSR (Compressed Sparse Row) or CSC (Compressed Sparse Column) format depending on the chosen storage order.

# Features
* Access elements: the non-const call operator can access elements and add them (in uncompressed format); whereas the const version returns 0 if elements is within the range of the matrix but not present.

* Compression and Uncompression: Convert between uncompressed and compressed formats.

* Matrix-vector multiplication: Operator*, extended also to the case where the vector is a matrix with just one column.

* Matrix Market File I/O: Read and write matrices in Matrix Market format from files.

* Norm Computations: Compute various norms (One norm, Infinity norm, Frobenius norm) for the matrix.

* Concepts and Traits: Utilizes C++ concepts and traits to handle numeric and complex types.




# How to install
Type

```
git clone git@github.com:jasmine-elga/PACS-Challenge-2.git
```

The Makefile contains the variable `PACS_ROOT`.
Before compiling the code, set it to the path of the course Examples as following:
```
export PACS_ROOT=/complete_path_to_pacs-examples/pacs-examples/Examples/ 
```

# How to compile and run the code

```
make
```

will compile the code. 


To run the code, type:

```
./main
```


In the main, 3 different test cases are implemented:
- A simple test case with a small matrix
- A more complex test case, where the matrix is read from the file `lnsp_131.mtx`
- A test case with a complex matrix 
