#ifndef SparseMatrixTraits_HPP
#define SparseMatrixTraits_HPP
#include<complex>
#include <iostream>
#include <type_traits>

namespace algebra{
//Concept for understand if a datum is numeric
template<typename T>
concept Numeric = std::is_arithmetic_v<T>;
//Concept to understand if a datum is complex
template<typename T>
concept Complex = requires(T t) {
    { t.real() } -> Numeric;
    { t.imag() } -> Numeric;
};

//The elements of the matrix can now be either a real or a complex value
template<typename T>
concept RealOrComplex = Numeric<T> || Complex<T>;
}
#endif /* SparseMatrixTraits_HPP */