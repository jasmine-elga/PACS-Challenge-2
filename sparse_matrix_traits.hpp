/**
 * @file sparse_matrix_traits.hpp
 * @brief Contains traits and concepts for sparse matrix elements.
 */


#ifndef SparseMatrixTraits_HPP
#define SparseMatrixTraits_HPP
#include<complex>
#include <iostream>
#include <type_traits>

namespace algebra{
    /**
     * @brief Concept to check if a type is numeric.
     * @tparam T Type to check.
     */
    template<typename T>
    concept Numeric = std::is_arithmetic_v<T>;


    /**
     * @brief Concept to check if a type is complex.
     * @tparam T Type to check.
     */
    template<typename T>
    concept Complex = requires(T t) {
        { t.real() } -> Numeric;
        { t.imag() } -> Numeric;
    };

    /**
     * @brief Concept to check if a type is either real or complex.
     * @tparam T Type to check.
     */
    template<typename T>
    concept RealOrComplex = Numeric<T> || Complex<T>;

} // namespace algebra

#endif /* SparseMatrixTraits_HPP */