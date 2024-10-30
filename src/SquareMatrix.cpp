#include "../include/SquareMatrix.hpp"

// namespace p {
// template <typename T, size_t Size>
// SquareMatrix<T, Size>::SquareMatrix(const SquareMatrix<T, Size>& other) = default;
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>::SquareMatrix(SquareMatrix<T, Size>&& other) = default;
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator=(const SquareMatrix<T, Size>& other) =
//     default;
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator=(SquareMatrix<T, Size>&& other) = default;
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>::SquareMatrix(size_t value) {
//     matrix_.resize(Size);
//     for (auto& row : matrix_) {
//         row.resize(Size);
//     }
//     for (size_t i = 0; i < Size; ++i) {
//         matrix_[i][i] = value;
//     }
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>::SquareMatrix(std::vector<std::vector<T>> matrix) : matrix_(matrix) {
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator+=(const SquareMatrix<T, Size>& rhs) {
//     for (size_t i = 0; i < Size; ++i) {
//         for (size_t j = 0; j < Size; ++j) {
//             matrix_[i][j] += rhs.matrix_[i][j];
//         }
//     }
//     return *this;
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator-=(const SquareMatrix<T, Size>& rhs) {
//     for (size_t i = 0; i < Size; ++i) {
//         for (size_t j = 0; j < Size; ++j) {
//             matrix_[i][j] -= rhs.matrix_[i][j];
//         }
//     }
//     return *this;
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator*=(const SquareMatrix<T, Size>& rhs) {
//     std::vector<std::vector<T>> result(Size, std::vector<T>(Size));
//     for (size_t i = 0; i < Size; ++i) {
//         for (size_t j = 0; j < Size; ++j) {
//             for (size_t k = 0; k < Size; ++k) {
//                 result[i][j] += matrix_[i][k] * rhs.matrix_[k][j];
//             }
//         }
//     }
//     matrix_ = result;
//     return *this;
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size> operator+(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
//     lhs += rhs;
//     return lhs;
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size> operator-(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
//     lhs -= rhs;
//     return lhs;
// }
//
// template <typename T, size_t Size>
// SquareMatrix<T, Size> operator*(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
//     lhs *= rhs;
//     return lhs;
// }
//
// template <typename T, size_t Size>
// void SquareMatrix<T, Size>::Print() const {
//     for (size_t i = 0; i < Size; ++i) {
//         for (size_t j = 0; j < Size; ++j) {
//             std::cout << matrix_[i][j] << ' ';
//         }
//         std::cout << std::endl;
//     }
// }
// }  // namespace p
