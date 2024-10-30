#pragma once

#include <iostream>
#include <cstddef>
#include <vector>

namespace p {
/**
 * A simple square matrix class that supports as little operations as possible to be supported by
 * `Polynomial` class.
 *
 * @tparam T type of values in matrix
 * @tparam Size size of the matrix
 */
template <typename T, size_t Size>
class SquareMatrix {
public:
    SquareMatrix(const SquareMatrix<T, Size>& other);
    SquareMatrix(SquareMatrix<T, Size>&& other);
    SquareMatrix<T, Size>& operator=(const SquareMatrix<T, Size>& other);
    SquareMatrix<T, Size>& operator=(SquareMatrix<T, Size>&& other);

    SquareMatrix(size_t value = 0);
    SquareMatrix(std::vector<std::vector<T>> matrix);

    // Below is a minimal set of overloads that allow `SquareMatrix` class to be supported by
    // `Polynomial` class.

    SquareMatrix<T, Size>& operator+=(const SquareMatrix<T, Size>& rhs);
    SquareMatrix<T, Size>& operator-=(const SquareMatrix<T, Size>& rhs);
    SquareMatrix<T, Size>& operator*=(const SquareMatrix<T, Size>& rhs);

    friend SquareMatrix<T, Size> operator+(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
        lhs += rhs;
        return lhs;
    }

    friend SquareMatrix<T, Size> operator-(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend SquareMatrix<T, Size> operator*(SquareMatrix<T, Size> lhs, SquareMatrix<T, Size> rhs) {
        lhs *= rhs;
        return lhs;
    }

    // Debugging functions.

    const std::vector<T>& operator[](size_t idx) const;

    void Print() const;

    const std::vector<std::vector<T>>& GetMatrix() const;

private:
    std::vector<std::vector<T>> matrix_;
};

template <typename T, size_t Size>
SquareMatrix<T, Size>::SquareMatrix(const SquareMatrix<T, Size>& other) = default;

template <typename T, size_t Size>
SquareMatrix<T, Size>::SquareMatrix(SquareMatrix<T, Size>&& other) = default;

template <typename T, size_t Size>
SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator=(const SquareMatrix<T, Size>& other) =
    default;

template <typename T, size_t Size>
SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator=(SquareMatrix<T, Size>&& other) = default;

template <typename T, size_t Size>
SquareMatrix<T, Size>::SquareMatrix(size_t value) {
    matrix_.resize(Size);
    for (auto& row : matrix_) {
        row.resize(Size);
    }
    for (size_t i = 0; i < Size; ++i) {
        matrix_[i][i] = value;
    }
}

template <typename T, size_t Size>
SquareMatrix<T, Size>::SquareMatrix(std::vector<std::vector<T>> matrix) : matrix_(matrix) {
}

template <typename T, size_t Size>
SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator+=(const SquareMatrix<T, Size>& rhs) {
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = 0; j < Size; ++j) {
            matrix_[i][j] += rhs.matrix_[i][j];
        }
    }
    return *this;
}

template <typename T, size_t Size>
SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator-=(const SquareMatrix<T, Size>& rhs) {
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = 0; j < Size; ++j) {
            matrix_[i][j] -= rhs.matrix_[i][j];
        }
    }
    return *this;
}

template <typename T, size_t Size>
SquareMatrix<T, Size>& SquareMatrix<T, Size>::operator*=(const SquareMatrix<T, Size>& rhs) {
    std::vector<std::vector<T>> result(Size, std::vector<T>(Size));
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = 0; j < Size; ++j) {
            for (size_t k = 0; k < Size; ++k) {
                result[i][j] += matrix_[i][k] * rhs.matrix_[k][j];
            }
        }
    }
    matrix_ = result;
    return *this;
}

template <typename T, size_t Size>
const std::vector<T>& SquareMatrix<T, Size>::operator[](size_t idx) const {
    return matrix_[idx];
}

template <typename T, size_t Size>
void SquareMatrix<T, Size>::Print() const {
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = 0; j < Size; ++j) {
            std::cout << matrix_[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

template <typename T, size_t Size>
const std::vector<std::vector<T>>& SquareMatrix<T, Size>::GetMatrix() const {
    return matrix_;
}
}  // namespace p
