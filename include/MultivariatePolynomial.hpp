#pragma once

#include <array>
#include <cxxabi.h>
#include <map>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

namespace p {
/* A class that stores a multivariate polynomial.
 *
 * @tparam T type of the variable
 * WARNING: same constraints as for `Polynomial` class must apply to T.
 */
template <typename T, size_t N>
class MVPolynomial {
public:
    // Using rule of zero (explicitly) as no special/unusual constructing is required.

    MVPolynomial();
    MVPolynomial(const MVPolynomial<T, N>&);
    MVPolynomial(MVPolynomial<T, N>&&);
    MVPolynomial<T, N>& operator=(const MVPolynomial<T, N>&);
    MVPolynomial<T, N>& operator=(MVPolynomial<T, N>&&);

    // Constructs a `MVPolynomial` object from a vector of pairs {powers_array, coefficient}.
    explicit MVPolynomial(std::vector<std::pair<std::array<T, N>, T>> list);

    // Mathematical operator overloads

    MVPolynomial<T, N> operator-();
    MVPolynomial<T, N>& operator+=(const MVPolynomial<T, N>& rhs);
    MVPolynomial<T, N>& operator-=(const MVPolynomial<T, N>& rhs);
    MVPolynomial<T, N>& operator*=(const MVPolynomial<T, N>& rhs);

    // Static methods that return predetermined multivariate polynomials that don't rely on any
    // other multivariate polynomial.

    // Returns a multiplicative identity multivariate polynomial.
    static MVPolynomial<T, N> Identity();
    // Returns a constant multivariate polynomial in form `value * x_1^0 * ... * x_n^0`.
    static MVPolynomial<T, N> Constant(const T& value);
    // Returns an additive identity multivariate polynomial.
    static MVPolynomial<T, N> Zero();

    friend MVPolynomial<T, N> operator+(MVPolynomial<T, N> lhs, MVPolynomial<T, N> rhs) {
        lhs += rhs;
        return lhs;
    }

    friend MVPolynomial<T, N> operator-(MVPolynomial<T, N> lhs, MVPolynomial<T, N> rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend MVPolynomial<T, N> operator*(MVPolynomial<T, N> lhs, MVPolynomial<T, N> rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend MVPolynomial<T, N> operator+(T value, MVPolynomial<T, N> rhs) {
        std::array<T, N> constant_coefficient = {0};
        rhs.coefficients_[constant_coefficient] += value;
        return rhs;
    }

    friend MVPolynomial<T, N> operator+(MVPolynomial<T, N> lhs, T value) {
        return value + lhs;
    }

    friend MVPolynomial<T, N> operator-(MVPolynomial<T, N> lhs, T value) {
        std::array<T, N> constant_coefficient = {0};
        if (lhs.coefficients_[constant_coefficient] < value) {
            lhs.UnsignednessCheck();
        }
        lhs.coefficients_[constant_coefficient] -= value;
        return lhs;
    }

    friend MVPolynomial<T, N> operator-(T value, MVPolynomial<T, N> rhs) {
        return -(rhs - value);
    }

    friend MVPolynomial<T, N> operator*(T value, MVPolynomial<T, N> rhs) {
        for (auto& [powers, coeff] : rhs.coefficients_) {
            rhs.coefficients_[powers] = value * rhs.coefficients_[powers];
        }
        return rhs;
    }

    friend MVPolynomial<T, N> operator*(MVPolynomial<T, N> lhs, T value) {
        for (auto& [powers, coeff] : lhs.coefficients_) {
            lhs.coefficients_[powers] *= value;
        }
        return lhs;
    }

    // Other useful functions.

    /**
     * Raise `MVPolynomial` object to `power`th power. For user's safety, this is not done inplace.
     *
     * @param power The power `MVPolynomial` is raised to.
     * @return resulting `MVPolynomial` object.
     */
    MVPolynomial<T, N> Power(size_t power) const;

    /**
     * Evaluate multivariate polynomial at point x = (x_0, ..., x_N).
     *
     * @param point The point polynomial is being evaluated at.
     * @return resulting value.
     */
    T At(const std::array<T, N>& point) const;

    /**
     * Evaluate multivariate polynomial at point x = (x_0, ..., x_N) (copy of At() function).
     *
     * @param point The point polynomial is being evaluated at.
     * @return resulting value.
     */
    T operator()(const std::array<T, N>& point);

    /**
     * Substitute another multivariate polynomials in the current one. For user's safety, this is
     * not done inplace.
     *
     * @param other Multivariate polynomials that are being substituted into the current one. Given
     * as an std::array of N multivariate polynomials, where other[i] is the multivariate polynomial
     * that will be substituted for x_i in the current polynomial.
     * @return new multivariate polynomial.
     */
    MVPolynomial<T, N> Substitute(const std::array<MVPolynomial<T, N>, N>& other) const;

    /**
     * Returns the coefficients of the multivariate polynomial as they are stored.
     *
     * @return coefficients in form of std::map<std::array<T, N>, T>.
     */
    const std::map<std::array<T, N>, T>& GetCoefficients() const;

private:
    std::map<std::array<T, N>, T> coefficients_;

    // Automatic unsigned overflow/underflow check.
    void UnsignednessCheck();

    // Utility function for exponentiating `value` to `power`.
    T BinPower(const T& value, size_t power) const;
};

template <typename T, size_t N>
MVPolynomial<T, N>::MVPolynomial() = default;

template <typename T, size_t N>
MVPolynomial<T, N>::MVPolynomial(const MVPolynomial<T, N>&) = default;

template <typename T, size_t N>
MVPolynomial<T, N>::MVPolynomial(MVPolynomial<T, N>&&) = default;

template <typename T, size_t N>
MVPolynomial<T, N>& MVPolynomial<T, N>::operator=(const MVPolynomial<T, N>&) = default;

template <typename T, size_t N>
MVPolynomial<T, N>& MVPolynomial<T, N>::operator=(MVPolynomial<T, N>&&) = default;

template <typename T, size_t N>
MVPolynomial<T, N>::MVPolynomial(std::vector<std::pair<std::array<T, N>, T>> list) {
    for (const auto& pair : list) {
        coefficients_.insert(pair);
    }
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::operator-() {
    UnsignednessCheck();
    std::vector<std::pair<std::array<T, N>, T>> coeffs;
    for (const auto& [powers, coeff] : coefficients_) {
        coeffs.push_back({powers, -coeff});
    }
    return MVPolynomial<T, N>(coeffs);
}

template <typename T, size_t N>
MVPolynomial<T, N>& MVPolynomial<T, N>::operator+=(const MVPolynomial<T, N>& rhs) {
    for (const auto& [powers, coeff] : rhs.coefficients_) {
        coefficients_[powers] += coeff;
    }
    return *this;
}

template <typename T, size_t N>
MVPolynomial<T, N>& MVPolynomial<T, N>::operator-=(const MVPolynomial<T, N>& rhs) {
    for (const auto& [powers, coeff] : rhs.coefficients_) {
        coefficients_[powers] -= coeff;
    }
    return *this;
}

template <typename T, size_t N>
MVPolynomial<T, N>& MVPolynomial<T, N>::operator*=(const MVPolynomial<T, N>& rhs) {
    std::map<std::array<T, N>, T> new_coefficients;
    for (const auto& m1 : coefficients_) {
        for (const auto& m2 : rhs.coefficients_) {
            std::array<T, N> new_powers;
            for (size_t i = 0; i < N; ++i) {
                new_powers[i] = m1.first[i] + m2.first[i];
            }
            new_coefficients[new_powers] += m1.second * m2.second;
        }
    }
    coefficients_ = new_coefficients;
    return *this;
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::Identity() {
    std::array<T, N> constant_coefficient = {0};
    return MVPolynomial<T, N>({{constant_coefficient, 1}});
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::Constant(const T& value) {
    std::array<T, N> constant_coefficient = {0};
    return MVPolynomial<T, N>({{constant_coefficient, value}});
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::Zero() {
    std::array<T, N> constant_coefficient = {0};
    return MVPolynomial<T, N>({{constant_coefficient, 0}});
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::Power(size_t power) const {
    if (power == 0) {
        return MVPolynomial<T, N>::Identity();
    }
    MVPolynomial<T, N> result = *this;
    for (size_t i = 1; i < power; ++i) {
        result *= *this;
    }
    return result;
}

template <typename T, size_t N>
T MVPolynomial<T, N>::At(const std::array<T, N>& point) const {
    T result = 0;
    for (const auto& [powers, coeff] : coefficients_) {
        T intermediate = coeff;
        for (size_t i = 0; i < N; ++i) {
            intermediate *= BinPower(point[i], powers[i]);
        }
        result += intermediate;
    }
    return result;
}

template <typename T, size_t N>
T MVPolynomial<T, N>::operator()(const std::array<T, N>& point) {
    return At(point);
}

template <typename T, size_t N>
MVPolynomial<T, N> MVPolynomial<T, N>::Substitute(
    const std::array<MVPolynomial<T, N>, N>& other) const {
    MVPolynomial<T, N> result;
    for (const auto& [powers, coeff] : coefficients_) {
        MVPolynomial<T, N> intermediate = coeff;
        for (size_t i = 0; i < N; ++i) {
            intermediate *= other[i].Power(powers[i]);
        }
        result += intermediate;
    }
    return result;
}

template <typename T, size_t N>
const std::map<std::array<T, N>, T>& MVPolynomial<T, N>::GetCoefficients() const {
    return coefficients_;
}

template <typename T, size_t N>
void MVPolynomial<T, N>::UnsignednessCheck() {
    if constexpr (std::is_unsigned_v<T>) {
        int status = 0;
        std::string T_type = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
        throw std::runtime_error(
            "An operation on Polynomial<" + (status == 0 ? T_type : typeid(T).name()) +
            "> of unsigned type requires one of the operands to become negative.");
    }
}

template <typename T, size_t N>
T MVPolynomial<T, N>::BinPower(const T& value, size_t power) const {
    if (power == 0) {
        return 1;
    }
    T result = BinPower(value, power / 2);
    if (power % 2 == 1) {
        return result * result * value;
    } else {
        return result * result;
    }
}
}  // namespace p
