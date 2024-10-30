#pragma once

#include <cxxabi.h>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace p {
// Forward declaration, required for `friend class Polynomial<T>;` to work in Variable class.
template <typename T>
class Polynomial;

template <typename T, size_t Size>
class SquareMatrix;

/**
 * `Variable` class allows user to construct a `Polynomial` object in a way that us available in
some calculation-related Python libraries: a user creates an object `Variable<T> x("x")` and then
writes a human-readable mathematical expression of the polynomial in form `Polynomial<int> p = x^2 +
1`.
 * WARNING: due to the fact that basically any basic mathematical operation has a higher priority
than `^`, user must put the exponentiated variable in braces.
 *
 * @tparam T type of the variable
 */
template <typename T>
class Variable {
public:
    Polynomial<T> operator^(size_t power);

    // A workaround that allows to write something like `2 * x` instead of `2 * (x^1)`.
    template <typename U>
    friend Polynomial<U> operator*(U value, Variable<U> v);

    // Respectively, a workaround that allows user to write something like `x * 2` instead of `(x^1)
    // * 2`.
    template <typename U>
    friend Polynomial<U> operator*(Variable<U> v, U value);

private:
    template <typename U>
    friend class Polynomial;
};

/** A class that stores the polynomial.
 *
 * @tparam T type of the variable.
 * WARNING: constraints for T templated type:
 *  - T must support addition, subtraction and multiplication, hence respective overloaded
 * operators;
 *  - T must have an ability to be *implicitly* constructed from 0 (additive identity) and 1
 * (multiplicative identity).
 */
template <typename T>
class Polynomial {
public:
    // Using rule of zero (explicitly) as no special/unusual constructing is required.

    Polynomial();
    Polynomial(const Polynomial<T>&);
    Polynomial(Polynomial<T>&&);
    Polynomial<T>& operator=(const Polynomial<T>&);
    Polynomial<T>& operator=(Polynomial<T>&&);

    /**
     * Construct a `Polynomial` object from a vector of terms.
     *
     * @param list A vector of terms, where on $i$-th position stands a coefficient before x^i.
     */
    Polynomial(std::vector<T> list);

    // Mathematical operator overloads.

    Polynomial<T> operator-();
    Polynomial<T>& operator+=(const Polynomial<T>& rhs);
    Polynomial<T>& operator-=(const Polynomial<T>& rhs);
    Polynomial<T>& operator*=(const Polynomial<T>& rhs);

    // Static methods that return predetermined polynomials that don't rely on any other polynomial.

    /**
     * Returns a multiplicative identity polynomial.
     */
    static Polynomial<T> Identity();
    /**
     * Returns a polynomial in form `value * x^0`.
     */
    static Polynomial<T> Constant(const T& value);
    /**
     * Returns an additive identity polynomial.
     */
    static Polynomial<T> Zero();

    // Friend operator overloads
    // NOTE: in the following implementations it is considered that arithmetical operations on T
    // type are non-commutative (which is indeed true for matrices), hence a lot of seemingly
    // unnecessary copy and pasted code.

    friend Polynomial<T> operator+(Polynomial<T> lhs, Polynomial<T> rhs) {
        lhs += rhs;
        return lhs;
    }

    friend Polynomial<T> operator-(Polynomial<T> lhs, Polynomial<T> rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend Polynomial<T> operator*(Polynomial<T> lhs, Polynomial<T> rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Polynomial<T> operator+(T value, Polynomial<T> rhs) {
        rhs.coefficients_[0] += value;
        return rhs;
    }

    friend Polynomial<T> operator+(Polynomial<T> lhs, T value) {
        return value + lhs;
    }

    friend Polynomial<T> operator-(Polynomial<T> lhs, T value) {
        if (lhs.coefficients_[0] < value) {
            lhs.UnsignednessCheck();
        }
        lhs.coefficients_[0] -= value;
        return lhs;
    }

    friend Polynomial<T> operator-(T value, Polynomial<T> rhs) {
        return -(rhs - value);
    }

    friend Polynomial<T> operator*(T value, Polynomial<T> rhs) {
        for (size_t i = 0; i < rhs.coefficients_.size(); ++i) {
            rhs.coefficients_[i] = value * rhs.coefficients_[i];
        }
        return rhs;
    }

    friend Polynomial<T> operator*(Polynomial<T> lhs, T value) {
        for (size_t i = 0; i < lhs.coefficients_.size(); ++i) {
            lhs.coefficients_[i] *= value;
        }
        return lhs;
    }

    // Other useful functions.

    /**
     * Raises the polynomial to `power`th power. For user's safety exponentiation is not done
     * inplace.
     *
     * @param power Power the polynomial is raised to.
     * @return the resulting polynomial.
     */
    Polynomial<T> Power(size_t power) const;

    /**
     * Evaluates polynomial at point `point`.
     *
     * @param point A point the polynomial is begin evaluated at.
     * @return value at point `point`.
     */
    const T At(const T& point) const;

    /**
     * Evaluates polynomial at point `point` (copy of `At`).
     *
     * @param point A point the polynomial is begin evaluated at.
     * @return value at point `point`.
     */
    const T operator()(const T& point) const;

    /**
     * Substitutes polynomial `other` into the current polynomial. For user's safety substitution is
     * not done inplace.
     *
     * @param other The polynomial that is being substituted.
     * @return the resulting polynomial.
     */
    Polynomial<T> Substitute(const Polynomial<T>& other) const;

    /**
     * Returns coefficients of the polynomial as they are stored inside of the class.
     *
     * @return a map of coefficients.
     */
    const std::map<size_t, T>& GetCoefficients() const;

private:
    std::map<size_t, T> coefficients_;

    // Allows to construct a monomial when using Variable object.
    template <typename U>
    friend class Variable;

    // This automatic (i.e. user should not have an ability to call it) check is quite obviously
    // only done for explicitly unsigned types.
    void UnsignednessCheck();

    // Utility function for (rather quickly) exponentiating `value` to `power` (not to overload
    // standard `std::pow`).
    T BinPower(const T& value, size_t power) const;
};

template <typename T>
Polynomial<T> Variable<T>::operator^(size_t power) {
    Polynomial<T> result;
    result.coefficients_[power] = 1;
    return result;
}

template <typename T>
Polynomial<T> operator*(T value, Variable<T> v) {
    return value * (v ^ 1);
}

template <typename T>
Polynomial<T> operator*(Variable<T> v, T value) {
    return (v ^ 1) * value;
}

template <typename T>
Polynomial<T>::Polynomial() = default;

template <typename T>
Polynomial<T>::Polynomial(const Polynomial<T>&) = default;

template <typename T>
Polynomial<T>::Polynomial(Polynomial<T>&&) = default;

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>&) = default;

template <typename T>
Polynomial<T>& Polynomial<T>::operator=(Polynomial<T>&&) = default;

template <typename T>
Polynomial<T>::Polynomial(std::vector<T> list) {
    size_t i = 0;
    for (const T& c : list) {
        coefficients_[i] = c;
        ++i;
    }
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-() {
    UnsignednessCheck();
    std::vector<T> coeff;
    for (const auto& [key, value] : coefficients_) {
        coeff.resize(key + 1);
        coeff[key] = -value;
    }
    return Polynomial<T>(coeff);
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& rhs) {
    for (const auto& [power, c] : rhs.coefficients_) {
        coefficients_[power] += c;
    }
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& rhs) {
    for (const auto& [power, c] : rhs.coefficients_) {
        coefficients_[power] -= c;
    }
    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& rhs) {
    std::map<size_t, T> new_coefficients;
    for (const auto& m1 : coefficients_) {
        for (const auto& m2 : rhs.coefficients_) {
            new_coefficients[m1.first + m2.first] += m1.second * m2.second;
        }
    }
    coefficients_ = new_coefficients;
    return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::Identity() {
    return Polynomial<T>({1});
}

template <typename T>
Polynomial<T> Polynomial<T>::Constant(const T& value) {
    return Polynomial<T>({value});
}

template <typename T>
Polynomial<T> Polynomial<T>::Zero() {
    return Polynomial<T>({0});
}

template <typename T>
Polynomial<T> Polynomial<T>::Power(size_t power) const {
    if (power == 0) {
        return Polynomial<T>::Identity();
    }
    Polynomial<T> result = *this;
    for (size_t i = 1; i < power; ++i) {
        result *= *this;
    }
    return result;
}

template <typename T>
const T Polynomial<T>::At(const T& point) const {
    T result = 0;
    for (const auto& [key, value] : coefficients_) {
        result += value * BinPower(point, key);
    }
    return result;
}

template <typename T>
const T Polynomial<T>::operator()(const T& point) const {
    return At(point);
}

template <typename T>
Polynomial<T> Polynomial<T>::Substitute(const Polynomial<T>& other) const {
    Polynomial<T> result;
    for (size_t i = 0; i < coefficients_.size(); ++i) {
        result += coefficients_.at(i) * other.Power(i);
    }
    return result;
}

template <typename T>
const std::map<size_t, T>& Polynomial<T>::GetCoefficients() const {
    return coefficients_;
}

template <typename T>
void Polynomial<T>::UnsignednessCheck() {
    if constexpr (std::is_unsigned_v<T>) {
        int status = 0;
        std::string T_type = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
        throw std::runtime_error(
            "An operation on Polynomial<" + (status == 0 ? T_type : typeid(T).name()) +
            "> of unsigned type requires one of the operands to become negative.");
    }
}

template <typename T>
T Polynomial<T>::BinPower(const T& value, size_t power) const {
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
