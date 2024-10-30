#include "../include/Polynomial.hpp"

// namespace p {
// template <typename T>
// class Polynomial;
//
// template <typename T>
// Polynomial<T> Variable<T>::operator^(size_t power) {
//     Polynomial<T> result;
//     result.coefficients_[power] = 1;
//     return result;
// }
//
// template <typename T>
// Polynomial<T> operator*(T value, Variable<T> v) {
//     return value * (v ^ 1);
// }
//
// template <typename T>
// Polynomial<T> operator*(Variable<T> v, T value) {
//     return (v ^ 1) * value;
// }
//
// template <typename T>
// Polynomial<T>::Polynomial() = default;
//
// template <typename T>
// Polynomial<T>::Polynomial(const Polynomial<T>&) = default;
//
// template <typename T>
// Polynomial<T>::Polynomial(Polynomial<T>&&) = default;
//
// template <typename T>
// Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>&) = default;
//
// template <typename T>
// Polynomial<T>& Polynomial<T>::operator=(Polynomial<T>&&) = default;
//
// template <typename T>
// Polynomial<T>::Polynomial(std::vector<T> list) {
//     size_t i = 0;
//     for (const T& c : list) {
//         coefficients_[i] = c;
//         ++i;
//     }
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::operator-() {
//     UnsignednessCheck();
//     std::vector<T> coeff;
//     for (const auto& [key, value] : coefficients_) {
//         coeff.resize(key + 1);
//         coeff[key] = -value;
//     }
//     return Polynomial<T>(coeff);
// }
//
// template <typename T>
// Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& rhs) {
//     for (const auto& [power, c] : rhs.coefficients_) {
//         coefficients_[power] += c;
//     }
//     return *this;
// }
//
// template <typename T>
// Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& rhs) {
//     for (const auto& [power, c] : rhs.coefficients_) {
//         coefficients_[power] -= c;
//     }
//     return *this;
// }
//
// template <typename T>
// Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& rhs) {
//     std::map<size_t, T> new_coefficients;
//     for (const auto& m1 : coefficients_) {
//         for (const auto& m2 : rhs.coefficients_) {
//             new_coefficients[m1.first + m2.first] += m1.second * m2.second;
//         }
//     }
//     coefficients_ = new_coefficients;
//     return *this;
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::Identity() {
//     return Polynomial<T>({1});
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::Constant(const T& value) {
//     return Polynomial<T>({value});
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::Zero() {
//     return Polynomial<T>({});
// }
//
// template <typename T>
// Polynomial<T> operator+(Polynomial<T> lhs, Polynomial<T> rhs) {
//     lhs += rhs;
//     return lhs;
// }
//
// template <typename T>
// Polynomial<T> operator-(Polynomial<T> lhs, Polynomial<T> rhs) {
//     lhs -= rhs;
//     return lhs;
// }
//
// template <typename T>
// Polynomial<T> operator*(Polynomial<T> lhs, Polynomial<T> rhs) {
//     lhs *= rhs;
//     return lhs;
// }
//
// template <typename T>
// Polynomial<T> operator+(T value, Polynomial<T> rhs) {
//     rhs.coefficients_[0] += value;
//     return rhs;
// }
//
// template <typename T>
// Polynomial<T> operator+(Polynomial<T> lhs, T value) {
//     return value + lhs;
// }
//
// template <typename T>
// Polynomial<T> operator-(Polynomial<T> lhs, T value) {
//     if (lhs.coefficients_[0] < value) {
//         lhs.UnsignednessCheck();
//     }
//     lhs.coefficients_[0] -= value;
//     return lhs;
// }
//
// template <typename T>
// Polynomial<T> operator-(T value, Polynomial<T> rhs) {
//     return -(rhs - value);
// }
//
// template <typename T>
// Polynomial<T> operator*(T value, Polynomial<T> rhs) {
//     for (size_t i = 0; i < rhs.coefficients_.size(); ++i) {
//         rhs.coefficients_[i] = value * rhs.coefficients_[i];
//     }
//     return rhs;
// }
//
// template <typename T>
// Polynomial<T> operator*(Polynomial<T> lhs, T value) {
//     for (size_t i = 0; i < lhs.coefficients_.size(); ++i) {
//         lhs.coefficients_[i] *= value;
//     }
//     return lhs;
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::Power(size_t power) const {
//     if (power == 0) {
//         return Polynomial<T>::Identity();
//     }
//     Polynomial<T> result = *this;
//     for (size_t i = 1; i < power; ++i) {
//         result *= *this;
//     }
//     return result;
// }
//
// template <typename T>
// const T Polynomial<T>::At(const T& point) const {
//     T result = 0;
//     for (const auto& [key, value] : coefficients_) {
//         result += value * BinPower(point, key);
//     }
//     return result;
// }
//
// template <typename T>
// const T Polynomial<T>::operator()(const T& point) const {
//     return At(point);
// }
//
// template <typename T>
// Polynomial<T> Polynomial<T>::Substitute(const Polynomial<T>& other) const {
//     Polynomial<T> result;
//     for (size_t i = 0; i < coefficients_.size(); ++i) {
//         result += coefficients_.at(i) * other.Power(i);
//     }
//     return result;
// }
//
// template <typename T>
// const std::map<size_t, T>& Polynomial<T>::GetCoefficients() const {
//     return coefficients_;
// }
//
// template <typename T>
// void Polynomial<T>::UnsignednessCheck() {
//     if constexpr (std::is_unsigned_v<T>) {
//         int status = 0;
//         std::string T_type = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
//         throw std::runtime_error(
//             "An operation on Polynomial<" + (status == 0 ? T_type : typeid(T).name()) +
//             "> of unsigned type requires one of the operands to become negative.");
//     }
// }
//
// template <typename T>
// T Polynomial<T>::BinPower(const T& value, size_t power) const {
//     if (power == 0) {
//         return 1;
//     }
//     T result = BinPower(value, power / 2);
//     if (power % 2 == 1) {
//         return result * result * value;
//     } else {
//         return result * result;
//     }
// }
// }  // namespace p
