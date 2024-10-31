#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>
#include <utility>

#include "doctest/doctest.h"
#include "../include/SquareMatrix.hpp"
#include "../include/Polynomial.hpp"
#include "../include/MultivariatePolynomial.hpp"

#define ll long long

std::mt19937 rnd(time(nullptr));

std::pair<p::Polynomial<double>, std::vector<std::pair<double, size_t>>> GeneratePolynomial(
    size_t power_limit) {
    p::Variable<double> x;
    p::Polynomial<double> poly;
    std::vector<std::pair<double, size_t>> coeffs = {};
    for (size_t i = 0; i < 10; ++i) {
        double coeff = static_cast<double>(rnd() % 10000 - 5000) / 1000;
        size_t power = rnd() % power_limit;
        poly += coeff * (x ^ power);
        coeffs.push_back({coeff, power});
    }
    return {poly, coeffs};
}

p::SquareMatrix<ll, 3> GenerateMatrix() {
    std::vector<std::vector<ll>> matrix(3, std::vector<ll>(3));
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            matrix[i][j] = rnd() % 10 - 5;
        }
    }
    return p::SquareMatrix<ll, 3>(matrix);
}

std::pair<p::Polynomial<p::SquareMatrix<ll, 3>>,
          std::vector<std::pair<p::SquareMatrix<ll, 3>, size_t>>>
GenerateSMPolynomial(size_t power_limit) {
    p::Variable<p::SquareMatrix<ll, 3>> x;
    p::Polynomial<p::SquareMatrix<ll, 3>> poly;
    std::vector<std::pair<p::SquareMatrix<ll, 3>, size_t>> coeffs = {};
    for (size_t i = 0; i < 10; ++i) {
        p::SquareMatrix<ll, 3> coeff = GenerateMatrix();
        size_t power = rnd() % power_limit;
        poly += (coeff * (x ^ power));
        coeffs.push_back({coeff, power});
    }
    return {poly, coeffs};
}

template <typename T>
T BinPower(T value, size_t power) {
    if (power == 0) {
        return 1;
    }
    T res = BinPower(value, power / 2);
    if (power % 2 == 1) {
        return res * res * value;
    } else {
        return res * res;
    }
}

template <typename T>
T Evaluate(std::vector<std::pair<T, size_t>> coeffs, T at) {
    T result = 0;
    for (const auto &p : coeffs) {
        result += p.first * BinPower(at, p.second);
    }
    return result;
}

p::SquareMatrix<double, 3> EvaluateSMPolynomial(
    std::vector<std::pair<p::SquareMatrix<double, 3>, size_t>> coeffs,
    p::SquareMatrix<double, 3> at) {
    p::SquareMatrix<double, 3> result;
    for (const auto &p : coeffs) {
        result += p.first * BinPower(at, p.second);
    }
    return result;
}

void CheckCorrectness(double a, double b) {
    if (std::isinf(a) || std::isinf(b) || std::isnan(a) || std::isnan(b)) {
        WARN("Skipping test due to NAN or INF in evaluations.");
        return;
    }
    if (std::abs(std::max(a, b)) < 1e-10) {
        CHECK(std::abs(a - b) < 1e-10);
    } else {
        CHECK(std::min(a, b) / std::max(a, b) > (1 - 1e-10));
    }
}

void CheckMatrixCorrectness(p::SquareMatrix<ll, 3> a, p::SquareMatrix<ll, 3> b) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            CHECK_EQ(a[i][j], b[i][j]);
        }
    }
}

void RunTest(p::Polynomial<double> poly, std::vector<std::pair<double, size_t>> coeffs,
             size_t iters) {
    for (size_t dummy = 0; dummy < 100; ++dummy) {
        double value = static_cast<double>(rnd() % 10000 - 5000) / 1000;
        double poly_evaluation = poly(value);
        double actual_evaluation = Evaluate(coeffs, value);
        CheckCorrectness(poly_evaluation, actual_evaluation);
    }
}

void RunMatrixTest(p::Polynomial<p::SquareMatrix<ll, 3>> poly,
                   std::vector<std::pair<p::SquareMatrix<ll, 3>, size_t>> coeffs, size_t iters) {
    using Matrix = p::SquareMatrix<ll, 3>;
    for (size_t dummy2 = 0; dummy2 < iters; ++dummy2) {
        Matrix value = GenerateMatrix();
        Matrix poly_evaluation = poly(value);
        Matrix actual_evaluation = Evaluate(coeffs, value);
        CheckMatrixCorrectness(poly_evaluation, actual_evaluation);
    }
}

TEST_CASE("Warning") {
    WARN(
        "Different amount of test ran every time is caused by ignoring test where at least one of "
        "the operands is either undefined or infinite.");
}

TEST_CASE("SquareMatrix is implemented correctly") {
    p::SquareMatrix<ll, 2> m1({{1, 2}, {3, 4}});
    p::SquareMatrix<ll, 2> m1c = m1;
    p::SquareMatrix<ll, 2> m2({{5, 6}, {7, 8}});
    std::vector<std::vector<ll>> multiplication_result = {{19, 22}, {43, 50}};
    std::vector<std::vector<ll>> addition_result = {{6, 8}, {10, 12}};
    std::vector<std::vector<ll>> subtraction_result = {{-4, -4}, {-4, -4}};
    std::vector<std::vector<ll>> power_result = {{7, 10}, {15, 22}};

    CHECK_EQ((m1 * m2).GetMatrix(), multiplication_result);
    CHECK_EQ((m1 + m2).GetMatrix(), addition_result);
    CHECK_EQ((m1 - m2).GetMatrix(), subtraction_result);
    m1 += m2;
    CHECK_EQ(m1.GetMatrix(), addition_result);
    m1 -= m2;
    CHECK_EQ(m1.GetMatrix(), m1c.GetMatrix());
    m1 *= m1;
    CHECK_EQ(m1.GetMatrix(), power_result);
}

TEST_CASE("Polynomial class works with rather standard coefficients") {
    SUBCASE("Polynomial's static methods work correctly") {
        SUBCASE("Zero") {
            for (size_t dummy = 0; dummy < 1000; ++dummy) {
                auto poly_coeff = GeneratePolynomial(5);
                p::Polynomial<double> poly = poly_coeff.first + p::Polynomial<double>::Zero();
                std::vector<std::pair<double, size_t>> coeffs = poly_coeff.second;

                RunTest(poly, coeffs, 100);
            }
        }
        SUBCASE("Identity") {
            for (size_t dummy = 0; dummy < 1000; ++dummy) {
                auto poly_coeff = GeneratePolynomial(5);
                p::Polynomial<double> poly = poly_coeff.first * p::Polynomial<double>::Identity();
                std::vector<std::pair<double, size_t>> coeffs = poly_coeff.second;

                RunTest(poly, coeffs, 100);
            }
        }
        SUBCASE("Constant") {
            for (size_t dummy = 0; dummy < 1000; ++dummy) {
                auto poly_coeff = GeneratePolynomial(5);
                p::Polynomial<double> poly =
                    poly_coeff.first + p::Polynomial<double>::Constant(666);
                std::vector<std::pair<double, size_t>> coeffs = poly_coeff.second;
                for (size_t dummy2 = 0; dummy2 < 100; ++dummy2) {
                    double value = static_cast<double>(rnd() % 10000 - 5000) / 1000;
                    double poly_evaluation = poly(value);
                    double actual_evaluation = Evaluate(coeffs, value);
                    CheckCorrectness(poly_evaluation, actual_evaluation + 666);
                }
            }
        }
    }

    SUBCASE("Polynomial stress test") {
        for (size_t dummy = 0; dummy < 100; ++dummy) {
            auto poly_coeff = GeneratePolynomial(10);
            p::Polynomial<double> poly = poly_coeff.first;
            std::vector<std::pair<double, size_t>> coeffs = poly_coeff.second;

            RunTest(poly, coeffs, 100);
        }
    }
}

TEST_CASE("Polynomial class works with SquareMatrix class") {
    using Matrix = p::SquareMatrix<ll, 3>;
    SUBCASE("Polynomial's static methods work correctly") {
        SUBCASE("Zero") {
            for (size_t dummy = 0; dummy < 100; ++dummy) {
                auto poly_coeff = GenerateSMPolynomial(5);
                p::Polynomial<Matrix> poly = poly_coeff.first + p::Polynomial<Matrix>::Zero();
                std::vector<std::pair<Matrix, size_t>> coeffs = poly_coeff.second;

                RunMatrixTest(poly, coeffs, 10);
            }
        }
        SUBCASE("Identity") {
            for (size_t dummy = 0; dummy < 100; ++dummy) {
                auto poly_coeff = GenerateSMPolynomial(5);
                p::Polynomial<Matrix> poly = poly_coeff.first * p::Polynomial<Matrix>::Identity();
                std::vector<std::pair<Matrix, size_t>> coeffs = poly_coeff.second;

                RunMatrixTest(poly, coeffs, 10);
            }
        }
        SUBCASE("Constant") {
            for (size_t dummy = 0; dummy < 100; ++dummy) {
                auto poly_coeff = GenerateSMPolynomial(5);
                p::Polynomial<Matrix> poly =
                    poly_coeff.first + p::Polynomial<Matrix>::Constant(666);
                std::vector<std::pair<Matrix, size_t>> coeffs = poly_coeff.second;

                for (size_t dummy2 = 0; dummy2 < 10; ++dummy2) {
                    Matrix value = GenerateMatrix();
                    Matrix poly_evaluation = poly(value);
                    Matrix actual_evaluation = Evaluate(coeffs, value);
                    CheckMatrixCorrectness(poly_evaluation, actual_evaluation + 666);
                }
            }
        }
    }

    SUBCASE("Polynomial stress test") {
        for (size_t dummy = 0; dummy < 100; ++dummy) {
            auto poly_coeff = GenerateSMPolynomial(10);
            p::Polynomial<Matrix> poly = poly_coeff.first;
            std::vector<std::pair<Matrix, size_t>> coeffs = poly_coeff.second;

            RunMatrixTest(poly, coeffs, 100);
        }
    }
}

TEST_CASE("Important warning") {
    WARN(
        "Polynomial substitution and exponentiation are not tested as these two operations in the "
        "Polynomial class (apparently) only support commutative variables.");
}
