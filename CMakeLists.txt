cmake_minimum_required(VERSION 3.15)
project(polynomial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)

add_library(polynomial
    src/SquareMatrix.cpp
    src/Polynomial.cpp
)

enable_testing()

include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
)
FetchContent_MakeAvailable(doctest)

add_executable(test_main tests/test_main.cpp)
target_link_libraries(test_main PRIVATE polynomial doctest::doctest)

add_test(NAME polynomial_tests COMMAND test_main)
