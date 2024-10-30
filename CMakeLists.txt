cmake_minimum_required(VERSION 3.15)
project(polynomial LANGUAGES CXX)

# Set C++20 standard explicitly
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include the directory for header files
include_directories(include)

# Library sources
add_library(polynomial
    src/SquareMatrix.cpp
    src/Polynomial.cpp
)

# Enable testing
enable_testing()

# Fetch doctest (https://github.com/doctest/doctest)
include(FetchContent)
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG v2.4.9 # You can change the version if needed
)
FetchContent_MakeAvailable(doctest)

# Add the test executable
add_executable(test_main tests/test_main.cpp)
target_link_libraries(test_main PRIVATE polynomial doctest::doctest)

# Add tests
add_test(NAME polynomial_tests COMMAND test_main)
