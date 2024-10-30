# polynomial_library

Small library for working with polynomials in C++

## Dependencies

- `cmake`;
- C++20-ready compiler;
- `doctest` testing framework will be downloaded automatically.

## Usage pipeline

### Initialization and testing

From the root directory:

```sh
mkdir build/
cd build/
cmake ..
make test_main
./test_main
```

### Actual usage

Include headers from `include/`, invoke stored classes and functions with `p::` namespace prefix. Example `p::Polynomial` usage:

```cpp
    ...
    p::Variable<int> x;
    p::Polynomial<int> poly = (x ^ 2) - 1;
    std::cout << poly(3);  // prints out `8`
    ...
```

## P.S.

Everything was tested on latest updated Arch Linux distribution.
