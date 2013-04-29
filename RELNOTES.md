# Release Notes
2013/04/29
vmmlib - a templatized C++ vector and matrix math library

Its basic functionality includes a vector and a matrix class, with additional functionality for the often-used 3d and 4d vectors and 3x3 and 4x4 matrices.
More advanced features include solvers, frustum computations and frustum culling classes, and spatial data structures. vmmlib also offers support for manipulating 3rd-order and 4th-order tensors, as well as several algorithms for tensor approximation.

Vmmlib is implemented using C++ templates, making it versatile. Being a header library, it is very easy to integrate into other (your) libraries and programs. There is no need to build and install a library, just include the headers and you’re set.
The BSD license allows the usage both in open source and commercial closed source software.

# Features

# New in this release


###CMake
* use now only CMake to create a Makefile

###Use Tensor Approximation Classes (based on LAPACK, BLAS)
* define VMMLIB_USE_LAPACK in Makefile (changed from VMMLIB_BASIC_ONLY)

### Matrix
* added matrix validators

##Tensor Approximation Classes
* memory mapping classes for tensor data structures (tensor_mmapper.hpp) 
* created t3_converter: used to import/export tensors; converting/reading/writing files; quantization
* moved tensor3 reading/writing/quantization/tensor times matrix multiplications to t3_converter.hpp
* added tensor times matrix multiplications (ttm) and tensor times vector multiplication (ttv): multiply tensors along various modes
* 
* use OpenMP for some parallel computeations (tensor3, ttm)
* added tensor4 data structure (including t4_converter for converting/reading/writing files)
* added tensor_stats TODO RAFA
* fixed bug in CP3 tensor decomposition 
* updated incremental methods (added ihooi; updated ihopm)

###BLAS WRAPPERS
* added BLAS DAXPY wrapper



## New Features

## Enhancements

## Optimizations

## Tools

## Documentation

## Bug Fixes

## Known Bugs

* SVD LAPACK wrapper: TODO RAFA
* Memory Mapping for windows

## Planned Future Extensions

* tensor 4 decomposition and reconstruction

# About

## Platform Support

## Documentation

## Support

* Linux
* Mac OS X
* only supports WIN64 operating systems

# Errata
