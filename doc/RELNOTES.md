Release Notes {#mainpage}
============

[TOC]

# Introduction {#Introduction}
vmmlib - a templatized C++ vector and matrix math library

VMMLib includes a vector and a matrix class, with additional
functionality for the often-used 3d and 4d vectors and 3x3 and 4x4
matrices. More advanced features include solvers, frustum computations
and frustum culling classes and spatial data structures

VMMLib is implemented using C++ templates, making it versatile. Being a
header library, it is very easy to integrate into other (your) libraries
and programs. There is no need to build and install a library, just
include the headers and you’re set. The BSD license allows the usage
both in open source and commercial closed source software.

# New in this release {#NewInThisRelease}

vmmlib 1.8 includes several fixes and enhancements over the last release, such as:

##Bug Fixes {#BugFixes}
* Sanitization of Matrix::get_translation API
* Compilation options and warnings cleanup
* Fix for shadowing member variables

##Enhancements {#Enhancements}
* Provide option to find project via find_file
* Added support for snapshot module
* Remove OSS package targets from CommonCPack to speed up CMake run and
  use less disk space for non-OSS packages
* Added TUVOK library
* Added various C++11 features

## Unit Tests {#UnitTests}
* Added test for C++11 template aliases

## Known Issues {#KnownIssues}
* Memory mapping for windows to be tested
* Test for slerp and matrix validators are not yet implemented
* Tests that depend on rand() are deactivated: they may break with different stdlib versions
* Tests that depend on LAPACK and BLAS are not fully supported for Windows

## Planned Future Extensions {#PlannedFutureExtensions}
* Decomposition and reconstruction algorithms for 4D tensors

# Documentation {#Documentation}

# About {#About}

# Platform Support {#PlatformSupport}
* Linux
* Mac OS X
* WIN64 operating systems

# Support {#Support}
* Comments and requests can be issued at http://github.com/VMML/vmmlib/issues
* Contributions can be merged into vmmlib via a pull request

# Errata {#Errata}
