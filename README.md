# erlynum

A library for Erlang/OTP (and Elixir language) of numerical routines on   
single and double-precision real and complex number vectors and
matrices.

## Overview
The library provides a generic plain arrays support for Erlang. These
arrays might have a 'view' of vector (1-dimension) or matrix (2-dimension).

The most common operations to be supported {TODO}: summ, axpy, scale, cdot etc.

These operations mostly implemented using BLAST backend library (required at
runtime) to achieve possible utilization of SIMD instructions and gain the best performance.
 
The following BLAS backends are supported:
 * Open Source:
   - [ATLAS Library](http://math-atlas.sourceforge.net/) - one of the fastest 
   open source implementation
   - [Netlib BLAS Library](http://www.netlib.org/blas/) - mature library written in 
   Fortran with C-bindings
 * Non-Open Source:
   - [Intel Math Kernel Library](https://software.intel.com/en-us/intel-mkl/) - 
   extremely optimized for modern Intel processors. There is free version available 
   for non-commercial use.

## Status
The software is at early development stage now.
At present it allows to create vectors and matrices, and perform some BLAS Level 1
operations.

The future development (in order of priority):
  * Complete set of vector operations conforming to BLAS Level 1
  * Matrix to vector and vector to matrix operations (BLAS Level 2)
  * Matrix to matrix operations (BLAS Level 3)
  * Element-wise common math functions
  * Random values vector generations for several distributions
  * Common statistical functions (1-dimension and 2-dimension)
  * Something else. The final goal is to gain about NumPy functionality.

## Prerequirements

  * Erlang/OTP >= 18
  * `rebar3` (preferred) or `rebar` make tool
  * Modern GNU C Compliler (`gcc`) to build module itself
  * One of the cBLAS-interface libraries available at runtime 
  (not required to build)
         

## Build

```$ rebar3 compile edoc```
    
or

```$ rebar compile``` 
    

## API Reference

Build with edoc (supported by rebar3) and use it!