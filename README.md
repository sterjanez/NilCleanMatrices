# Nil-Clean Matrices
An element a of an associative ring is called a nil-clean element if it is the sum of an idempotent e and a nilpotent q. The least index of nilpotence of q satisfying this property is called the nil-clean index of a.
This repository contains programs for finding nil-clean elements of nil-clean index k in the full matrix ring M_n(F_2), for certain k. The most interesting are cases k = 3 and k = 4. See my math paper [Šter, On expressing matrices over Z_2 as the sum of an idempotent and a nilpotent, Linear Algebra and Its Applications, 2018] for more information.
## Nil-clean1.pas
This Pascal code contains a series of commands for finding nil-clean matrices of index 3 and 4 in M_n(F_2). The program is in Slovene and the algorithms are not the fastest possible, but using it one can see that there exist nil-clean elements of index 4 if n = 4 or n = 8. For n = 9 the program is too slow to solve the problem. It takes aprroximately 10 hours to handle the case n = 8.
## Nil-clean2.cpp
This C++ code contains a program to find, for a fixed matrix A over F_2, an idempotent matrix E such that (A - E)^3 = 0, or to prove that such an idempotent does not exist. The program uses an improved algorithm and is extremely fast. Using it one can show that in the matrix ring M_n(F_2), there exist nil-clean elements of index 4 if n = 9, and that every element is nil-clean of index (at most) 3 if n = 10 or n = 11. (It takes about 10 minutes to handle the case n = 11, for a given input matrix.) By a recent paper of Shitov [Shitov, The ring M_{8k+4}(Z_2) is nil-clean of index four, Indagationes Mathematicae, 2019], there exist nil-clean elements of index 4 if n = 12. Hence for n <= 12, M_n(F_2) is of index four exactly when n = 4, 8, 9, 12.
## Nil-clean2.exe
Executable file for the C++ code in Nil-clean2.cpp.
## Nil-clean3.cpp
An improved version of Nil-clean2.cpp. Uses the same algorithm with bit-wise parallelization, where matrices of size 8×8 over F_2 are presented as single 64-bit integers. Uses some fast bit-wise computations, for example, a Gaussian elimination takes about 100 instructions (compared to 8^3 which is usual for this operation), and matrix multiplication is exectuded in about 50 operations (compared to 8^3 for the usual multiplication). Important: use 64-bit processor and compiler for optimal performance. The compiled program is about 10-times faster than the code in the Nil-clean2.cpp file. (The case n = 11 is handled in about 1 minute, for a fixed input matrix. It is estimated that handling the n = 13 case for a fixed matrix would take a couple of weeks.)

Usage: for an input n×n matrix A, write the list of last columns of companion matrices in the Frobenius canonical form of A. Flag "m" is reserved for modified companion blocks (see my LAA (2018) paper). Then enter the starting set, which may be any strictly increasing sequence of {0,...,n-1}. The set may also be empty (hit enter to enter the empty set).
## Nil-clean3.exe
Executable file for the file Nil-clean3.cpp.
## Nil-clean4.cpp
An extended version of Nil-clean3.cpp. Contains various additional methods: random heuristics, searching for matrices with common entries or common rows, finding irreducible polynomials over F_2 etc.
## Nil-clean4.exe
Executable file for the file Nil-clean4.cpp.
## NilClean5/
This version also includes parallelized algorithm (for multiple CPU). To build the program, either:
- run with PowerShell the file *build.ps1* (Windows system), or
- run with CMake the file *CMakeLists.txt* in the main directory.

Last update: 12th February 2022
