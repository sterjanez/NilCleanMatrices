# Nil-clean
An element a of an associative ring is called a nil-clean element if it is the sum of an idempotent e and a nilpotent q. The least index of nilpotence of q satisfying this property is called the nil-clean index of a.
This repository contains programs for finding nil-clean elements in the full matrix ring M_n(F_2), for certain k. The most
interesting are cases k=3 and k=4. See my math paper [Šter, On expressing matrices over Z_2 as the sum of an idempotent and a nilpotent, Linear Algebra and Its Applications, 2018] for more information.
## Nil-clean.pas
This Pascal code contains a series of commands for finding nil-clean matrices of index 3 and 4 in M_n(F_2). The program is in Slovene and the algorithms are not the fastest possible, but using it one can see that there exist nil-clean elements of index 4 if n=4 or n=8. For n=9 the program is too slow to solve the problem.
## Nil-clean.cpp
This C++ code contains a program to find, for a fixed matrix A over F_2, an idempotent matrix E such that (A-E)^3=0, or to prove that such an idempotent does not exist. The program uses an improved algorithm and is extremely fast. Using it one can show that in the matrix ring M_n(F_2), every element is nil-clean of index (at most) 3 if n=10 or n=11, and there exist elements which are nil-clean of index 4 if n=9. By recent paper of Shitov [Shitov, The ring M_{8k+4}(Z_2) is nil-clean of index four, Indagationes Mathematicae, 2019], the latter holds for n=12 as well. Hence for n<=12, M_n(F_2) is of index four exactly for n=4,8,9,12.
