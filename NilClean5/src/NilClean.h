#ifndef NIL_CLEAN_H
#define NIL_CLEAN_H

#include <cstdint>
#include <vector>
#include "Matrix.h"
#include "IntegerSet.h"

// Methods for finding nil-clean decompositions of index 3 of matrices over F_2.
namespace NilClean {

// Given n x n matrix A and 0 < p < n, find idempotents E of rank p in block form
// E = [I0;XI][IY;00][I0;XI] (with the upper left block of size p x p),
// such that (A - E) ^ 3 = 0. Heights of columns in block X are limited by heights[].
// Example: heights[] = {1, 1, ..., 1} means that X can only have the top nonzero row.
// If find_one is set, find only one idempotent.
// Append idempotents to the vector idempotents.
// Use Matrix8 matrices. Only for use if p, n - p <= 8.
// If rand > 0 then use random heuristic with rand attempts.
void findIdempotents8(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand);

// Given n x n matrix A and 0 < p < n, find idempotents E of rank p in block form
// E = [I0;XI][IY;00][I0;XI] (with the upper left block of size p x p),
// such that (A - E) ^ 3 = 0. Heights of columns in block X are limited by heights[].
// Example: heights[] = {1, 1, ..., 1} means that X can only have the top nonzero row.
// If find_one is set, find only one idempotent.
// Append idempotents to the vector idempotents.
// Use MatrixInt matrices. Only for use if p, n - p <= 8.
// If rand > 0 then use random heuristic with rand attempts.
void findIdempotentsInt(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand);

// Identical function as findIdempotents8, except that it uses Matrix class
// and works for arbitrary n <= MAX_DIM.
void findIdempotents(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand);

// Given n x n matrix A and a subset M of {0, ..., n - 1}, find idempotents E
// such that (A - E) ^ 3 = 0 and such that P^T * E * P is of the form [I0;XI][IY;00][I0;XI],
// where P is a permuation matrix associated with M. If find_one is set,
// find only one idempotent. Append idempotents to the vector idempotents.
// If rand > 0, use random heuristic with rand attempts.
void findIdempotents(Matrix const &A, IntegerSet const &M,
	bool find_one, std::vector<Matrix> &idempotents, std::uint64_t rand);

// Given n x n matrix A, find idempotents E such that (A - E) ^ 3 = 0.
// If show_steps is set, write steps in the procedure. Begin with the subset M0.
// If find_one is set, then find just one E.
// Append idempotents to the vector idempotents.
// If rand > 0, use random heuristic with rand attempts.
// If rand_sets > 0 then rand > 0 and use random heuristic with rand_sets sets.
void findIdempotents(Matrix const &A, bool find_one,
	IntegerSet const &M0, std::vector<Matrix> &idempotents,
	std::uint64_t rand, std::uint64_t rand_sets, bool show_steps = true);

// Same as findIdempotents. Use multi-threading. No option "random", no starting subset.
void findIdempotentsParallel(Matrix const &A, bool find_one,
	std::vector<Matrix> &idempotents, bool show_steps = true);

}

#endif