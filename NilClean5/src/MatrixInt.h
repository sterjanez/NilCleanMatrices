#ifndef MATRIX_INT_H
#define MATRIX_INT_H

#include <cstdint>
#include "Matrix.h"

class Matrix;

// Matrix entries, where each byte represents one row.
// E.g., the matrix unit E(1,1) is 00000001(bin) followed by seven zero bytes.
using MatrixInt = std::uint64_t;

namespace MI {

	// Convert m x n Matrix (where m,n <= 8) to Matrix8.
	MatrixInt create(Matrix const &A);

	// Diagonal matrix [I0;00] of rank n.
	MatrixInt create(std::uint8_t n = 0);

	// Matrix unit.
	MatrixInt create(std::uint8_t i, std::uint8_t j);
	
	// Product.
	MatrixInt product(MatrixInt A, MatrixInt B);

	// Right multiply by E_{ij} * A.
	MatrixInt rightMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j, MatrixInt B);

	// Left multiply with E_{ij}.
	MatrixInt leftMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j);

	// Right multiply with E_{ij}.
	MatrixInt rightMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j);

	// Operation involving LU decomposition:

	// Check if the matrix C equals AX for some X. Begin with start_row-th row.
	bool isAX(MatrixInt C, MatrixInt A, std::uint8_t start_row = 0);

	// Check if the matrix equals XA for some X. Begin with start_col-th column.
	bool isXA(MatrixInt C, MatrixInt A, std::uint8_t start_col = 0);

	// Find rank of a m x n matrix A and invertible
	// U (m x m) and V (n x n) such that UAV=[I0;00].
	std::uint8_t rank(MatrixInt A,
		std::uint8_t m, std::uint8_t n, MatrixInt &U, MatrixInt &V);

	// Inverse of n x n matrix.
	MatrixInt inv(MatrixInt A, std::uint8_t n);

}

#endif