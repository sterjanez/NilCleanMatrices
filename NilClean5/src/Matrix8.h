#ifndef MATRIX8_H
#define MATRIX8_H

#include <cstdint>
#include "Matrix.h"

class Matrix;

// Matrix over F_2 of size 8 x 8 as a single 64-bit integer.
class Matrix8 {

	// Multiples of 8.
	static constexpr int multiples8[9] = {0, 8, 16, 24, 32, 40, 48, 56, 64};

	// First row and column selector constants.
	static constexpr std::uint64_t firstRow = 0xFF;
	static constexpr std::uint64_t firstCol = 0x0101010101010101;

	// Matrix entries, where each byte represents one row.
	// E.g., the matrix unit E(1,1) is 00000001(bin) followed by seven zero bytes.
	std::uint64_t m_mat;

public:

	std::uint64_t getEntries() const;

	void setEntries(std::uint64_t entries);

	// Set matrix to zero or 8 x 8 identity.
	Matrix8& operator=(std::int8_t value);

	// Constructors:

	// Convert m x n Matrix (where m,n <= 8) to Matrix8.
	Matrix8(Matrix const &A);

	// Diagonal matrix [I0;00] of rank n.
	Matrix8(std::uint8_t n = 0);

	// Matrix unit.
	Matrix8(std::uint8_t i, std::uint8_t j);

	// Relations:

	// Matrix equality.
	bool operator==(Matrix8 A) const;

	// Matrix inequality.
	bool operator!=(Matrix8 A) const;

	// Algebra operations:

	Matrix8 operator+(Matrix8 A) const;

	Matrix8& operator+=(Matrix8 A);

	Matrix8 operator*(Matrix8 A) const;

	// Matrix product.
	Matrix8& operator*=(Matrix8 A);

	// Right multiply by E_{ij} * A.
	Matrix8 rightMultiply(std::uint8_t i, std::uint8_t j, Matrix8 A) const;

	// Left multiply with E_{ij}.
	Matrix8 leftMultiply(std::uint8_t i, std::uint8_t j) const;

	// Right multiply with E_{ij}.
	Matrix8 rightMultiply(std::uint8_t i, std::uint8_t j) const;

	// Operation involving LU decomposition:

	// Check if the matrix equals AX for some X. Begin with start_row-th row.
	bool isAX(Matrix8 A, std::uint8_t start_row = 0) const;

	// Check if the matrix equals XA for some X. Begin with start_col-th column.
	bool isXA(Matrix8 A, std::uint8_t start_col = 0) const;

	// Find rank of a m x n matrix A and invertible
	// U (m x m) and V (n x n) such that UAV=[I0;00].
	std::uint8_t rank(std::uint8_t m, std::uint8_t n, Matrix8 &U, Matrix8 &V) const;

	// Inverse of n x n matrix.
	Matrix8 inv(std::uint8_t n);

};

#endif