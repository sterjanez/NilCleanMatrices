#ifndef MATRIX_H
#define MATRIX_H

#include <cstdint>
#include <string>
#include <vector>
#include "Matrix8.h"
#include "MatrixInt.h"

// Maximal matrix size, should not exceed 32.
static constexpr std::uint8_t MAX_DIM = 25;

class Matrix8;

using MatrixInt = std::uint64_t;

// Matrix of size m x n over F2, with m, n <= MAX_DIM.
class Matrix {

	std::uint8_t m_height;

	std::uint8_t m_width;

	// Matrix rows, where matrix entries are bits in each row.
	// For example, matrix [100;011] has rows
	// m_rows[0] = 001(bin) = 1 and m_rows[1] = 110(bin) = 6.
	// Other rows are set to zero.
	std::uint32_t m_rows[MAX_DIM];

public:

	std::uint8_t getHeight() const;

	std::uint8_t getWidth() const;

	bool getEntry(std::uint8_t i, std::uint8_t j) const;

	void setEntry(std::uint8_t i, std::uint8_t j, bool value);

	std::uint32_t getRow(std::uint8_t i) const;

	void setRow(std::uint8_t i, std::uint32_t row);

	std::uint32_t getCol(std::uint8_t j) const;

	void setCol(std::uint8_t j, std::uint32_t column);

	// Set matrix to zero or to identity (when m = n).
	Matrix& operator=(std::int8_t value);

	// Constructors:

	// Zero matrix.
	Matrix(std::uint8_t m = 0, std::uint8_t n = 0);

	// Convert string to matrix.
	// Direct summands are separated by a semicolon ';'.
	// If the string of a direct summand begins with "c",
	// create companion matrix with given last column.
	// If the string begins with "m", create a modified companion
	// matrix C + diag(1, 0, 1, 0, ...), where C is a companion matrix
	// with given last column and the diagonal matrix ends either with 1 or 0.
	// Otherwise the direct summand string defines rows separated by comma ','.
	// Example:
	// Matrix("10,01;c110;m110") -> [10;01] + [001;101;010] + [101;101;011]
	Matrix(std::string const &inputStr);

	// Convert Matrix8 to m x n Matrix.
	Matrix(std::uint8_t m, std::uint8_t n, Matrix8 A);

	// Convert MatrixInt to m x n Matrix.
	Matrix(std::uint8_t m, std::uint8_t n, MatrixInt A);

	// Assemble four Matrix8 blocks into a single m x n matrix,
	// with the upper left block of size i x j.
	Matrix(std::uint8_t m, std::uint8_t n, std::uint8_t i, std::uint8_t j,
		Matrix8 A11, Matrix8 A12, Matrix8 A21, Matrix8 A22);

	// Assemble four MatrixInt blocks into a single m x n matrix,
	// with the upper left block of size i x j.
	Matrix(std::uint8_t m, std::uint8_t n, std::uint8_t i, std::uint8_t j,
		MatrixInt A11, MatrixInt A12, MatrixInt A21, MatrixInt A22);

	// Conversion to string:

	// Convert matrix to string.
	// Example: toString([10;01]) -> "10,01".
	std::string toString() const;

	// Convert matrix to a user-friendly multi-line string.
	std::string toMultiLineString() const;

	// Relations:

	// Matrix equality. Matrices must be of equal sizes.
	bool operator==(Matrix const &A) const;

	// Check if the matrix equals zero or identity.
	bool operator==(std::int8_t value) const;

	// Matrix inequality. Matrices must be of equal sizes.
	bool operator!=(Matrix const &A) const;

	// Check if the matrix equals zero or identity.
	bool operator!=(std::int8_t value) const;

	// A well-ordering on the set of matrices of the same size.
	// Rows are ordered as binary numbers,
	// with the first row most significant.

	bool operator<(Matrix const &A) const;

	bool operator>(Matrix const &A) const;

	// Algebra operations:

	// Matrix addition.
	Matrix operator+(Matrix const &A) const;

	// Add zero or identity matrix.
	Matrix operator+(std::int8_t value) const;

	Matrix& operator+=(Matrix const &A);

	Matrix& operator+=(std::int8_t value);

	// Matrix multiplication.
	Matrix operator*(Matrix const &A) const;

	Matrix& operator*=(Matrix const &A);

	// Other functions:

	// Check if the product is zero.
	bool zeroProduct(Matrix const &A) const;

	Matrix transpose() const;

	// Trace of a square matrix.
	bool trace() const;

	// Count zero entries.
	std::uint16_t countZeros() const;

	// Matrix building and decomposing:

	// Submatrix with given upper left and lower right corner.
	Matrix submatrix(std::uint8_t i1, std::uint8_t j1,
		std::uint8_t i2, std::uint8_t j2) const;

	// Decompose into blocks, with the upper left block of size m x n.
	void decompose(std::uint8_t m, std::uint8_t n,
		Matrix &A11, Matrix &A12, Matrix &A21, Matrix &A22) const;

	// Decompose into Matrix8 blocks, with the upper left block of size m x n.
	void decompose(std::uint8_t m, std::uint8_t n,
		Matrix8 &A11, Matrix8 &A12, Matrix8 &A21, Matrix8 &A22) const;

	// Decompose into MatrixInt blocks, with the upper left block of size m x n.
	void decompose(std::uint8_t m, std::uint8_t n,
		MatrixInt &A11, MatrixInt &A12, MatrixInt &A21, MatrixInt &A22) const;

	// Paste a smaller matrix A into the given matrix
	// with the upper left corner position (i, j).
	void paste(Matrix const &A, std::uint8_t i, std::uint8_t j);

	// Block diagonal matrix with upper left block
	// the given matrix and lower right block A.
	Matrix blockDiagonal(Matrix const &A) const;

	// Row and column manipulation:

	void swapRows(std::uint8_t i1, std::uint8_t i2);

	void swapCols(std::uint8_t j1, std::uint8_t j2);

	// Add i1-th row to i2-th row.
	void addRow(std::uint8_t i1, std::uint8_t i2);

	// Add i1-th row of A to i2-th row.
	void addRow(std::uint8_t i1, Matrix const &A, std::uint8_t i2);

	// Add j1-th column to j2-th column.
	void addCol(std::uint8_t j1, std::uint8_t j2);

	// Add j1-th column of A to j2-th column.
	void addCol(std::uint8_t j1, Matrix const &A, std::uint8_t j2);

	// Add j-th column to columns determined by the variable columns.
	void addColToCols(std::uint8_t j, std::uint32_t columns);

	// LU decomposition procedures:

	// Check if the matrix can be written as AX for a given matrix A.
	// Verify this property for blocks starting with start_row-th row.
	bool isAX(Matrix const &A, std::uint8_t start_row = 0) const;

	// Check if the matrix can be written as XA for a given matrix A.
	// Verify this property for blocks starting with start_col-th column.
	bool isXA(Matrix const &A, std::uint8_t start_col = 0) const;

	// Find rank of the calling matrix A and
	// invertible square matrices U, V such that UAV = [I0;00].
	std::uint8_t rank(Matrix &U, Matrix &V) const;

	// Matrix inverse.
	Matrix inv() const;

	// Matrix power.
	Matrix operator^(std::int32_t power) const;

	// Other:

	// Return the list of all direct sums with summands companion matrices
	// of prescribed dimensions. For example, sizes can be "3,2,3,1".
	static std::vector<Matrix> companionMatrices(std::string const &sizes);

	// Given lists l1, ..., ln (with n >= 1) of matrices and a list of column indices,
	// find all matrices A1 in l1 such that there exist matrices A2, ..., An
	// in the other lists having the same columns with given indices as A1.
	static std::vector<Matrix> commonColumns(
		std::vector<std::vector<Matrix>> const &lists,
		std::vector<std::uint8_t> const &columns);

	// Given nonempty lists l1, ..., lk of m x n matrices (with k >= 1),
	// find an instance of indices i1, ..., ik
	// such that matrices l1[i1], ..., lk[ik] have maximal number of common entries.
	// For example, if l1 = {[01;11], [10;10]} and l2 = {[10,00]} then
	// i1 = 1 and i2 = 0 since [10;10] and [10;00] have three common entries
	// (which is maximal in this case).
	static std::vector<std::size_t> maxCommonEntries(
		std::vector<std::vector<Matrix>> const &lists, std::uint16_t &countCommon);

};

Matrix operator*(std::int8_t value, Matrix const &A);

#endif