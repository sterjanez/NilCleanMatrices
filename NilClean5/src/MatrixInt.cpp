#include <cstdint>
#include "MatrixInt.h"
#include "Matrix.h"

// Multiples of 8.
static constexpr int multiples8[9] ={0, 8, 16, 24, 32, 40, 48, 56, 64};

// First row and column selector constants.
static constexpr std::uint64_t firstRow = 0xFF;
static constexpr std::uint64_t firstCol = 0x0101010101010101;

MatrixInt MI::create(Matrix const &A)
{
	MatrixInt result = 0;
	for (std::uint8_t i = 0; i < A.getHeight(); i++) {
		result ^= static_cast<std::uint64_t>(A.getRow(i)) << multiples8[i];
	}
	return result;
}

MatrixInt MI::create(std::uint8_t n)
{
	MatrixInt result = 0;
	for (std::uint8_t i = 0; i < n; i++) {
		result ^= static_cast<std::uint64_t>(1) << (multiples8[i] + i);
	}
	return result;
}

MatrixInt MI::create(std::uint8_t i, std::uint8_t j)
{
	return static_cast<std::uint64_t>(1) << (multiples8[i] + j);
}

MatrixInt MI::product(MatrixInt A, MatrixInt B)
{
	MatrixInt result = 0;
	for (std::uint8_t i = 0; i < 8; i++) {
		result ^= ((A >> i) & firstCol) * ((B >> multiples8[i]) & firstRow);
	}
	return result;
}

MatrixInt MI::rightMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j, MatrixInt B)
{
	return ((A >> i) & firstCol) * ((B >> multiples8[j]) & firstRow);
}

MatrixInt MI::leftMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j)
{
	return ((A >> multiples8[j]) & firstRow) << multiples8[i];
}

MatrixInt MI::rightMultiply(MatrixInt A, std::uint8_t i, std::uint8_t j)
{
	return ((A >> i) & firstCol) << j;
}

bool MI::isAX(MatrixInt C, MatrixInt A, std::uint8_t start_row)
{
	if (start_row) {
		A >>= multiples8[start_row];
		C >>= multiples8[start_row];
	}
	while (A) {
		std::uint64_t rows = (A & firstCol) * firstRow; // rows containing nonzero first elements
		std::uint64_t pivotRow = (rows & -static_cast<std::int64_t>(rows)) * firstRow; // lowest (pivot) row selector
		A ^= rows & ((pivotRow & A) * firstCol); // eliminate pivot row and reduce other nonzero rows
		C ^= rows & ((pivotRow & C) * firstCol); // corresponding operation in C
		A >>= 1;
	}
	return !C;
}

bool MI::isXA(MatrixInt C, MatrixInt A, std::uint8_t start_col)
{
	if (start_col) {
		A >>= start_col;
		C >>= start_col;
		std::uint64_t selector = firstCol *
			((static_cast<std::uint64_t>(1) << (8 - start_col)) - 1);
		A &= selector;
		C &= selector;
	}
	while (A) {
		std::uint64_t cols = (A & firstRow) * firstCol; // columns containing nonzero first elements
		std::uint64_t pivotCol = (cols & -static_cast<std::int64_t>(cols)) * firstCol; // lowest (pivot) column selector
		A ^= cols & ((pivotCol & A) * firstRow);
		C ^= cols & ((pivotCol & C) * firstRow);
		A >>= 8;
	}
	return !C;
}

std::uint8_t MI::rank(MatrixInt A, std::uint8_t m, std::uint8_t n, MatrixInt &U, MatrixInt &V)
{
	U = create(m);
	V = create(n);

	// We maintain U, V and the matrix UAV = [I00;0X0;000] with I of size rank x rank
	// and X of size (m - rank) x (maxrank - rank).
	std::uint8_t rank = 0;
	std::uint8_t maxrank = n;

	while (rank < maxrank) {
		std::uint64_t col = (A >> rank) & firstCol; // rank-th column of A
		if (col) {
			std::uint64_t rankEntry = static_cast<std::uint64_t>(1) << multiples8[rank];
			if ((col & rankEntry) == 0) {
				std::uint64_t pivotRow =
					(col & -static_cast<std::int64_t>(col)) * firstRow; // pivot row selector
				std::uint64_t rowA = pivotRow & A;
				std::uint64_t rowU = pivotRow & U;
				std::uint64_t pivotRankRow = firstRow << multiples8[rank];
				A ^= rowA == pivotRow ? pivotRankRow : pivotRankRow & (rowA / firstRow);
				U ^= rowU == pivotRow ? pivotRankRow : pivotRankRow & (rowU / firstRow); // move pivot to (rank,rank)
			}
			else {
				col ^= rankEntry;
			}
			A ^= ((A >> multiples8[rank]) & firstRow) * col;
			U ^= ((U >> multiples8[rank]) & firstRow) * col; // rank-th column of A is now (0,...,0,1,0,...,0)
			std::uint64_t row =
				((A >> multiples8[rank]) & firstRow) ^
				(static_cast<std::uint64_t>(1) << rank);
			A ^= row << multiples8[rank];
			V ^= ((V >> rank) & firstCol) * row;
			rank++;
		}
		else {
			maxrank--;
			if (rank != maxrank) {
				A ^= ((A >> maxrank) & firstCol) << rank;
				A ^= ((A >> rank) & firstCol) << maxrank;
				V ^= ((V >> maxrank) & firstCol) << rank;
				V ^= ((V >> rank) & firstCol) << maxrank;
			}
		}
	}
	return rank;
}

MatrixInt MI::inv(MatrixInt A, std::uint8_t n)
{
	MatrixInt U, V;
	std::uint8_t matrixRank = rank(A, n, n, U, V);
	return product(V, U);
}