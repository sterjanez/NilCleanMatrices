#include <cstdint>
#include "Matrix8.h"
#include "Matrix.h"

std::uint64_t Matrix8::getEntries() const
{
	return m_mat;
}

void Matrix8::setEntries(std::uint64_t entries)
{
	m_mat = entries;
}

Matrix8& Matrix8::operator=(std::int8_t value)
{
	*this = Matrix8(value ? 8 : 0);
	return *this;
}

Matrix8::Matrix8(Matrix const &A)
{
	m_mat = 0;
	for (std::uint8_t i = 0; i < A.getHeight(); i++) {
		m_mat ^= static_cast<std::uint64_t>(A.getRow(i)) << multiples8[i];
	}
}

Matrix8::Matrix8(std::uint8_t n)
{
	m_mat = 0;
	for (std::uint8_t i = 0; i < n; i++) {
		m_mat ^= static_cast<std::uint64_t>(1) << (multiples8[i] + i);
	}
}

Matrix8::Matrix8(std::uint8_t i, std::uint8_t j)
{
	m_mat = static_cast<std::uint64_t>(1) << (multiples8[i] + j);
}

bool Matrix8::operator==(Matrix8 A) const
{
	return m_mat == A.m_mat;
}

bool Matrix8::operator!=(Matrix8 A) const
{
	return m_mat != A.m_mat;
}

Matrix8 Matrix8::operator+(Matrix8 A) const
{
	Matrix8 result;
	result.m_mat = m_mat ^ A.m_mat;
	return result;
}

Matrix8& Matrix8::operator+=(Matrix8 A)
{
	m_mat ^= A.m_mat;
	return *this;
}

Matrix8 Matrix8::operator*(Matrix8 A) const
{
	Matrix8 result;
	for (std::uint8_t i = 0; i < 8; i++) {
		result.m_mat ^= ((m_mat >> i) & firstCol) *
			((A.m_mat >> multiples8[i]) & firstRow);
	}
	return result;
}

Matrix8& Matrix8::operator*=(Matrix8 A)
{
	std::uint64_t newMat = 0;
	for (std::uint8_t i = 0; i < 8; i++) {
		newMat ^= ((m_mat >> i) & firstCol) *
			((A.m_mat >> multiples8[i]) & firstRow);
	}
	m_mat = newMat;
	return *this;
}

Matrix8 Matrix8::rightMultiply(std::uint8_t i, std::uint8_t j, Matrix8 A) const
{
	Matrix8 result;
	result.m_mat = ((m_mat >> i) & firstCol) *
		((A.m_mat >> multiples8[j]) & firstRow);
	return result;
}

Matrix8 Matrix8::leftMultiply(std::uint8_t i, std::uint8_t j) const
{
	Matrix8 result;
	result.m_mat = ((m_mat >> multiples8[j]) & firstRow) << multiples8[i];
	return result;
}

Matrix8 Matrix8::rightMultiply(std::uint8_t i, std::uint8_t j) const
{
	Matrix8 result;
	result.m_mat = ((m_mat >> i) & firstCol) << j;
	return result;
}

bool Matrix8::isAX(Matrix8 A, std::uint8_t start_row) const
{
	std::uint64_t matA = A.m_mat;
	std::uint64_t matB = m_mat;
	if (start_row) {
		matA = matA >> multiples8[start_row];
		matB = matB >> multiples8[start_row];
	}
	while (matA) {
		std::uint64_t rows = (matA & firstCol) * firstRow; // rows containing nonzero first elements
		std::uint64_t pivotRow = (rows & -static_cast<std::int64_t>(rows)) * firstRow; // lowest (pivot) row selector
		matA ^= rows & ((pivotRow & matA) * firstCol); // eliminate pivot row and reduce other nonzero rows
		matB ^= rows & ((pivotRow & matB) * firstCol); // corresponding operation in B
		matA >>= 1;
	}
	return !matB;
}

bool Matrix8::isXA(Matrix8 A, std::uint8_t start_col) const
{
	std::uint64_t matA = A.m_mat;
	std::uint64_t matB = m_mat;
	if (start_col) {
		matA >>= start_col;
		matB >>= start_col;
		std::uint64_t selector = firstCol *
			((static_cast<std::uint64_t>(1) << (8 - start_col)) - 1);
		matA &= selector;
		matB &= selector;
	}
	while (matA) {
		std::uint64_t cols = (matA & firstRow) * firstCol; // columns containing nonzero first elements
		std::uint64_t pivotCol = (cols & -static_cast<std::int64_t>(cols)) * firstCol; // lowest (pivot) column selector
		matA ^= cols & ((pivotCol & matA) * firstRow);
		matB ^= cols & ((pivotCol & matB) * firstRow);
		matA >>= 8;
	}
	return !matB;
}

std::uint8_t Matrix8::rank(std::uint8_t m, std::uint8_t n, Matrix8 &U, Matrix8 &V) const
{
	U = Matrix8(m);
	V = Matrix8(n);

	// We maintain U, V and the matrix UAV = [I00;0X0;000] with I of size rank x rank
	// and X of size (m - rank) x (maxrank - rank).
	std::uint8_t rank = 0;
	std::uint8_t maxrank = n;
	std::uint64_t matA = m_mat;

	while (rank < maxrank) {
		std::uint64_t col = (matA >> rank) & firstCol; // rank-th column of A
		if (col) {
			std::uint64_t rankEntry = static_cast<std::uint64_t>(1) << multiples8[rank];
			if ((col & rankEntry) == 0) {
				std::uint64_t pivotRow =
					(col & -static_cast<std::int64_t>(col)) * firstRow; // pivot row selector
				std::uint64_t rowA = pivotRow & matA;
				std::uint64_t rowU = pivotRow & U.m_mat;
				std::uint64_t pivotRankRow = firstRow << multiples8[rank];
				matA ^= rowA == pivotRow ? pivotRankRow : pivotRankRow & (rowA / firstRow);
				U.m_mat ^= rowU == pivotRow ? pivotRankRow : pivotRankRow & (rowU / firstRow); // move pivot to (rank,rank)
			}
			else {
				col ^= rankEntry;
			}
			matA ^= ((matA >> multiples8[rank]) & firstRow) * col;
			U.m_mat ^= ((U.m_mat >> multiples8[rank]) & firstRow) * col; // rank-th column of A is now (0,...,0,1,0,...,0)
			std::uint64_t row =
				((matA >> multiples8[rank]) & firstRow) ^
				(static_cast<std::uint64_t>(1) << rank);
			matA ^= row << multiples8[rank];
			V.m_mat ^= ((V.m_mat >> rank) & firstCol) * row;
			rank++;
		}
		else {
			maxrank--;
			if (rank != maxrank) {
				matA ^= ((matA >> maxrank) & firstCol) << rank;
				matA ^= ((matA >> rank) & firstCol) << maxrank;
				V.m_mat ^= ((V.m_mat >> maxrank) & firstCol) << rank;
				V.m_mat ^= ((V.m_mat >> rank) & firstCol) << maxrank;
			}
		}
	}
	return rank;
}

Matrix8 Matrix8::inv(std::uint8_t n)
{
	Matrix8 U, V;
	std::uint8_t matrixRank = rank(n, n, U, V);
	return V * U;
}