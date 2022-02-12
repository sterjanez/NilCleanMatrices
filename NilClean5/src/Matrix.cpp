#include <cstdint>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include "Matrix.h"
#include "Matrix8.h"
#include "MatrixInt.h"

std::uint8_t Matrix::getHeight() const
{
	return m_height;
}

std::uint8_t Matrix::getWidth() const
{
	return m_width;
}

bool Matrix::getEntry(std::uint8_t i, std::uint8_t j) const
{
	return m_rows[i] & (static_cast<std::uint32_t>(1) << j);
}

void Matrix::setEntry(std::uint8_t i, std::uint8_t j, bool value)
{
	if (value) {
		m_rows[i] |= static_cast<std::uint32_t>(1) << j;
	}
	else {
		m_rows[i] &= ~(static_cast<std::uint32_t>(1) << j);
	}
}

std::uint32_t Matrix::getRow(std::uint8_t i) const
{
	return m_rows[i];
}

void Matrix::setRow(std::uint8_t i, std::uint32_t row)
{
	m_rows[i] = row;
}

std::uint32_t Matrix::getCol(std::uint8_t j) const
{
	std::uint32_t result = 0;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] & (static_cast<std::uint32_t>(1) << j)) {
			result |= static_cast<std::uint32_t>(1) << i;
		}
	}
	return result;
}

void Matrix::setCol(std::uint8_t j, std::uint32_t column)
{
	std::uint32_t selector = static_cast<std::uint32_t>(1) << j;
	std::uint32_t negSelector = ~selector;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (column & (static_cast<std::uint32_t>(1) << i)) {
			m_rows[i] |= selector;
		}
		else {
			m_rows[i] &= negSelector;
		}
	}
}

Matrix& Matrix::operator=(std::int8_t value)
{
	if (value) {
		for (std::uint8_t i = 0; i < m_height; i++) {
			m_rows[i] = static_cast<std::uint32_t>(1) << i;
		}
	}
	else {
		for (std::uint8_t i = 0; i < m_height; i++) {
			m_rows[i] = 0;
		}
	}
	return *this;
}

Matrix::Matrix(std::uint8_t m, std::uint8_t n) :
	m_height{m},
	m_width{n}
{
	for (std::uint8_t i = 0; i < MAX_DIM; i++) {
		m_rows[i] = 0;
	}
}

Matrix::Matrix(std::string const &inputStr)
{
	std::size_t pos = inputStr.find(';');
	if (pos != std::string::npos) {
		*this = Matrix(inputStr.substr(0, pos)).blockDiagonal(inputStr.substr(pos + 1));
		return;
	}
	if (!inputStr.empty() && (inputStr[0] == 'c' || inputStr[0] == 'm')) {
		bool modified = inputStr[0] == 'm';
		std::string coefficients = inputStr.substr(1);
		std::uint8_t n = static_cast<std::uint8_t>(coefficients.length());
		*this = Matrix(n, n);
		for (std::uint8_t i = 1; i < n; i++) {
			setEntry(i, i - 1, true);
		}
		for (std::uint8_t i = 0; i < n; i++) {
			setEntry(i, n - 1, coefficients[i] == '1');
		}
		if (inputStr[0] == 'm') {
			for (std::uint8_t i = 0; 2 * i < n; i++) {
				setEntry(2 * i, 2 * i, !getEntry(2 * i, 2 * i));
			}
		}
		return;
	}
	if (inputStr.empty()) {
		*this = Matrix();
		return;
	}
	std::size_t posComma = inputStr.find(",");
	if (posComma == std::string::npos) {
		*this = Matrix(1, static_cast<std::uint8_t>(inputStr.length()));
	}
	else {
		*this = Matrix(
			static_cast<std::uint8_t>((inputStr.length() + 1) / (posComma + 1)),
			static_cast<std::uint8_t>(posComma));
	}
	std::uint8_t i = 0;
	std::uint8_t j = 0;
	for (char entry : inputStr) {
		if (entry == ',') {
			i++;
			j = 0;
		}
		else {
			setEntry(i, j, entry == '1');
			j++;
		}
	}
}

Matrix::Matrix(std::uint8_t m, std::uint8_t n, Matrix8 A)
{
	Matrix B(m, n);
	std::uint32_t selector = (static_cast<std::uint32_t>(1) << n) - 1;
	for (std::uint8_t i = 0; i < m; i++) {
		B.setRow(i, (A.getEntries() >> (i << 3)) & selector);
	}
	*this = B;
}

Matrix::Matrix(std::uint8_t m, std::uint8_t n, MatrixInt A)
{
	Matrix B(m, n);
	std::uint32_t selector = (static_cast<std::uint32_t>(1) << n) - 1;
	for (std::uint8_t i = 0; i < m; i++) {
		B.setRow(i, (A >> (i << 3)) & selector);
	}
	*this = B;
}

Matrix::Matrix(std::uint8_t m, std::uint8_t n, std::uint8_t i, std::uint8_t j,
	Matrix8 A11, Matrix8 A12, Matrix8 A21, Matrix8 A22)
{
	Matrix mA11(i, j, A11);
	Matrix mA12(i, n - j, A12);
	Matrix mA21(m - i, j, A21);
	Matrix mA22(m - i, n - j, A22);
	Matrix A(m, n);
	A.paste(mA11, 0, 0);
	A.paste(mA12, 0, j);
	A.paste(mA21, i, 0);
	A.paste(mA22, i, j);
	*this = A;
}

Matrix::Matrix(std::uint8_t m, std::uint8_t n, std::uint8_t i, std::uint8_t j,
	MatrixInt A11, MatrixInt A12, MatrixInt A21, MatrixInt A22)
{
	Matrix mA11(i, j, A11);
	Matrix mA12(i, n - j, A12);
	Matrix mA21(m - i, j, A21);
	Matrix mA22(m - i, n - j, A22);
	Matrix A(m, n);
	A.paste(mA11, 0, 0);
	A.paste(mA12, 0, j);
	A.paste(mA21, i, 0);
	A.paste(mA22, i, j);
	*this = A;
}


std::string Matrix::toString() const
{
	std::string result;
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			result += getEntry(i, j) ? "1" : "0";
		}
		if (i + 1 != m_height) {
			result += ",";
		}
	}
	return result;
}

std::string Matrix::toMultiLineString() const
{
	std::string emptyLine(static_cast<std::size_t>(m_width) * 2 + 1, ' ');
	std::string result;
	result += char(218);
	result += emptyLine;
	result += char(191);
	result += "\n";
	for (std::uint8_t i = 0; i < m_height; i++) {
		result += char(179);
		result += ' ';
		for (std::uint8_t j = 0; j < m_width; j++) {
			result += getEntry(i, j) ? "1 " : "0 ";
		}
		result += char(179);
		result += "\n";
	}
	result += char(192);
	result += emptyLine;
	result += char(217);
	return result;
}

bool Matrix::operator==(Matrix const &A) const
{
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] != A.m_rows[i]) {
			return false;
		}
	}
	return true;
}

bool Matrix::operator==(std::int8_t value) const
{
	if (value) {
		for (std::uint8_t i = 0; i < m_height; i++) {
			if (m_rows[i] != static_cast<std::uint32_t>(1) << i) {
				return false;
			}
		}
	}
	else {
		for (std::uint8_t i = 0; i < m_height; i++) {
			if (m_rows[i]) {
				return false;
			}
		}
	}
	return true;
}

bool Matrix::operator!=(Matrix const &A) const
{
	return !(*this == A);
}

bool Matrix::operator!=(std::int8_t value) const
{
	return !(*this == value);
}

bool Matrix::operator<(Matrix const &A) const
{
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] < A.m_rows[i]) {
			return true;
		}
		else if (m_rows[i] > A.m_rows[i]) {
			return false;
		}
	}
	return false;
}

bool Matrix::operator>(Matrix const &A) const
{
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] > A.m_rows[i]) {
			return true;
		}
		else if (m_rows[i] < A.m_rows[i]) {
			return false;
		}
	}
	return false;
}

Matrix Matrix::operator+(Matrix const &A) const
{
	Matrix B = *this;
	for (std::uint8_t i = 0; i < m_height; i++) {
		B.m_rows[i] ^= A.m_rows[i];
	}
	return B;
}

Matrix Matrix::operator+(std::int8_t value) const
{
	if (!value) {
		return *this;
	}
	Matrix A = *this;
	for (std::uint8_t i = 0; i < m_height; i++) {
		A.m_rows[i] ^= static_cast<std::uint32_t>(1) << i;
	}
	return A;
}

Matrix& Matrix::operator+=(Matrix const &A)
{
	*this = (*this) + A;
	return *this;
}

Matrix& Matrix::operator+=(std::int8_t value)
{
	*this = (*this) + value;
	return *this;
}

Matrix Matrix::operator*(Matrix const &A) const
{
	Matrix B(m_height, A.m_width);
	for (std::uint8_t j = 0; j < m_width; j++) {
		std::uint32_t selector = static_cast<std::uint32_t>(1) << j;
		for (std::uint8_t i = 0; i < m_height; i++) {
			if (m_rows[i] & selector) {
				B.m_rows[i] ^= A.m_rows[j];
			}
		}
	}
	return B;
}

Matrix& Matrix::operator*=(Matrix const &A)
{
	*this = (*this) * A;
	return *this;
}

bool Matrix::zeroProduct(Matrix const &A) const
{
	for (std::uint8_t i = 0; i < m_height; i++) {
		std::uint32_t sum = 0;
		for (std::uint8_t j = 0; j < m_width; j++) {
			if (m_rows[i] & (static_cast<std::uint32_t>(1) << j)) {
				sum ^= A.m_rows[j];
			}
		}
		if (sum) {
			return false;
		}
	}
	return true;
}

Matrix Matrix::transpose() const
{
	Matrix A(m_width, m_height);
	for (std::uint8_t i = 0; i < m_height; i++) {
		A.setCol(i, getRow(i));
	}
	return A;
}

bool Matrix::trace() const
{
	bool sum = false;
	for (std::uint8_t i = 0; i < m_height; i++) {
		sum ^= getEntry(i, i);
	}
	return sum;
}

std::uint16_t Matrix::countZeros() const
{
	std::uint16_t result = 0;
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			if (!getEntry(i, j)) {
				result++;
			}
		}
	}
	return result;
}

Matrix Matrix::submatrix(std::uint8_t i1, std::uint8_t j1,
	std::uint8_t i2, std::uint8_t j2) const
{
	Matrix A(i2 - i1, j2 - j1);
	for (std::uint8_t i = i1; i < i2; i++) {
		A.m_rows[i - i1] = (m_rows[i] &
			((static_cast<std::uint32_t>(1) << j2) - 1)) >> j1;
	}
	return A;
}

void Matrix::decompose(std::uint8_t m, std::uint8_t n,
	Matrix &A11, Matrix &A12, Matrix &A21, Matrix &A22) const
{
	A11 = submatrix(0, 0, m, n);
	A12 = submatrix(0, n, m, m_width);
	A21 = submatrix(m, 0, m_height, n);
	A22 = submatrix(m, n, m_height, m_width);
}

void Matrix::decompose(std::uint8_t m, std::uint8_t n,
	Matrix8 &A11, Matrix8 &A12, Matrix8 &A21, Matrix8 &A22) const
{
	Matrix mA11, mA12, mA21, mA22;
	decompose(m, n, mA11, mA12, mA21, mA22);
	A11 = Matrix8(mA11);
	A12 = Matrix8(mA12);
	A21 = Matrix8(mA21);
	A22 = Matrix8(mA22);
}

void Matrix::decompose(std::uint8_t m, std::uint8_t n,
	MatrixInt &A11, MatrixInt &A12, MatrixInt &A21, MatrixInt &A22) const
{
	Matrix mA11, mA12, mA21, mA22;
	decompose(m, n, mA11, mA12, mA21, mA22);
	A11 = Matrix8(mA11).getEntries();
	A12 = Matrix8(mA12).getEntries();
	A21 = Matrix8(mA21).getEntries();
	A22 = Matrix8(mA22).getEntries();
}

void Matrix::paste(Matrix const &A, std::uint8_t i, std::uint8_t j)
{
	for (std::uint8_t iRel = 0; iRel < A.m_height; iRel++) {
		m_rows[i + iRel] &= ~(((static_cast<std::uint32_t>(1) << A.m_width) - 1) << j);
		m_rows[i + iRel] |= A.m_rows[iRel] << j;
	}
}

Matrix Matrix::blockDiagonal(Matrix const &A) const
{
	Matrix B(m_height + A.m_height, m_width + A.m_width);
	B.paste(*this, 0, 0);
	B.paste(A, m_height, m_width);
	return B;
}

void Matrix::swapRows(std::uint8_t i1, std::uint8_t i2)
{
	std::swap(m_rows[i1], m_rows[i2]);
}

void Matrix::swapCols(std::uint8_t j1, std::uint8_t j2)
{
	std::uint32_t selector1 = static_cast<std::uint32_t>(1) << j1;
	std::uint32_t selector2 = static_cast<std::uint32_t>(1) << j2;
	std::uint32_t selector = selector1 ^ selector2;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] & selector1) {
			if (!(m_rows[i] & selector2)) {
				m_rows[i] ^= selector;
			}
		}
		else if (m_rows[i] & selector2) {
			m_rows[i] ^= selector;
		}
	}
}

void Matrix::addRow(std::uint8_t i1, std::uint8_t i2)
{
	m_rows[i2] ^= m_rows[i1];
}

void Matrix::addRow(std::uint8_t i1, Matrix const &A, std::uint8_t i2)
{
	m_rows[i2] ^= A.m_rows[i1];
}

void Matrix::addCol(std::uint8_t j1, std::uint8_t j2)
{
	std::uint32_t select1 = static_cast<std::uint32_t>(1) << j1;
	std::uint32_t select2 = static_cast<std::uint32_t>(1) << j2;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] & select1) {
			m_rows[i] ^= select2;
		}
	}
}

void Matrix::addCol(std::uint8_t j1, Matrix const &A, std::uint8_t j2)
{
	std::uint32_t select1 = static_cast<std::uint32_t>(1) << j1;
	std::uint32_t select2 = static_cast<std::uint32_t>(1) << j2;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (A.m_rows[i] & select1) {
			m_rows[i] ^= select2;
		}
	}
}

void Matrix::addColToCols(std::uint8_t j, std::uint32_t columns)
{
	std::uint32_t select = static_cast<std::uint32_t>(1) << j;
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (m_rows[i] & select) {
			m_rows[i] ^= columns;
		}
	}
}

bool Matrix::isAX(const Matrix& A, std::uint8_t start_row) const
{
	Matrix B = A;
	Matrix C = *this;
	for (std::uint8_t j0 = 0; j0 < B.m_width; j0++) {
		std::uint32_t selector = static_cast<std::uint32_t>(1) << j0;
		for (std::uint8_t i0 = start_row; i0 < m_height; i0++) {
			if (B.m_rows[i0] & selector) { // found pivot
				if (i0 != start_row) {
					B.swapRows(i0, start_row);
					C.swapRows(i0, start_row);
				}
				for (std::uint8_t i = i0 + 1; i < B.m_height; i++)
					if (B.m_rows[i] & selector) {
						B.addRow(start_row, i);
						C.addRow(start_row, i);
					}
				start_row++;
				break;
			}
		}
	}
	for (std::uint8_t i = start_row; i < m_height; i++) {
		if (C.m_rows[i]) {
			return false;
		}
	}
	return true;
}

bool Matrix::isXA(const Matrix& A, std::uint8_t start_col) const
{
	Matrix B = A;
	Matrix C = *this;
	for (std::uint8_t i0 = 0; i0 < B.m_height; i0++) {
		for (std::uint8_t j0 = start_col; j0 < m_width; j0++) {
			if (B.m_rows[i0] & (static_cast<std::uint32_t>(1) << j0)) { // found pivot
				if (j0 != start_col) {
					B.swapCols(j0, start_col);
					C.swapCols(j0, start_col);
				}
				std::uint32_t columns = B.m_rows[i0] &
					~(static_cast<std::uint32_t>(1) << start_col);
				B.addColToCols(start_col, columns);
				C.addColToCols(start_col, columns);
				start_col++;
				break;
			}
		}
	}
	std::uint32_t selector = (static_cast<std::uint32_t>(1) << m_width) -
		(static_cast<std::uint32_t>(1) << start_col);
	for (std::uint8_t i = 0; i < m_height; i++) {
		if (C.m_rows[i] & selector) {
			return false;
		}
	}
	return true;
}

std::uint8_t Matrix::rank(Matrix& U, Matrix& V) const
{
	Matrix B = *this;
	U = Matrix(m_height, m_height) + 1;
	V = Matrix(m_width, m_width) + 1;

	// We maintain equality UAV = B, where B is of the form
	// B = [I00;0X0] with I of size rank x rank
	// and X of size (m - rank) x (maxrank - rank).
	std::uint8_t rank = 0;
	std::uint8_t maxrank = m_width;
	std::uint32_t selector = static_cast<std::uint32_t>(1) << rank;

	while (rank < maxrank) {
		bool foundPivot = false;
		for (std::uint8_t i0 = rank; i0 < m_height; i0++) {
			if (B.m_rows[i0] & selector) {
				foundPivot = true;
				if (i0 != rank) {
					B.swapRows(i0, rank);
					U.swapRows(i0, rank);
				} // move pivot to position (rank, rank)
				for (std::uint8_t i = i0 + 1; i < m_height; i++) {
					if (B.m_rows[i] & selector) {
						B.addRow(rank, i);
						U.addRow(rank, i);
					}
				}
				V.addColToCols(rank, B.m_rows[rank] ^ selector);
				B.m_rows[rank] = selector;
				rank++;
				selector = static_cast<std::uint32_t>(1) << rank;
				break;
			}
		}
		if (!foundPivot) {
			maxrank--;
			if (rank != maxrank) {
				B.swapCols(rank, maxrank);
				V.swapCols(rank, maxrank);
			}
		}
	}
	return rank;
}

Matrix Matrix::inv() const
{
	Matrix U, V;
	std::uint8_t matrixRank = rank(U, V);
	return V * U;
}

Matrix Matrix::operator^(std::int32_t power) const
{
	if (power == 0) {
		return Matrix(m_height, m_height) + 1;
	}
	if (power < 0) {
		return inv() ^ (-power);
	}
	return ((*this) ^ (power - 1)) * (*this);
}

std::vector<Matrix> Matrix::companionMatrices(std::string const &sizes)
{
	std::size_t pos = sizes.find(',');
	if (pos != std::string::npos) {
		std::vector<Matrix> list1 = companionMatrices(sizes.substr(0, pos));
		std::vector<Matrix> list2 = companionMatrices(sizes.substr(pos + 1));
		std::vector<Matrix> list3;
		list3.reserve(list1.size() * list2.size());
		for (Matrix const &A : list1) {
			for (Matrix const &B : list2) {
				list3.push_back(A.blockDiagonal(B));
			}
		}
		return list3;
	}
	std::uint8_t n = std::stoi(sizes);
	std::vector<Matrix> list;
	list.reserve(static_cast<std::size_t>(1) << n);
	for (std::uint64_t counter = 0;
		counter < (static_cast<std::uint64_t>(1) << n); counter++)
	{
		std::string s = "c";
		for (std::uint8_t i = 0; i < n; i++) {
			s += (counter & (static_cast<std::uint64_t>(1) << i)) == 0 ? "0" : "1";
		}
		list.push_back(s);
	}
	return list;
}

std::vector<Matrix> Matrix::commonColumns(
	std::vector<std::vector<Matrix>> const &lists,
	std::vector<std::uint8_t> const &columns)
{
	std::vector<std::set<Matrix>> listSets;
	for (std::size_t k = 1; k < lists.size(); k++) {
		std::set<Matrix> listSet;
		for (Matrix const &A : lists[k]) {
			Matrix subA(A.getHeight(), static_cast<std::uint8_t>(columns.size()));
			for (std::uint8_t i = 0; i < subA.getHeight(); i++) {
				for (std::uint8_t j = 0; j < subA.getWidth(); j++) {
					subA.setEntry(i, j, A.getEntry(i, columns[j]));
				}
			}		
			listSet.insert(subA);
		}
		listSets.push_back(listSet);
	}
	std::vector<Matrix> solutions;
	for (Matrix const &A : lists[0]) {
		Matrix subA(A.getHeight(), static_cast<std::uint8_t>(columns.size()));
		for (std::uint8_t i = 0; i < subA.getHeight(); i++) {
			for (std::uint8_t j = 0; j < subA.getWidth(); j++) {
				subA.setEntry(i, j, A.getEntry(i, columns[j]));
			}
		}
		bool success = true;
		for (std::size_t i = 0; i < listSets.size() && success; i++) {
			if (listSets[i].find(subA) == listSets[i].end()) {
				success = false;
				break;
			}
		}
		if (success) {
			solutions.push_back(A);
		}
	}
	return solutions;
}

std::vector<std::size_t> Matrix::maxCommonEntries(
	std::vector<std::vector<Matrix>> const &lists, std::uint16_t &countCommon)
{
	std::uint8_t m = lists[0][0].getHeight();
	std::uint8_t n = lists[0][0].getWidth();
	std::uint16_t mn = static_cast<std::uint16_t>(m) * static_cast<std::uint16_t>(n);
	std::vector<std::size_t> maxIndices(lists.size(), 0);
	std::vector<std::size_t> indices(lists.size(), 0);
	countCommon = 0;
	bool finished = false;
	while (!finished) {
		Matrix commonOnes(m, n);
		Matrix commonZeros(m, n);
		for (std::uint8_t i = 0; i < n; i++) {
			std::uint32_t rowOnes = lists[0][indices[0]].m_rows[i];
			std::uint32_t rowZeros = rowOnes;
			for (std::size_t k = 1; k < lists.size(); k++) {
				rowOnes &= lists[k][indices[k]].m_rows[i];
				rowZeros |= lists[k][indices[k]].m_rows[i];
			}
			commonOnes.m_rows[i] = rowOnes;
			commonZeros.m_rows[i] = rowZeros;
		}
		std::uint16_t common = mn - commonOnes.countZeros() + commonZeros.countZeros();
		if (common > countCommon) {
			countCommon = common;
			maxIndices = indices;
		}
		finished = true;
		for (std::size_t k = 0; k < lists.size(); k++) {
			if (indices[k] == lists[k].size() - 1) {
				indices[k] = 0;
			}
			else {
				indices[k]++;
				finished = false;
				break;
			}
		}
	}
	return maxIndices;
}

Matrix operator*(std::int8_t value, Matrix const &A)
{
	return value ? A : Matrix(A.getHeight(), A.getWidth());
}