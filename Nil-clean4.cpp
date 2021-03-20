// Program: Nil-clean matrices over F_2
// Author: Janez Ster
// Last update: February 2021

// Description: Given a n×n matrix A over F_2 (for small n), the program finds
// idempotent matrices E over F_2 such that (A-E)^3=0.

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>
#include <random>
#include <set>

const int MAX_DIM = 25; // maximal matrix size, must not exceed 32

// Class describing polynomial over F_2.
class Polynomial {

public:

	// Polynomial coefficients. If the polynomial is zero, the coefficient vector is empty.
	// Otherwise the leading term of Coef is 1.
	std::vector<bool> Coef;

	// Create zero polynomial.
	Polynomial() {}

	// Monomial of degree n.
	Polynomial(int n) {
		for (int i = 0; i < n; i++)
			Coef.push_back(false);
		Coef.push_back(true);
	}

	// Convert string to polynomial. Example: "01101" is converted to x+x^2+x^4.
	Polynomial(const std::string& s) {
		for (char c : s)
			Coef.push_back(c == '1');
	}

	// Polynomial degree. Degree of 0 is 0.
	int Degree() const {
		if (Coef.empty())
			return 0;
		return Coef.size() - 1;
	}

	// Convert polynomial to string. Example: 1+x^2 is converted to "1 + x^2".
	std::string ToString() const {
		if (Coef.empty())
			return "0";
		std::string s;
		for (int i = 0; i < Coef.size(); i++)
			if (Coef.at(i)) {
				if (s != "")
					s += " + ";
				if (i == 0)
					s += "1";
				else if (i == 1)
					s += "x";
				else {
					s += "x^";
					s += std::to_string(i);
				}
			}
		return s;
	}

	// Polynomial sum.
	Polynomial operator +(const Polynomial& p) const {
		int maxdegree = Coef.size();
		if (p.Coef.size() > maxdegree)
			maxdegree = p.Coef.size();
		Polynomial sum;
		for (int i = 0; i < maxdegree; i++) {
			bool value = false;
			if (i < Coef.size())
				value ^= Coef.at(i);
			if (i < p.Coef.size())
				value ^= p.Coef.at(i);
			sum.Coef.push_back(value);
		}
		while (!sum.Coef.empty() && !sum.Coef.back())
			sum.Coef.pop_back();
		return sum;
	}

	// Polynomial product.
	Polynomial operator *(const Polynomial& p) const {
		Polynomial sum;
		for (int i = 0; i < Coef.size(); i++)
			for (int j = 0; j < p.Coef.size(); j++)
				if (Coef.at(i) && p.Coef.at(j))
					sum = sum + Polynomial(i + j);
		return sum;
	}

	// Polynomial equality.
	bool operator ==(const Polynomial& p) const {
		if (Coef.size() != p.Coef.size())
			return false;
		for (int i = 0; i < Coef.size(); i++)
			if (Coef.at(i) != p.Coef.at(i))
				return false;
		return true;
	}

	// Polynomial inequality.
	bool operator !=(const Polynomial& p) const {
		return !(*this == p);
	}

	// Check if the polynomial is equal to a constant polynomial.
	bool operator ==(int k) const {
		return Coef.size() == k;
	}

	// Check if the polynomial is not equal to a constant polynomial.
	bool operator !=(int k) const {
		return !(*this == k);
	}

	// Polynomial quotient.
	Polynomial operator /(const Polynomial& p) const {
		Polynomial quotient;
		Polynomial remainder = *this;
		while (remainder.Coef.size() >= p.Coef.size()) {
			int n = remainder.Coef.size() - p.Coef.size();
			quotient = quotient + Polynomial(n);
			remainder = remainder + Polynomial(n) * p;
		}
		return quotient;
	}

	// Division remainder.
	Polynomial Mod(const Polynomial& p) const {
		return *this + (*this / p) * p;
	}

	// Replace a polynomial with its successor in the well-ordered set of all polynomials.
	// Polynomials are ordered as binary numbers, for example,
	// 0 < 1 < x^3 < 1 + x^3 < x + x^3 < x^4.
	void Inc() {
		for (int i = 0; i < Coef.size(); i++)
			if (Coef.at(i))
				Coef.at(i) = false;
			else {
				Coef.at(i) = true;
				return;
			}
		Coef.push_back(true);
	}

	// Check if the polynomial is irreducible. Slow algorithm.
	bool Irreducible() const {
		if (Degree() == 0)
			return false;
		Polynomial divisor_candidate = Polynomial(1);
		while (2 * divisor_candidate.Degree() <= Degree()) {
			if (Mod(divisor_candidate) == 0)
				return false;
			divisor_candidate.Inc();
		}
		return true;
	}

};

// Class describing m×n matrix over F_2.
class Matrix {

public:

	// Matrix dimensions.
	int Height, Width;

	// Matrix rows. Matrix entries are bits in each row.
	// For example, matrix [100;011] has rows
	// Row[0] = 001(bin) = 1 and Row[1] = 110(bin) = 6.
	int Row[MAX_DIM];

	// Create 0×0 matrix.
	Matrix() {
		Height = 0;
		Width = 0;
	}

	// Create zero m×n matrix.
	Matrix(int m, int n) {
		Height = m;
		Width = n;
		for (int i = 0; i < Height; i++)
			Row[i] = 0;
	}

	// Get matrix entry.
	inline bool GetEntry(int i, int j) const {
		return Row[i] & (1 << j);
	}

	// Set matrix entry.
	inline void SetEntry(int i, int j, bool value) {
		if (value)
			Row[i] |= 1 << j;
		else
			Row[i] &= ~(1 << j);
	}

	// Set a m×n matrix either to zero or a n×n matrix to identity.
	Matrix& operator =(int value) {
		if (value)
			for (int i = 0; i < Height; i++)
				Row[i] = 1 << i;
		else
			for (int i = 0; i < Height; i++)
				Row[i] = 0;
		return *this;
	}

	// Create matrix from a string. Example: "10,01" presents 2×2 identity.
	Matrix(const std::string& s) {
		std::size_t pos = s.find(",");
		if (pos == std::string::npos) {
			if (s == "") {
				Height = 0;
				Width = 0;
			}
			else {
				Height = 1;
				Width = s.length();
			}
		}
		else {
			Height = (s.length() + 1) / (pos + 1);
			Width = pos;
		}
		for (int i = 0; i < Height; i++) {
			Row[i] = 0;
			for (int j = 0; j < Width; j++)
				SetEntry(i, j, s[i * (Width + 1) + j] == '1');
		}
	}

	// Convert matrix to string. Example: 2×2 identity is converted to "10,01".
	std::string ToString() const {
		std::string s;
		for (int i = 0; i < Height; i++) {
			for (int j = 0; j < Width; j++)
				s += GetEntry(i, j) ? "1" : "0";
			if (i != Height - 1)
				s += ",";
		}
		return s;
	}

	// Convert matrix to a user-friendly multi-line string.
	std::string ToMultiLineString() const {
		std::string s;
		s += (char)218;
		for (int i = 0; i < 2 * Width + 1; i++)
			s += " ";
		s += (char)191;
		s += "\n";
		for (int i = 0; i < Height; i++) {
			s += (char)179;
			s += " ";
			for (int j = 0; j < Width; j++)
				s += GetEntry(i, j) ? "1 " : "0 ";
			s += (char)179;
			s += "\n";
		}
		s += (char)192;
		for (int i = 0; i < 2 * Width + 1; i++)
			s += " ";
		s += (char)217;
		return s;
	}

	// Matrix sum.
	Matrix operator +(const Matrix& A) const {
		Matrix B = *this;
		for (int i = 0; i < Height; i++)
			B.Row[i] ^= A.Row[i];
		return B;
	}

	// Add a scalar matrix to a n×n matrix.
	Matrix operator +(int value) const {
		if (!value)
			return *this;
		Matrix A = *this;
		for (int i = 0; i < Height; i++)
			A.Row[i] ^= 1 << i;
		return A;
	}

	// Matrix product.
	Matrix operator *(const Matrix& A) const {
		Matrix B(Height, A.Width);
		for (int j = 0; j < Width; j++) {
			int selector = 1 << j;
			for (int i = 0; i < Height; i++)
				if (Row[i] & selector)
					B.Row[i] ^= A.Row[j];
		}
		return B;
	}

	// Matrix equality. Matrices must be of equal sizes.
	inline bool operator ==(const Matrix& A) const {
		for (int i = 0; i < Height; i++)
			if (Row[i] != A.Row[i]) return false;
		return true;
	}

	// Matrix inequality. Matrices must be of equal sizes.
	inline bool operator !=(const Matrix& A) const {
		return !(*this == A);
	}

	// Check if a n×n matrix equals a scalar matrix.
	bool operator ==(int value) const {
		if (value) {
			for (int i = 0; i < Height; i++)
				if (Row[i] != 1 << i) return false;
		}
		else {
			for (int i = 0; i < Height; i++)
				if (Row[i]) return false;
		}
		return true;
	}

	// Check if a n×n matrix equals a scalar matrix.
	bool operator !=(int value) const {
		return !(*this == value);
	}

	// Check if the product is zero.
	bool ZeroProduct(const Matrix& A) const {
		for (int i = 0; i < Height; i++) {
			int sum = 0;
			for (int j = 0; j < Width; j++) {
				if (Row[i] & (1 << j))
					sum ^= A.Row[j];
			}
			if (sum) return false;
		}
		return true;
	}

	// Swap rows.
	inline void RowSwap(int i1, int i2) {
		int row = Row[i1];
		Row[i1] = Row[i2];
		Row[i2] = row;
	}

	// Swap columns.
	inline void ColSwap(int j1, int j2) {
		int selector1 = 1 << j1;
		int selector2 = 1 << j2;
		int selector = selector1 ^ selector2;
		for (int i = 0; i < Height; i++)
			if (Row[i] & selector1) {
				if (!(Row[i] & selector2))
					Row[i] ^= selector;
			}
			else if (Row[i] & selector2)
				Row[i] ^= selector;
	}

	// Add i1-th row to i2-th row.
	inline void AddRow(int i1, int i2) {
		Row[i2] ^= Row[i1];
	}

	// Add i1-th row of A to i2-th row.
	inline void AddRow(int i1, const Matrix& A, int i2) {
		Row[i2] ^= A.Row[i1];
	}

	// Add j1-th column to j2-th column.
	inline void AddCol(int j1, int j2) {
		int select1 = 1 << j1;
		int select2 = 1 << j2;
		for (int i = 0; i < Height; i++)
			if (Row[i] & select1) Row[i] ^= select2;
	}

	// Add j1-th column of A to j2-th column.
	inline void AddCol(int j1, const Matrix& A, int j2) {
		int select1 = 1 << j1;
		int select2 = 1 << j2;
		for (int i = 0; i < Height; i++)
			if (A.Row[i] & select1) Row[i] ^= select2;
	}

	// Add j-th column to other (different) columns.
	inline void AddColToCols(int j, int columns) {
		int selector = 1 << j;
		for (int i = 0; i < Height; i++)
			if (Row[i] & selector)
				Row[i] ^= columns;
	}

	// Check if the matrix can be written as AX for a given matrix A.
	// Verify this property for blocks starting with start_row-th row.
	bool IsAX(const Matrix& A, int start_row = 0) const {
		Matrix B = A;
		Matrix C = *this;
		for (int j0 = 0; j0 < B.Width; j0++) {
			int selector = 1 << j0;
			for (int i0 = start_row; i0 < Height; i0++)
				if (B.Row[i0] & selector) { // found pivot
					if (i0 != start_row) {
						B.RowSwap(i0, start_row);
						C.RowSwap(i0, start_row);
					}
					for (int i = i0 + 1; i < B.Height; i++)
						if (B.Row[i] & selector) {
							B.AddRow(start_row, i);
							C.AddRow(start_row, i);
						}
					start_row++;
					break;
				}
		}
		for (int i = start_row; i < Height; i++)
			if (C.Row[i]) return false;
		return true;
	}

	// Check if the matrix can be written as XA for a given matrix A.
	// Verify this property for blocks starting with start_col-th column.
	bool IsXA(const Matrix& A, int start_col = 0) const {
		Matrix B = A;
		Matrix C = *this;
		for (int i0 = 0; i0 < B.Height; i0++) {
			for (int j0 = start_col; j0 < Width; j0++)
				if (B.Row[i0] & (1 << j0)) { // found pivot
					if (j0 != start_col) {
						B.ColSwap(j0, start_col);
						C.ColSwap(j0, start_col);
					}
					int columns = B.Row[i0] & ~(1 << start_col);
					B.AddColToCols(start_col, columns);
					C.AddColToCols(start_col, columns);
					start_col++;
					break;
				}
		}
		int selector = (1 << Width) - (1 << start_col);
		for (int i = 0; i < Height; i++)
			if (C.Row[i] & selector) return false;
		return true;
	}

	// Find rank of the calling matrix A and
	// invertible square matrices U, V such that UAV = [I0;00].
	void Diagonalize(int& rank, Matrix& U, Matrix& V) const {
		Matrix B = *this;
		U = Matrix(Height, Height) + 1;
		V = Matrix(Width, Width) + 1;

		// We maintain equality UAV = B, where B is of the form
		// B = [I00;0X0] with I of size rank × rank
		// and X of size (m - rank) × (maxrank - rank).

		rank = 0;
		int maxrank = Width;
		int selector = 1 << rank;

		while (rank < maxrank) {
			bool foundpivot = false;
			for (int i0 = rank; i0 < Height; i0++)
				if (B.Row[i0] & selector) {
					foundpivot = true;
					if (i0 != rank) {
						B.RowSwap(i0, rank);
						U.RowSwap(i0, rank);
					} // move pivot to position (rank, rank)
					for (int i = i0 + 1; i < Height; i++)
						if (B.Row[i] & selector) {
							B.AddRow(rank, i);
							U.AddRow(rank, i);
						}
					V.AddColToCols(rank, B.Row[rank] ^ selector);
					B.Row[rank] = selector;
					rank++;
					selector = 1 << rank;
					break;
				}
			if (!foundpivot) {
				maxrank--;
				if (rank != maxrank) {
					B.ColSwap(rank, maxrank);
					V.ColSwap(rank, maxrank);
				}
			}
		}
	}

	// Matrix inverse.
	Matrix Inv() const {
		int rank;
		Matrix U, V;
		Diagonalize(rank, U, V);
		return V * U;
	}

	// Matrix transpose.
	Matrix Transpose() const {
		Matrix A(Width, Height);
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				A.SetEntry(j, i, GetEntry(i, j));
		return A;
	}

	// Matrix trace.
	bool Trace() const {
		bool sum = false;
		for (int i = 0; i < Height; i++)
			sum ^= GetEntry(i, i);
		return sum;
	}

	// Submatrix with given upper left and lower right corner.
	Matrix SubMatrix(int i1, int j1, int i2, int j2) const {
		Matrix A(i2 - i1, j2 - j1);
		for (int i = i1; i < i2; i++)
			A.Row[i - i1] = (Row[i] & ((1 << j2) - 1)) >> j1;
		return A;
	}

	// Decompose into blocks, with the upper left block of size m×n.
	void Decompose(int m, int n, Matrix& A11, Matrix& A12, Matrix& A21, Matrix& A22) const {
		A11 = SubMatrix(0, 0, m, n);
		A12 = SubMatrix(0, n, m, Width);
		A21 = SubMatrix(m, 0, Height, n);
		A22 = SubMatrix(m, n, Height, Width);
	}

	// Paste a smaller matrix A into the given matrix with the upper left position (i, j).
	void Paste(const Matrix& A, int i, int j) {
		for (int k = 0; k < A.Height; k++) {
			Row[i + k] &= ~(((1 << A.Width) - 1) << j);
			Row[i + k] |= A.Row[k] << j;
		}
	}

	// Block diagonal matrix with upper left block the given matrix and lower right block A.
	Matrix BlockDiag(const Matrix& A) const {
		Matrix B(Height + A.Height, Width + A.Width);
		B.Paste(*this, 0, 0);
		B.Paste(A, Height, Width);
		return B;
	}

	// A well-ordering on the set of matrices of the same size.
	// Rows are ordered as binary numbers, with the first row most significant.
	bool operator <(const Matrix& A) const {
		for (int i = 0; i < Height; i++)
			if (Row[i] < A.Row[i])
				return true;
			else if (Row[i] > A.Row[i])
				return false;
		return false;
	}

	// Count number of zero entries.
	int CountZeros() const {
		int n = 0;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				if (!GetEntry(i, j))
					n++;
		return n;
	}

};

typedef unsigned long long Matrix64;

// Constants and methods to be used in manipulations with 8×8 matrices
// presented by a single 64-bit unsigned integer (matrices of type Matrix64).
namespace Mat64 {

	// Multiples of 8.
	constexpr int N[9] = { 0, 8, 16, 24, 32, 40, 48, 56, 64 };

	// First row and column selector constants.
	constexpr Matrix64 FirstRow = 0xFF, FirstCol = 0x0101010101010101;

	// Convert m×n Matrix (where m,n <= 8) to Matrix64.
	Matrix64 ToMatrix64(const Matrix& A) {
		Matrix64 matrix = 0;
		for (int i = 0; i < A.Height; i++)
			matrix ^= ((Matrix64)A.Row[i]) << N[i];
		return matrix;
	}

	// Convert the upper left m×n block of Matrix64 to Matrix.
	Matrix ToMatrix(Matrix64 A, int m, int n) {
		Matrix B(m, n);
		int selector = (1 << n) - 1;
		for (int i = 0; i < m; i++)
			B.Row[i] = (A >> N[i]) & selector;
		return B;
	}

	// Identity n×n matrix.
	Matrix64 Identity(int n) {
		Matrix64 A = 0;
		for (int i = 0; i < n; i++)
			A ^= 1ull << (i + N[i]);
		return A;
	}

	// Matrix product.
	Matrix64 Product(Matrix64 A, Matrix64 B) {
		Matrix64 C = 0;
		for (int i = 0; i < 8; i++)
			C ^= ((A >> i) & FirstCol) * ((B >> N[i]) & FirstRow);
		return C;
	}

	// Check if AX = B for some X. Begin with start_row-th row.
	bool AXisB(Matrix64 A, Matrix64 B, int start_row = 0) {
		if (start_row) {
			A = A >> N[start_row];
			B = B >> N[start_row];
		}
		while (A) {
			Matrix64 rows = (A & FirstCol) * FirstRow; // rows containing nonzero first elements
			Matrix64 pivotRow = (rows & -(signed long long)rows) * FirstRow; // lowest (pivot) row selector
			A ^= rows & ((pivotRow & A) * FirstCol); // eliminate pivot row and reduce other nonzero rows
			B ^= rows & ((pivotRow & B) * FirstCol); // corresponding operation in B
			A >>= 1;
		}
		return !B;
	}

	// Check if XA = B for some X. Begin with start_col-th column.
	bool XAisB(Matrix64 A, Matrix64 B, int start_col = 0) {
		if (start_col) {
			A >>= start_col;
			B >>= start_col;
			Matrix64 selector = FirstCol * ((1ull << (8 - start_col)) - 1);
			A &= selector;
			B &= selector;
		}
		while (A) {
			Matrix64 cols = (A & FirstRow) * FirstCol; // columns containing nonzero first elements
			Matrix64 pivotCol = (cols & -(signed long long)cols) * FirstCol; // lowest (pivot) column selector
			A ^= cols & ((pivotCol & A) * FirstRow);
			B ^= cols & ((pivotCol & B) * FirstRow);
			A >>= 8;
		}
		return !B;
	}

	// Find rank of a m×n matrix A and invertible U (m×m) and V (n×n) such that UAV=[I0;00].
	void Diagonalize(Matrix64 A, int m, int n, int& rank, Matrix64& U, Matrix64& V) {
		U = 0;
		for (int i = 0; i < m; i++)
			U ^= 1ull << (i + N[i]);
		V = 0;
		for (int i = 0; i < n; i++)
			V ^= 1ull << (i + N[i]);

		// We maintain U, V and the matrix UAV = [I00;0X0;000] with I of size rank × rank
		// and X of size (m - rank) × (maxrank - rank).

		rank = 0;
		int maxrank = n;

		while (rank < maxrank) {
			Matrix64 col = (A >> rank) & FirstCol; // rank-th column of A
			if (col) {
				Matrix64 rankEntry = 1ull << N[rank];
				if ((col & rankEntry) == 0) {
					Matrix64 pivotRow = (col & -(signed long long)col) * FirstRow; // pivot row selector
					Matrix64 rowA = pivotRow & A;
					Matrix64 rowU = pivotRow & U;
					Matrix64 pivotRankRow = FirstRow << N[rank];
					if (rowA == pivotRow)
						A ^= pivotRankRow;
					else
						A ^= pivotRankRow & (rowA / FirstRow);
					if (rowU == pivotRow)
						U ^= pivotRankRow;
					else
						U ^= pivotRankRow & (rowU / FirstRow); // move pivot to (rank,rank)
				}
				else
					col ^= rankEntry;
				A ^= ((A >> N[rank])& FirstRow) * col;
				U ^= ((U >> N[rank])& FirstRow) * col; // rank-th column of A is now (0,...,0,1,0,...,0)

				Matrix64 row = ((A >> N[rank]) & FirstRow) ^ (1ull << rank);
				A ^= row << N[rank];
				V ^= ((V >> rank) & FirstCol) * row;
				rank++;
			}
			else {
				maxrank--;
				if (rank != maxrank) {
					A ^= ((A >> maxrank) & FirstCol) << rank;
					A ^= ((A >> rank) & FirstCol) << maxrank;
					V ^= ((V >> maxrank) & FirstCol) << rank;
					V ^= ((V >> rank) & FirstCol) << maxrank;
				}
			}
		}
	}

	// Inverse of n×n matrix.
	Matrix64 Inv(Matrix64 A, int n) {
		int rank;
		Matrix64 U, V;
		Diagonalize(A, n, n, rank, U, V);
		return Product(V, U);
	}

	// Decompose Matrix into four Matrix64 matrices, with the upper left block of size m×n.
	void Decompose(const Matrix& A, int m, int n, Matrix64& A11, Matrix64& A12, Matrix64& A21, Matrix64& A22) {
		Matrix mA11, mA12, mA21, mA22;
		A.Decompose(m, n, mA11, mA12, mA21, mA22);
		A11 = ToMatrix64(mA11);
		A12 = ToMatrix64(mA12);
		A21 = ToMatrix64(mA21);
		A22 = ToMatrix64(mA22);
	}

	// Assemble four blocks into a single m×n matrix, where the upper left block is of size i×j.
	Matrix Assemble(int m, int n, int i, int j, Matrix64 A11, Matrix64 A12, Matrix64 A21, Matrix64 A22) {
		Matrix mA11, mA12, mA21, mA22;
		mA11 = Mat64::ToMatrix(A11, i, j);
		mA12 = Mat64::ToMatrix(A12, i, n - j);
		mA21 = Mat64::ToMatrix(A21, m - i, j);
		mA22 = Mat64::ToMatrix(A22, m - i, n - j);
		Matrix A(m, n);
		A.Paste(mA11, 0, 0);
		A.Paste(mA12, 0, j);
		A.Paste(mA21, i, 0);
		A.Paste(mA22, i, j);
		return A;
	}

}

// Polynomial matrix of degree at most 1 over
// variables x[i,j] for 0<=i,j<8. Coefficients of such a matrix are linear combinations
// x[i1,j1]+...+x[in,jn]+lambda, where lambda is a scalar (0 or 1).
// Entries in the matrix are presented as Matrix64 (denoting which monomials appear in the entry)
// together with the corresponding entry in another matrix of constants.
class PolyMatrix {

public:

	// Dimensions.
	int Height, Width;

	// Homogeneous part of matrix entries (homogeneous polynomials of degree 1).
	Matrix64 VarPart[MAX_DIM][MAX_DIM];

	// Constant part of matrix entries.
	Matrix ConstPart;

	// Create 0×0 matrix.
	PolyMatrix() {
		Height = 0;
		Width = 0;
	}

	// Create zero m×n matrix.
	PolyMatrix(int m, int n) {
		Height = m;
		Width = n;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				VarPart[i][j] = 0;
		ConstPart = Matrix(m, n);
	}

	// Set a m×n matrix either to zero or a n×n matrix to identity.
	PolyMatrix& operator =(int value) {
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				VarPart[i][j] = 0;
		ConstPart = value;
		return *this;
	}

	// Set a m×n matrix, where m,n<=8, to the canonical matrix with monomial entries.
	// The input string must be "X".
	PolyMatrix& operator =(const std::string& s) {
		if (s != "X")
			return *this;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				VarPart[i][j] = 1ull << (Mat64::N[i] + j);
		ConstPart = 0;
		return *this;
	}

	// Convert a constant matrix to PolyMatrix.
	PolyMatrix(const Matrix& A) {
		Height = A.Height;
		Width = A.Width;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				VarPart[i][j] = 0;
		ConstPart = A;
	}

	// Matrix sum.
	PolyMatrix operator +(const PolyMatrix& A) const {
		PolyMatrix B = *this;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				B.VarPart[i][j] ^= A.VarPart[i][j];
		B.ConstPart = B.ConstPart + A.ConstPart;
		return B;
	}

	// Add a scalar matrix to a n×n matrix.
	PolyMatrix operator +(int value) const {
		PolyMatrix A = *this;
		A.ConstPart = A.ConstPart + value;
		return A;
	}

	// Left multiply by a scalar matrix.
	PolyMatrix LeftMultiply(const Matrix& A) const {
		PolyMatrix B(A.Height, Width);
		for (int i = 0; i < A.Height; i++)
			for (int j = 0; j < Width; j++)
				for (int k = 0; k < Height; k++)
					if (A.GetEntry(i, k))
						B.VarPart[i][j] ^= VarPart[k][j];
		B.ConstPart = A * ConstPart;
		return B;
	}

	// Right multiply by a scalar matrix.
	PolyMatrix RightMultiply(const Matrix& A) const {
		PolyMatrix B(Height, A.Width);
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < A.Width; j++)
				for (int k = 0; k < Width; k++)
					if (A.GetEntry(k, j))
						B.VarPart[i][j] ^= VarPart[i][k];
		B.ConstPart = ConstPart * A;
		return B;
	}

	// Matrix equality. Matrices must be of equal sizes.
	inline bool operator ==(const PolyMatrix& A) const {
		if (ConstPart != A.ConstPart)
			return false;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				if (VarPart[i][j] != A.VarPart[i][j])
					return false;
		return true;
	}

	// Matrix inequality. Matrices must be of equal sizes.
	inline bool operator !=(const PolyMatrix& A) const {
		return !(*this == A);
	}

	// Check if a n×n matrix equals a scalar matrix.
	bool operator ==(int value) const {
		if (ConstPart != value)
			return false;
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++)
				if (VarPart[i][j])
					return false;
		return true;
	}

	// Check if a n×n matrix equals a scalar matrix.
	bool operator !=(int value) const {
		return !(*this == value);
	}

	// Create a list of equations defined by the equation X = 0.
	void EquationList(std::vector<Matrix64>& homParts, std::vector<bool>& constParts) const {
		homParts.clear();
		constParts.clear();
		for (int i = 0; i < Height; i++)
			for (int j = 0; j < Width; j++) {
				homParts.push_back(VarPart[i][j]);
				constParts.push_back(ConstPart.GetEntry(i, j));
			}
	}

	// Perform substitutions given in equations, by replacing the lowest bits with higher bits.
	void Substitution(const std::vector<Matrix64>& equations, const std::vector<bool>& consts) {
		for (int k = 0; k < equations.size(); k++) {
			Matrix64 pivot = equations.at(k) & -(signed long long)equations.at(k);
			for (int i = 0; i < Height; i++)
				for (int j = 0; j < Width; j++)
					if (VarPart[i][j] & pivot) {
						VarPart[i][j] ^= equations.at(k);
						ConstPart.SetEntry(i, j, ConstPart.GetEntry(i, j) ^ consts.at(k));
					}
		}
	}

};

// Class describing a subset of {0,...,N-1} for 1<=N<=MAX_DIM.
// Subsets are ordered by size, and same-sized subsets are ordered lexicographically.
// For example, {3} < {1,4} < {1,5} < {2,3}.
class Subset {
public:

	int N;
	int Order; // size of the subset
	int Element[MAX_DIM]; // elements of the subset

	// Create empty set.
	Subset(int n) {
		N = n;
		Order = 0;
	}

	// The minimal subset of order k.
	Subset(int n, int k) {
		N = n;
		Order = k;
		for (int i = 0; i < k; i++)
			Element[i] = i;
	}

	// Create a set from a string. Example: "0,3,11" gives {0,3,11}.
	Subset(int n, std::string s) {
		N = n;
		Order = 0;
		if (s != "") s += ",";
		while (s.find(",") != std::string::npos) {
			std::string numberString = s.substr(0, s.find(","));
			Element[Order] = stoi(numberString);
			s = s.substr(s.find(",") + 1);
			Order++;
		}
	}

	// Increment subset by replacing it with its successor.
	// Return false if successor does not exist.
	bool Inc() {
		for (int i = Order - 1; i >= 0; i--)
			if (Element[i] != N - Order + i) {
				Element[i]++;
				for (int j = i + 1; j < Order; j++)
					Element[j] = Element[i] + j - i;
				return true;
			}
		if (Order < N) {
			Order++;
			for (int i = 0; i < Order; i++)
				Element[i] = i;
			return true;
		}
		return false;
	}

	// Convert a subset to a string. For example, {0,3,11} is converted to "0,3,11".
	std::string ToString() const {
		std::string s;
		for (int i = 0; i < Order; i++) {
			s += std::to_string(Element[i]);
			if (i != Order - 1) s += ",";
		}
		return s;
	}

	// Permutation matrix associated with a subset of {0,...,N-1} of order k.
	// Heights of first k columns are given by the subset in increasing order,
	// and the remaining columns form an increasing sequence as well.
	// Example: subset {1,3} with n=4 is associated with P=[0010;1000;0001;0100]
	// (1st column has height 2, second height 4,
	// and the remaining two columns have heights 1 and 3.
	Matrix PermutationMatrix() const {
		Matrix P(N, N);
		for (int j = 0; j < Order; j++)
			P.SetEntry(Element[j], j, true);
		for (int j = Order; j < N; j++) {
			int i = j - Order;
			for (int k = 0; k < Order; k++)
				if (Element[k] <= j - Order + k) i++;
			P.SetEntry(i, j, true);
		}
		return P;
	}

	// Check if the set contains an element.
	bool Contains(int value) const {
		for (int i = 0; i < Order; i++)
			if (Element[i] == value)
				return true;
		return false;
	}

	// Complement in the set {0,...,N-1}.
	Subset Complement() const {
		Subset A(N);
		int order = 0;
		for (int i = 0; i < N; i++)
			if (!Contains(i)) {
				A.Element[order] = i;
				order++;
			}
		A.Order = order;
		return A;
	}

};

typedef std::vector<Matrix> Matrices;

// Methods for finding nil-clean decompositions of index 3 of matrices over F_2.
namespace NilClean {

	// Given n×n matrix A and 0<p<n, find idempotents E of rank p in block form
	// E=[I0;XI][IY;00][I0;XI] (with the upper left block of size p×p),
	// such that (A-E)^3=0. Heights of columns in block X are limited by heights[].
	// (E.g. heights[] = {1,1,...,1} means that X can have only the top nonzero row.)
	// If find_one is set, find only one idempotent. Append idempotents to the vector idempotents.
	// Use matrices presented by single 64-bit integers. Only for use if p,n-p<=8.
	// If rand != 0 then use random heuristic with rand attempts.
	void FindIdempotents64(const Matrix& A, int p, int heights[MAX_DIM], bool find_one, Matrices& idempotents, int rand) {
		int n = A.Height;
		int q = n - p;

		Matrix B = (A + 1) * (A + 1) * (A + 1);
		Matrix C = A * A * A;

		// Allowed nonzero entries of X. Only in use if rand > 0.
		Matrix64 SelectorX = 0;
		for (int j = 0; j < p; j++)
			for (int i = 0; i < heights[j]; i++)
				SelectorX ^= 1ull << (Mat64::N[i] + j);

		// Random seed.
		std::random_device rd;

		// Random generator.
		std::default_random_engine generator(rd());

		// Random distribution.
		std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF);

		// The outer loop will run through all matrices P=[I0;XI] with X having given heights.
		// Initializing X in use only if rand = 0.
		Matrix64 X = 0;

		// We are looking for Y such that (PAP + [IY;00])^3 = 0.

		// CounterX determines which entry of X is about to change. Only in use if rand = 0.
		Matrix64 CounterX = 0;

		// Blocks of matrices A, B, C.
		Matrix64 A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
		Mat64::Decompose(A, p, p, A11, A12, A21, A22);
		Mat64::Decompose(B, p, p, B11, B12, B21, B22);
		Mat64::Decompose(C, p, p, C11, C12, C21, C22);

		// Blocks of matrices PAP, PBP, PCP. Maintained if rand = 0.
		// Initialization relevant only if rand = 0.
		Matrix64 XA11, XA12, XA21, XA22, XB11, XB12, XB21, XB22, XC11, XC12, XC21, XC22;
		Mat64::Decompose(A, p, p, XA11, XA12, XA21, XA22);
		Mat64::Decompose(B, p, p, XB11, XB12, XB21, XB22);
		Mat64::Decompose(C, p, p, XC11, XC12, XC21, XC22);

		// Counter of random attempts completed, irrelevant if rand = 0.
		int attemptsX = 0;

		bool counterXfound = true;
		while (counterXfound) {
			if (rand) {
				// Random choice of X.
				X = SelectorX & distribution(generator);

				// Compute blocks of PAP, PBP, PCP.
				XA11 = A11 ^ Mat64::Product(A12, X);
				XA21 = A21 ^ Mat64::Product(X, XA11) ^ Mat64::Product(A22, X);
				XA22 = A22 ^ Mat64::Product(X, A12);

				XB11 = B11 ^ Mat64::Product(B12, X);
				XB21 = B21 ^ Mat64::Product(X, XB11) ^ Mat64::Product(B22, X);
				XB22 = B22 ^ Mat64::Product(X, B12);

				XC11 = C11 ^ Mat64::Product(C12, X);
				XC21 = C21 ^ Mat64::Product(X, XC11) ^ Mat64::Product(C22, X);
				XC22 = C22 ^ Mat64::Product(X, C12);
			}
			if (Mat64::AXisB(XA21, XB21) && Mat64::XAisB(XA21, XC21)) {
				int rank;
				Matrix64 U, V;
				Mat64::Diagonalize(XA21, q, p, rank, U, V);
				if (Mat64::AXisB(Mat64::Product(U, XC21), Mat64::Product(U, XC22), rank) &&
					Mat64::XAisB(Mat64::Product(XB21, V), Mat64::Product(XB11, V), rank)) {

					Matrix64 Uinv, Vinv;
					Uinv = Mat64::Inv(U, q);
					Vinv = Mat64::Inv(V, p);

					// We maintain blocks of Q, R=Q^2 and S=Q^3.
					Matrix64 Q11, Q12, Q21, Q22, R11, R12, R21, R22, S11, S12, S21, S22;

					Q11 = Mat64::Product(Mat64::Product(Vinv, XA11), V) ^ Mat64::Identity(p);
					Q12 = Mat64::Product(Mat64::Product(Vinv, XA12), Uinv) ^ Mat64::Product(Mat64::Product(U, XB21 ^ Mat64::Product(XA21, XA11)), V);
					Q21 = Mat64::Identity(rank);
					Q22 = Mat64::Product(Mat64::Product(U, XA22), Uinv);

					R11 = Mat64::Product(Q11, Q11) ^ Mat64::Product(Q12, Q21);
					R12 = Mat64::Product(Q11, Q12) ^ Mat64::Product(Q12, Q22);
					R21 = Mat64::Product(Q21, Q11) ^ Mat64::Product(Q22, Q21);
					R22 = Mat64::Product(Q21, Q12) ^ Mat64::Product(Q22, Q22);

					S11 = Mat64::Product(Q11, R11) ^ Mat64::Product(Q12, R21);
					S12 = Mat64::Product(Q11, R12) ^ Mat64::Product(Q12, R22);
					S21 = Mat64::Product(Q21, R11) ^ Mat64::Product(Q22, R21);
					S22 = Mat64::Product(Q21, R12) ^ Mat64::Product(Q22, R22);

					// CounterY determines which entry of Q is about to change.
					Matrix64 CounterY = 0;

					bool counterYfound = true;
					while (counterYfound) {
						if (S11 == 0 && S12 == 0 && S21 == 0 && S22 == 0) {
							Matrix Q = Mat64::Assemble(n, n, p, p, Q11, Q12, Q21, Q22);
							Matrix W = Mat64::Assemble(n, n, p, p, Vinv, 0, 0, U);
							Matrix P = Mat64::Assemble(n, n, p, p, 0, 0, X, 0) + 1;
							idempotents.push_back(A + P * W.Inv() * Q * W * P);
							if (find_one) return;
						}
						counterYfound = false;
						for (int i = 0, j0 = rank; i < p && !counterYfound; i++) {
							if (i == rank) j0 = 0;
							for (int j = j0; j < q; j++) {
								Matrix64 selector = 1ull << (Mat64::N[i] + j);
								if (CounterY & selector)
									CounterY ^= selector;
								else {
									CounterY ^= selector;

									S11 ^= ((Q11 >> i)& Mat64::FirstCol)* ((Q21 >> Mat64::N[j])& Mat64::FirstRow);
									S12 ^= ((Q11 >> i)& Mat64::FirstCol)* ((Q22 >> Mat64::N[j])& Mat64::FirstRow);
									S21 ^= ((Q21 >> i)& Mat64::FirstCol)* ((Q21 >> Mat64::N[j])& Mat64::FirstRow);
									S22 ^= ((Q21 >> i)& Mat64::FirstCol)* ((Q22 >> Mat64::N[j])& Mat64::FirstRow);

									S11 ^= ((R21 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
									S12 ^= ((R22 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];

									R12 ^= ((Q11 >> i)& Mat64::FirstCol) << j;
									R22 ^= ((Q21 >> i)& Mat64::FirstCol) << j;

									R11 ^= ((Q21 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
									R12 ^= ((Q22 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];

									S12 ^= ((R11 >> i)& Mat64::FirstCol) << j;
									S22 ^= ((R21 >> i)& Mat64::FirstCol) << j;

									Q12 ^= selector;

									counterYfound = true;
									break;
								}
							}
						}
					}
				}
			}
			if (rand) {
				attemptsX++;
				counterXfound = attemptsX < rand;
			}
			else {
				counterXfound = false;
				for (int j = 0; j < p && !counterXfound; j++)
					for (int i = 0; i < heights[j]; i++)
						if (CounterX & (1ull << (Mat64::N[i] + j)))
							CounterX ^= 1ull << (Mat64::N[i] + j);
						else {
							CounterX ^= 1ull << (Mat64::N[i] + j);
							X ^= 1ull << (Mat64::N[i] + j);

							XA21 ^= ((XA11 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XA22 ^= ((XA12 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XA11 ^= ((XA12 >> i)& Mat64::FirstCol) << j;
							XA21 ^= ((XA22 >> i)& Mat64::FirstCol) << j;

							XB21 ^= ((XB11 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XB22 ^= ((XB12 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XB11 ^= ((XB12 >> i)& Mat64::FirstCol) << j;
							XB21 ^= ((XB22 >> i)& Mat64::FirstCol) << j;

							XC21 ^= ((XC11 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XC22 ^= ((XC12 >> Mat64::N[j])& Mat64::FirstRow) << Mat64::N[i];
							XC11 ^= ((XC12 >> i)& Mat64::FirstCol) << j;
							XC21 ^= ((XC22 >> i)& Mat64::FirstCol) << j;

							counterXfound = true;
							break;
						}
			}
		}
	}

	// Given n×n matrix A and 0<p<n, find idempotents E of rank p in block form
	// E=[I0;XI][IY;00][I0;XI] (with the upper left block of size p×p),
	// such that (A-E)^3=0. Heights of columns in block X are limited by heights[].
	// (E.g. heights[] = {1,1,...,1} means that X can have only the top nonzero row.)
	// If find_one is set, find only one idempotent. Append idempotents to the vector idempotents.
	// If rand != 0 then then use random heuristic with rand attempts.
	void FindIdempotents(const Matrix& A, int p, int heights[MAX_DIM], bool find_one, Matrices& idempotents, int rand) {
		int n = A.Height;
		int q = n - p;

		Matrix B = (A + 1) * (A + 1) * (A + 1);
		Matrix C = A * A * A;

		// The outer loop will run through all matrices P=[I0;XI] with X having given heights.
		// Maintained only if rand = 0.
		Matrix P = Matrix(n, n) + 1;

		// We are looking for Y such that (PAP + [IY;00])^3 = 0.

		// CounterX determines which entry of X is about to change. Relevant if rand = 0.
		Matrix CounterX(n - p, p);

		// Blocks of matrices A, B, C.
		Matrix A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
		A.Decompose(p, p, A11, A12, A21, A22);
		B.Decompose(p, p, B11, B12, B21, B22);
		C.Decompose(p, p, C11, C12, C21, C22);

		// Blocks of matrices PAP, PBP, PCP. Maintained only if rand = 0.
		Matrix XA11, XA12, XA21, XA22, XB11, XB12, XB21, XB22, XC11, XC12, XC21, XC22;
		A.Decompose(p, p, XA11, XA12, XA21, XA22);
		B.Decompose(p, p, XB11, XB12, XB21, XB22);
		C.Decompose(p, p, XC11, XC12, XC21, XC22);

		// Allowed nonzero entries of P. Only in use if rand > 0.
		Matrix SelectorP(n, n);
		for (int j = 0; j < p; j++)
			for (int i = 0; i < heights[j]; i++)
				SelectorP.SetEntry(i + p, j, true);

		// Random seed.
		std::random_device rd;

		// Random generator.
		std::default_random_engine generator(rd());

		// Random distribution.
		std::uniform_int_distribution<int> distribution(0, 0xFFFFFFFF);

		// Counter of random attempts completed, irrelevant if rand = 0.
		int attemptsX = 0;

		bool counterXfound = true;
		while (counterXfound) {
			if (rand) {
				for (int i = p; i < n; i++)
					P.Row[i] ^= SelectorP.Row[i] & distribution(generator);

				// Compute blocks of PAP, PBP, PCP.
				(P * A * P).Decompose(p, p, XA11, XA12, XA21, XA22);
				(P * B * P).Decompose(p, p, XB11, XB12, XB21, XB22);
				(P * C * P).Decompose(p, p, XC11, XC12, XC21, XC22);
			}
			if (XB21.IsAX(XA21) && XC21.IsXA(XA21)) {
				int rank;
				Matrix U, V;
				XA21.Diagonalize(rank, U, V);
				if ((U * XC22).IsAX(U * XC21, rank) && (XB11 * V).IsXA(XB21 * V, rank)) {
					Matrix Q(n, n);
					Q.Paste(V.Inv() * XA11 * V + 1, 0, 0);
					Q.Paste((V.Inv() * XA12 * U.Inv()).SubMatrix(0, 0, rank, rank) + (U * (XB21 + XA21 * XA11) * V).SubMatrix(0, 0, rank, rank), 0, p);
					Q.Paste(Matrix(rank, rank) + 1, p, 0);
					Q.Paste(U * XA22 * U.Inv(), p, p);

					// In the inner loop, we maintain Q and R=Q^2.
					Matrix R = Q * Q;

					// CounterY determines which entry of Q is about to change.
					Matrix CounterY(n, n);

					bool counterYfound = true;
					while (counterYfound) {
						if (Q.ZeroProduct(R)) {
							Matrix W = V.Inv().BlockDiag(U);
							idempotents.push_back(A + P * W.Inv() * Q * W * P);
							if (find_one) return;
						}
						counterYfound = false;
						for (int i = 0, j0 = p + rank; i < p && !counterYfound; i++) {
							if (i == rank) j0 = p;
							for (int j = j0; j < n; j++)
								if (CounterY.GetEntry(i, j)) CounterY.SetEntry(i, j, false);
								else {
									CounterY.SetEntry(i, j, true);
									R.AddCol(i, Q, j);
									R.AddRow(j, Q, i);
									Q.Row[i] ^= 1 << j;
									counterYfound = true;
									break;
								}
						}
					}
				}
			}
			if (rand) {
				attemptsX++;
				counterXfound = attemptsX < rand;
			}
			else {
				counterXfound = false;
				for (int j = 0; j < p && !counterXfound; j++)
					for (int i = 0; i < heights[j]; i++)
						if (CounterX.GetEntry(i, j))
							CounterX.SetEntry(i, j, false);
						else {
							CounterX.SetEntry(i, j, true);
							P.Row[i + p] ^= 1 << j;

							XA21.AddRow(j, XA11, i);
							XA22.AddRow(j, XA12, i);
							XA11.AddCol(i, XA12, j);
							XA21.AddCol(i, XA22, j);

							XB21.AddRow(j, XB11, i);
							XB22.AddRow(j, XB12, i);
							XB11.AddCol(i, XB12, j);
							XB21.AddCol(i, XB22, j);

							XC21.AddRow(j, XC11, i);
							XC22.AddRow(j, XC12, i);
							XC11.AddCol(i, XC12, j);
							XC21.AddCol(i, XC22, j);

							counterXfound = true;
							break;
						}
			}
		}
	}

	// Given n×n matrix A and a subset M of {0,...,n-1}, find idempotents E
	// such that (A-E)^3=0 and such that P^T*E*P is of the form [I0;XI][IY;00][I0;XI],
	// where P is a permuation matrix associated with M. If find_one is set,
	// it finds only one idempotent. Append idempotents to the vector idempotents.
	// If rand != 0, use random heuristic with rand attempts.
	void FindIdempotents(const Matrix& A, const Subset& M, bool find_one, Matrices& idempotents, int rand) {
		int n = M.N;
		int p = M.Order;
		if (((p & 1) == 1) ^ A.Trace())
			return;
		Matrix P = M.PermutationMatrix();
		int heights[MAX_DIM];
		for (int j = 0; j < p; j++)
			heights[j] = M.Element[j] - j;
		Matrices new_idempotents;
		Matrix PTAP = P.Transpose() * A * P;
		if (p <= 8 && n - p <= 8)
			FindIdempotents64(PTAP, p, heights, find_one, new_idempotents, rand);
		else
			FindIdempotents(PTAP, p, heights, find_one, new_idempotents, rand);
		for (Matrix E : new_idempotents)
			idempotents.push_back(P * E * P.Transpose());
	}

	// Given n×n matrix A, find idempotents E such that (A-E)^3=0.
	// If show_steps is set, write steps in the procedure. Begin with the subset M0.
	// If find_one is set, then find just one E.
	// Append idempotents to the vector idempotents.
	// If rand != 0, use random heuristic with rand attempts.
	// If rand_sets != 0 then rand != 0 and use random heuristic with rand_sets sets.
	void FindIdempotents(const Matrix& A, bool find_one, const Subset& M0, Matrices& idempotents, int rand, int rand_sets, bool show_steps = true) {
		if (rand_sets) {

			// Random seed.
			std::random_device rd;

			// Random generator.
			std::default_random_engine generator(rd());

			// Random distribution.
			std::uniform_int_distribution<int> distribution(0, (1 << A.Height) - 1);

			for (int i = 0; i < rand_sets; i++) {
				int set_binary = distribution(generator);
				std::string set_string = "";
				for (int k = 0; k < A.Height; k++)
					if (set_binary & (1 << k)) {
						if (set_string != "")
							set_string += ",";
						set_string += std::to_string(k);
					}
				Subset M(A.Height, set_string);
				if (show_steps)
					std::cout << "Set = {" << M.ToString() << "}\n";
				FindIdempotents(A, M, find_one, idempotents, rand);
				if (find_one && !idempotents.empty()) return;
			}
			return;
		}
		Subset M = M0;
		do {
			if (show_steps)
				std::cout << "Set = {" << M.ToString() << "}\n";
			FindIdempotents(A, M, find_one, idempotents, rand);
			if (find_one && !idempotents.empty()) return;
		} while (M.Inc());
	}

}

// Methods for creating companion matrices over F_2.
namespace Companion {

	// Companion matrix or modified companion matrix associated with coefficients in the string
	// coefficients. Example: CompanionMatrix("11100", false) = [00001;10001;01001;00100;00010],
	// CompanionMatrix("11100", true) = [10001;10001;01101;00100;00011].
	// Modified companion matrix is companion matrix + diag(1,0,1,0,...) (ending with 0 or 1).
	Matrix CompanionMatrix(const std::string& coefficients, bool modified = false) {
		int n = coefficients.length();
		Matrix A(n, n);
		for (int i = 1; i < n; i++)
			A.SetEntry(i, i - 1, true);
		for (int i = 0; i < n; i++)
			A.SetEntry(i, n - 1, coefficients[i] == '1');
		if (modified)
			for (int i = 0; 2 * i < n; i++)
				A.SetEntry(2 * i, 2 * i, !A.GetEntry(2 * i, 2 * i));
		return A;
	}

	// Block diagonal matrix with diagonal companion matrices.
	// Parameter "m" determines if a given block is modified or not.
	// Example: "11,m10" gives block diagonal 4×4 matrix diag([01;11],[11;10]).
	Matrix BlockCompanionMatrix(const std::string& coefficients) {
		std::size_t pos = coefficients.find(",");
		if (pos == std::string::npos) {
			if (coefficients[0] == 'm')
				return CompanionMatrix(coefficients.substr(1), true);
			else
				return CompanionMatrix(coefficients);
		}
		Matrix A = BlockCompanionMatrix(coefficients.substr(0, pos));
		Matrix B = BlockCompanionMatrix(coefficients.substr(pos + 1));
		return A.BlockDiag(B);
	}

	// Return the list of all block diagonal matrices with blocks of prescribed dimensions.
	// For example, sizes can be "3,2,3,1".
	Matrices BlockCompanionMatrices(const std::string& sizes) {
		std::size_t pos = sizes.find(",");
		if (pos != std::string::npos) {
			Matrices list1 = BlockCompanionMatrices(sizes.substr(0, pos));
			Matrices list2 = BlockCompanionMatrices(sizes.substr(pos + 1));
			Matrices list3;
			for (Matrix A : list1)
				for (Matrix B : list2)
					list3.push_back(A.BlockDiag(B));
			return list3;
		}
		int n = std::stoi(sizes);
		Matrices list;
		for (int counter = 0; counter < (1 << n); counter++) {
			std::string s;
			for (int i = 0; i < n; i++)
				s += (counter & (1 << i)) == 0 ? "0" : "1";
			list.push_back(BlockCompanionMatrix(s));
		}
		return list;
	}

}

// Methods for finding matrices over F_2 with common entries.
namespace CommonEntries {

	// Given lists l1,...,ln of matrices and column indices, find all matrices A1 in l1
	// such that there exist matrices A2,...,An in the other lists
	// having the same columns with given indices as A1. Append solutions to the list solutions.
	void CommonColumnsFinder(const std::vector<Matrices>& lists, const Subset& columns, Matrices& solutions) {
		if (lists.empty())
			return;
		std::vector<std::set<Matrix>> list_sets;
		for (int k = 1; k < lists.size(); k++) {
			std::set<Matrix> list_set;
			for (Matrix A : lists.at(k)) {
				Matrix B(A.Height, columns.Order);
				for (int i = 0; i < B.Height; i++)
					for (int j = 0; j < B.Width; j++)
						B.SetEntry(i, j, A.GetEntry(i, columns.Element[j]));
				list_set.insert(B);
			}
			list_sets.push_back(list_set);
		}
		for (int k = 0; k < lists.at(0).size(); k++) {
			Matrix A = lists.at(0).at(k);
			Matrix B(A.Height, columns.Order);
			for (int i = 0; i < B.Height; i++)
				for (int j = 0; j < B.Width; j++)
					B.SetEntry(i, j, A.GetEntry(i, columns.Element[j]));
			bool success = true;
			for (int i = 0; i < list_sets.size() && success; i++)
				if (list_sets.at(i).find(B) == list_sets.at(i).end()) {
					success = false;
					break;
				}
			if (success)
				solutions.push_back(A);
		}
	}

	// Given nonempty lists l1,...,lk of m×n matrices, find indices i1,...,ik such that
	// matrices l1.at(i1),...,lk.at(ik) have maximal number of common entries.
	// For example, if l1={[01;11],[10;10]} and l2={[10,00]} then
	// i1=1 and i2=0 since [10;10] and [10;00] have three common entries (which is maximal possible).
	std::vector<int> MaxCommonEntriesFinder(const std::vector<Matrices>& matrixLists, int& common_entries) {
		if (matrixLists.empty()) {
			common_entries = 0;
			std::vector<int> max_indices;
			return max_indices;
		}
		int m = matrixLists.at(0).at(0).Height;
		int n = matrixLists.at(0).at(0).Width;
		std::vector<int> max_indices;
		for (int i = 0; i < matrixLists.size(); i++)
			max_indices.push_back(0); // indices of matrices with max common entries
		std::vector<int> indices;
		for (int i = 0; i < matrixLists.size(); i++)
			indices.push_back(0); // looping indices
		common_entries = -1;
		bool finished = false;
		while (!finished) {
			int com_entries = 0;
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++) {
					bool value = matrixLists.at(0).at(indices.at(0)).GetEntry(i, j);
					bool is_common = true;
					for (int k = 1; k < matrixLists.size(); k++)
						if (matrixLists.at(k).at(indices.at(k)).GetEntry(i, j) != value) {
							is_common = false;
							break;
						}
					if (is_common)
						com_entries++; // count common entries
				}
			if (com_entries > common_entries) {
				common_entries = com_entries;
				for (int k = 0; k < matrixLists.size(); k++)
					max_indices.at(k) = indices.at(k);
			}
			finished = true;
			for (int k = 0; k < matrixLists.size(); k++)
				if (indices.at(k) == matrixLists.at(k).size() - 1)
					indices.at(k) = 0;
				else {
					indices.at(k)++;
					finished = false;
					break;
				}
		}
		return max_indices;
	}

}

// User interface.
namespace UI {

	// Send message.
	void Message(const std::string& text) {
		std::cout << text << "\n";
	}

	// Prompt for yes/no answer.
	bool YesNoQuestion(const std::string& question) {
		std::cout << question;
		std::string s;
		do {
			std::getline(std::cin, s);
			if (s != "y" && s != "n")
				std::cout << "Invalid command. Enter y/n: ";
		} while (s != "y" && s != "n");
		return s == "y";
	}

	// Prompt for a string.
	std::string EnterString(const std::string& query) {
		std::cout << query;
		std::string s;
		std::getline(std::cin, s);
		return s;
	}

	// Read a file containing a list of matrices. Return false if the file is not accessible.
	bool FileRead(const std::string& filename, Matrices& matrices) {
		std::ifstream file(filename);
		if (!file.is_open()) return false;
		std::string line;
		while (getline(file, line))
			matrices.push_back(Matrix(line));
		file.close();
		return true;
	}

	// Write matrices into a file. Return false if the file is not accessible.
	bool FileWrite(const std::string& filename, const Matrices& matrices, bool multi_line = false) {
		std::ofstream file(filename);
		if (!file.is_open()) return false;
		for (Matrix A : matrices)
			file << (multi_line ? A.ToMultiLineString() : A.ToString()) << "\n";
		file.close();
		return true;
	}

	// Prompt for a command. Return the index of the selected command beginning with 1.
	// The number of commands must not exceed 99.
	int Command(const std::vector<std::string>& commands) {
		std::cout << "\nCommands:\n";
		for (int i = 0; i < commands.size(); i++) {
			if (i < 9)
				std::cout << " ";
			std::cout << i + 1 << " " << commands.at(i) << "\n";
		}
		while (true) {
			std::string command = EnterString("Command: ");
			for (int i = 0; i < commands.size(); i++)
				if (std::to_string(i + 1) == command)
					return i + 1;
			Message("Invalid command.");
		}
	}

}

int main() {
	std::vector<std::string> commands = {
		"Find and save idempotents E such that (A-E)^3=0",
		"Display matrices from a file",
		"Convert matrices in a file to a user-friendly multi-line form",
		"Select matrices from a file",
		"Sort matrices with respect to number of zeros",
		"Compute chain dimensions",
		"Find all irreducible polynomials of given degree",
		"Bitwise AND or OR for a given set of matrices",
		"Find matrices in lists with maximal number of common entries",
		"Find maximal number of common columns in lists of matrices",
		"Find all matrices with common columns",
		"Find nil-clean decompositions for all companion matrices",
		"Exit"
	};
	UI::Message("Nil-clean matrices over F_2");
	UI::Message("Author: Janez Ster");
	while (true) {
		int command = UI::Command(commands);
		if (command == 13)
			return 0;
		else if (command == 12) {
			std::string stringA = UI::EnterString("Coefficients of companion matrices of A (example: m1xx1,xx1): ");
			int xcount = std::count(stringA.begin(), stringA.end(), 'x');
			int n = std::count(stringA.begin(), stringA.end(), '0') + std::count(stringA.begin(), stringA.end(), '1') + xcount;
			if (n == 0 || n > MAX_DIM) {
				UI::Message("Illegal matrix size.");
				continue;
			}
			std::string file_base_name = UI::EnterString("Base name of the output files: ");
			std::string starting_set = UI::EnterString("Starting set: ");
			bool find_one = UI::YesNoQuestion("Find one idempotent for each matrix (y/n)? ");
			bool random = UI::YesNoQuestion("Use random heuristic? ");
			int rand = 0;
			if (random)
				rand = std::stoi(UI::EnterString("Number of tries: "));
			std::vector<bool> coefficients;
			for (int i = 0; i < xcount; i++)
				coefficients.push_back(false);
			bool increased = true;
			int counter = 0;
			while (increased) {
				int index = 0;
				std::string stringB = stringA;
				for (int i = 0; i < stringA.length(); i++)
					if (stringA.at(i) == 'x')
						stringB[i] = coefficients.at(index++) ? '1' : '0';
				Matrix A = Companion::BlockCompanionMatrix(stringB);
				std::cout << "\n" << A.ToMultiLineString() << "\n";
				Subset M0(A.Height, starting_set);
				Matrices idempotents;
				NilClean::FindIdempotents(A, find_one, M0, idempotents, rand, false);
				std::cout << "Number of idempotents found: " << idempotents.size() << "\n";
				std::cout << "Checking ...";
				for (Matrix E : idempotents) {
					Matrix Q = A + E;
					if ((E * E != E) || (Q * Q * Q != 0)) {
						std::cout << " FAILED!!! Bug in the program!\n";
						return 0;
					}
				}
				UI::Message(" Done.");
				std::string file_name = file_base_name + std::to_string(counter++) + ".txt";
				UI::FileWrite(file_name, idempotents);
				std::cout << "Writing file " << file_name << " complete.\n";
				increased = false;
				for (int i = 0; i < coefficients.size(); i++)
					if (coefficients.at(i))
						coefficients.at(i) = false;
					else {
						coefficients.at(i) = true;
						increased = true;
						break;
					}
			}
		}
		else if (command == 11) {
			UI::Message("Given nonempty lists l1, ..., lk of mxn matrices and column indices, find all matrices A1 in l1");
			UI::Message("such that there exist matrices A2,....,Ak in l1,...,lk having common columns.");
			std::vector<Matrices> matrixLists;
			std::string s;
			do {
				s = UI::EnterString("File no. " + std::to_string(matrixLists.size() + 1) + " (Enter to end): ");
				if (s != "") {
					Matrices matrices;
					if (UI::FileRead(s, matrices))
						matrixLists.push_back(matrices);
					else
						UI::Message("File could not be read.");
				}
			} while (s != "");
			if (matrixLists.empty() || matrixLists.at(0).empty()) {
				UI::Message("Empty list.");
				continue;
			}
			int n = matrixLists.at(0).at(0).Width;
			Subset indices(n, UI::EnterString("Columns indices: "));
			Matrices solutions;
			CommonEntries::CommonColumnsFinder(matrixLists, indices, solutions);
			std::cout << "Number of solutions found: " << solutions.size();
			if (!solutions.empty() && UI::YesNoQuestion("Write solutions? ")) {
				if (UI::FileWrite(UI::EnterString("File name: "), solutions, false))
					UI::Message("Done.");
				else
					UI::Message("Error writing file.");
			}
		}
		else if (command == 10) {
			UI::Message("Given nonempty lists l1, ..., lk of mxn matrices, find maximal integer k such that");
			UI::Message("there exist matrices A1,....,Ak in l1,...,lk having k common columns.");
			std::vector<Matrices> matrixLists;
			std::string s;
			do {
				s = UI::EnterString("File no. " + std::to_string(matrixLists.size() + 1) + " (Enter to end): ");
				if (s != "") {
					Matrices matrices;
					if (UI::FileRead(s, matrices))
						matrixLists.push_back(matrices);
					else
						UI::Message("File could not be read.");
				}
			} while (s != "");
			if (matrixLists.empty() || matrixLists.at(0).empty()) {
				UI::Message("Empty list.");
				continue;
			}
			int n = matrixLists.at(0).at(0).Width;
			Subset MissingIndices(n);
			do {
				Subset indices = MissingIndices.Complement();
				Matrices solutions;
				CommonEntries::CommonColumnsFinder(matrixLists, indices, solutions);
				if (!solutions.empty()) {
					std::cout << "Found matrix:\n" << solutions.at(0).ToMultiLineString() << "\n";
					std::cout << "Indices of common columns: " << indices.ToString() << "\n";
					break;
				}
			} while (MissingIndices.Inc() && MissingIndices.Order != n);
			if (MissingIndices.Order == n)
				UI::Message("No common columns found.");
		}
		else if (command == 9) {
			UI::Message("Given lists l1, ..., lk of mxn matrices, find matrices Ai in li such that");
			UI::Message("the number of common entries in A1, ..., Ak is maximal possible.");
			std::vector<Matrices> matrixLists;
			std::string s;
			do {
				s = UI::EnterString("File no. " + std::to_string(matrixLists.size() + 1) + " (Enter to end): ");
				if (s != "") {
					Matrices matrices;
					if (UI::FileRead(s, matrices))
						matrixLists.push_back(matrices);
					else
						UI::Message("File could not be read.");
				}
			} while (s != "");
			if (matrixLists.empty()) {
				UI::Message("Empty list.");
				continue;
			}
			int common_entries;
			std::vector<int> indices = CommonEntries::MaxCommonEntriesFinder(matrixLists, common_entries);
			UI::Message("Found indices: ");
			for (int k = 0; k < indices.size(); k++) {
				if (k != 0)
					std::cout << ",";
				std::cout << indices.at(k) + 1;
			}
			std::cout << "\nCommon entries: " << common_entries << "\n";
			if (UI::YesNoQuestion("Write matrices into a file? ")) {
				Matrices matrices;
				for (int k = 0; k < indices.size(); k++)
					matrices.push_back(matrixLists.at(k).at(indices.at(k)));
				if (UI::FileWrite(UI::EnterString("File name: "), matrices))
					UI::Message("Done.");
				else
					UI::Message("Error writing file.");
			}
		}
		else if (command == 8) {
			Matrices matrices;
			if (!UI::FileRead(UI::EnterString("Input file name: "), matrices)) {
				UI::Message("Error reading file.");
				continue;
			}
			if (matrices.empty()) {
				UI::Message("Empty matrix list.");
				continue;
			}
			bool operation = UI::YesNoQuestion("Bitwise operation (Yes=AND, No=OR): ");
			Matrix A = matrices.at(0);
			for (Matrix X : matrices)
				for (int i = 0; i < A.Height; i++)
					for (int j = 0; j < A.Width; j++)
						if (operation)
							A.SetEntry(i, j, A.GetEntry(i, j) && X.GetEntry(i, j));
						else
							A.SetEntry(i, j, A.GetEntry(i, j) || X.GetEntry(i, j));
			std::cout << "Output matrix: " << A.ToMultiLineString() << "\n";
		}
		else if (command == 7) {
			int degree = std::stoi(UI::EnterString("Polynomial degree: "));
			Polynomial candidate = Polynomial(degree);
			while (candidate.Degree() == degree) {
				if (candidate.Irreducible())
					std::cout << candidate.ToString() << "\n";
				candidate.Inc();
			}
			UI::Message("Finished.");
		}
		else if (command == 6) {
			UI::Message("Given matrix A, a set of matrices E, and U = Lin(e1,...,ek),");
			UI::Message("compute dimensions of U, EU, AEU, A^2EU, A^3EU, EA^3EU.");
			UI::Message("");
			std::string stringA = UI::EnterString("Coefficients of companion matrices of A (example: m1111,101): ");
			int n = std::count(stringA.begin(), stringA.end(), '0') + std::count(stringA.begin(), stringA.end(), '1');
			if (n == 0 || n > MAX_DIM) {
				UI::Message("Illegal matrix size.");
				continue;
			}
			Matrix A = Companion::BlockCompanionMatrix(stringA);
			std::cout << "\n" << A.ToMultiLineString() << "\n";
			int k = std::stoi(UI::EnterString("k"));
			if (k > MAX_DIM) {
				UI::Message("Illegal subspace dimension.");
				continue;
			}
			Matrix U(n, k);
			U.Paste(Matrix(k, k) + 1, 0, 0);
			Matrices matrices;
			if (!UI::FileRead(UI::EnterString("File containing matrices E: "), matrices)) {
				UI::Message("Error reading file.");
				continue;
			}
			std::string filename = UI::EnterString("Output file name: ");
			std::ofstream file(filename);
			if (!file.is_open()) {
				UI::Message("Unable to open output file.");
				continue;
			}
			for (Matrix E : matrices) {
				Matrices chain;
				chain.push_back(E * U);
				chain.push_back(A * E * U);
				chain.push_back(A * A * E * U);
				chain.push_back(A * A * A * E * U);
				chain.push_back(E * A * A * A * E * U);
				Matrix M = U; // matrix with independent columns; will grow
				int rank = k; // rank of M
				std::string s = std::to_string(rank);
				for (Matrix T : chain) {
					for (int j = 0; j < k; j++) { // add column by column to M
						Matrix column = T.SubMatrix(0, j, n, j + 1);
						if (!column.IsAX(M)) {
							Matrix newM(n, rank + 1);
							newM.Paste(M, 0, 0);
							newM.Paste(column, 0, rank);
							M = newM;
							rank++;
						}
					}
					s += ",";
					s += std::to_string(rank);
				}
				file << s << "\n";
			}
			file.close();
			UI::Message("Done.");
		}
		else if (command == 5) {
			Matrices matrices;
			if (!UI::FileRead(UI::EnterString("Input file name: "), matrices)) {
				UI::Message("Error reading file.");
				continue;
			}
			if (matrices.empty()) {
				UI::Message("Empty matrix list.");
				continue;
			}
			std::cout << "Matrix dimensions: " << matrices.at(0).Height << " x " << matrices.at(0).Width << "\n";
			std::cout << "Number of matrices: " << matrices.size() << "\n";
			std::sort(
				matrices.begin(),
				matrices.end(),
				[](const Matrix& A, const Matrix& B) -> bool {
					int zerosA = A.CountZeros();
					int zerosB = B.CountZeros();
					if (zerosA > zerosB)
						return true;
					if (zerosA < zerosB)
						return false;
					return A < B;
				});
			std::cout << "Maximum number of zeros: " << matrices.at(0).CountZeros() << "\n";
			if (UI::YesNoQuestion("Write sorted list of matrices? ")) {
				if (UI::FileWrite(UI::EnterString("Output file name: "), matrices))
					UI::Message("Done.");
				else
					UI::Message("Error reading/writing file.");
			}
		}
		else if (command == 4) {
			Matrices matrices;
			if (!UI::FileRead(UI::EnterString("Input file name: "), matrices)) {
				UI::Message("Error reading file.");
				continue;
			}
			Matrix Ones = Matrix(UI::EnterString("Locations of ones: "));
			Matrix Zeros = Matrix(UI::EnterString("Locations of zeros: "));
			Matrices selected;
			for (Matrix A : matrices) {
				bool success = true;
				for (int i = 0; success && i < A.Height; i++)
					for (int j = 0; j < A.Width; j++) {
						if (Ones.GetEntry(i, j) && !A.GetEntry(i, j)) {
							success = false;
							break;
						}
						if (Zeros.GetEntry(i, j) && A.GetEntry(i, j)) {
							success = false;
							break;
						}
					}
				if (success)
					selected.push_back(A);
			}
			std::cout << "Number of matrices found: " << selected.size() << "\n";
			if (UI::YesNoQuestion("Write selected matrices? ")) {
				if (UI::FileWrite(UI::EnterString("Output file name: "), selected))
					UI::Message("Done.");
				else
					UI::Message("Error reading/writing file.");
			}
		}
		else if (command == 3) {
			Matrices matrices;
			if (UI::FileRead(UI::EnterString("Input file name: "), matrices) &&
				UI::FileWrite(UI::EnterString("Output file name: "), matrices, true))
				UI::Message("Done.");
			else
				UI::Message("Error reading/writing file.");
		}
		else if (command == 2) {
			Matrices matrices;
			if (UI::FileRead(UI::EnterString("File name: "), matrices))
				for (int i = 0; i < matrices.size(); i++)
					std::cout << "Matrix no. " << i + 1 << ":\n" << matrices.at(i).ToMultiLineString() << "\n";
			else
				UI::Message("Error opening file.");
		}
		else {
			std::string stringA = UI::EnterString("Coefficients of companion matrices of A (example: m1111,101): ");
			int n = std::count(stringA.begin(), stringA.end(), '0') + std::count(stringA.begin(), stringA.end(), '1');
			if (n == 0 || n > MAX_DIM) {
				UI::Message("Illegal matrix size.");
				continue;
			}
			Matrix A = Companion::BlockCompanionMatrix(stringA);
			std::cout << "\n" << A.ToMultiLineString() << "\n";
			bool find_one = UI::YesNoQuestion("Find one idempotent? ");
			int rand = 0, rand_sets = 0;
			if (UI::YesNoQuestion("Use random heuristic? ")) {
				rand = std::stoi(UI::EnterString("Number of random tries: "));
				if (UI::YesNoQuestion("Random set selection? "))
					rand_sets = std::stoi(UI::EnterString("Number of sets: "));
			}
			Subset M0(A.Height);
			if (rand_sets == 0)
				M0 = Subset(A.Height, UI::EnterString("Begin with the subset: "));

			Matrices idempotents;
			time_t t1 = time(0);
			NilClean::FindIdempotents(A, find_one, M0, idempotents, rand, rand_sets);
			time_t t2 = time(0);

			std::cout << "Finished in " << difftime(t2, t1) << " sec.\n";
			if (idempotents.empty()) {
				UI::Message("No idempotents found.");
				continue;
			}
			else if (idempotents.size() == 1)
				std::cout << "Found idempotent:\n" << idempotents.at(0).ToMultiLineString();
			else
				std::cout << "Number of idempotents found: " << idempotents.size();
			std::cout << "\nChecking ...";
			for (Matrix E : idempotents) {
				Matrix Q = A + E;
				if ((E * E != E) || (Q * Q * Q != 0)) {
					UI::Message(" FAILED!!! Bug in the program!");
					return 0;
				}
			}
			UI::Message(" Done.");
			if (UI::YesNoQuestion("Write idempotents into a file? ")) {
				if (UI::FileWrite(UI::EnterString("File name: "), idempotents))
					UI::Message("Done.");
				else
					UI::Message("Error writing file.");
			}
		}
	}
}