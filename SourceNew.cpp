// Program: Nil-clean matrices over F_2
// Author: Janez Ster
// Last update: December 28th 2020

// Description: Given a n×n matrix A over F_2 (for small n), the program finds
// idempotent matrices E over F_2 such that (A-E)^3=0.

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>

const int MAX_DIM = 13; // maximal matrix size

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
			C ^= ((A >> i)& FirstCol)* ((B >> N[i])& FirstRow);
		return C;
	}

	// Check if AX = B for some X.
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

	// Check if XA = B for some X.
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

};

typedef std::vector<Matrix> Matrices;

// Given n×n matrix A and 0<p<n, find idempotents E of rank p in block form
// E=[I0;XI][IY;00][I0;XI] (with the upper left block of size p×p),
// such that (A-E)^3=0. Heights of block X are limited by heights[].
// (E.g. heights[] = {0,0,...,0} means that X = 0 is the only option.)
// If find_one is set, find only one idempotent. Append idempotents to the vector idempotents.
// Use matrices presented by single 64-bit integers. Only for use if p,n-p<=8.
void FindIdempotents64(const Matrix& A, int p, int heights[MAX_DIM], bool find_one, Matrices& idempotents) {
	int n = A.Height;
	int q = n - p;

	Matrix B = (A + 1) * (A + 1) * (A + 1);
	Matrix C = A * A * A;

	// The outer loop will run through all matrices P=[I0;XI] with X having given heights.
	Matrix64 X = 0;

	// We are looking for Y such that (PAP + [IY;00])^3 = 0.

	// CounterX determines which entry of X is about to change
	Matrix64 CounterX = 0;

	// Blocks of matrices PAP, PBP, PCP. Will be maintained.
	Matrix64 A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
	Mat64::Decompose(A, p, p, A11, A12, A21, A22);
	Mat64::Decompose(B, p, p, B11, B12, B21, B22);
	Mat64::Decompose(C, p, p, C11, C12, C21, C22);

	bool counterXfound = true;
	while (counterXfound) {
		if (Mat64::AXisB(A21, B21) && Mat64::XAisB(A21, C21)) {
			int rank;
			Matrix64 U, V;
			Mat64::Diagonalize(A21, q, p, rank, U, V);

			if (Mat64::AXisB(Mat64::Product(U, C21), Mat64::Product(U, C22), rank) &&
				Mat64::XAisB(Mat64::Product(B21, V), Mat64::Product(B11, V), rank)) {
				Matrix64 Uinv, Vinv;
				Uinv = Mat64::Inv(U, q);
				Vinv = Mat64::Inv(V, p);

				// We maintain blocks of Q, R=Q^2 and S=Q^3.
				Matrix64 Q11, Q12, Q21, Q22, R11, R12, R21, R22, S11, S12, S21, S22;

				Q11 = Mat64::Product(Mat64::Product(Vinv, A11), V) ^ Mat64::Identity(p);
				Q12 = Mat64::Product(Mat64::Product(Vinv, A12), Uinv) ^ Mat64::Product(Mat64::Product(U, B21 ^ Mat64::Product(A21, A11)), V);
				Q21 = Mat64::Identity(rank);
				Q22 = Mat64::Product(Mat64::Product(U, A22), Uinv);

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

								S11 ^= ((Q11 >> i) & Mat64::FirstCol) * ((Q21 >> Mat64::N[j]) & Mat64::FirstRow);
								S12 ^= ((Q11 >> i) & Mat64::FirstCol) * ((Q22 >> Mat64::N[j]) & Mat64::FirstRow);
								S21 ^= ((Q21 >> i) & Mat64::FirstCol) * ((Q21 >> Mat64::N[j]) & Mat64::FirstRow);
								S22 ^= ((Q21 >> i) & Mat64::FirstCol) * ((Q22 >> Mat64::N[j]) & Mat64::FirstRow);

								S11 ^= ((R21 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
								S12 ^= ((R22 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];

								R12 ^= ((Q11 >> i) & Mat64::FirstCol) << j;
								R22 ^= ((Q21 >> i) & Mat64::FirstCol) << j;

								R11 ^= ((Q21 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
								R12 ^= ((Q22 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];

								S12 ^= ((R11 >> i) & Mat64::FirstCol) << j;
								S22 ^= ((R21 >> i) & Mat64::FirstCol) << j;

								Q12 ^= selector;

								counterYfound = true;
								break;
							}
						}
					}
				}
			}
		}
		counterXfound = false;
		for (int j = 0; j < p && !counterXfound; j++)
			for (int i = 0; i < heights[j]; i++)
				if (CounterX & (1ull << (Mat64::N[i] + j)))
					CounterX ^= 1ull << (Mat64::N[i] + j);
				else {
					CounterX ^= 1ull << (Mat64::N[i] + j);
					X ^= 1ull << (Mat64::N[i] + j);

					A21 ^= ((A11 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					A22 ^= ((A12 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					A11 ^= ((A12 >> i) & Mat64::FirstCol) << j;
					A21 ^= ((A22 >> i) & Mat64::FirstCol) << j;

					B21 ^= ((B11 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					B22 ^= ((B12 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					B11 ^= ((B12 >> i) & Mat64::FirstCol) << j;
					B21 ^= ((B22 >> i) & Mat64::FirstCol) << j;

					C21 ^= ((C11 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					C22 ^= ((C12 >> Mat64::N[j]) & Mat64::FirstRow) << Mat64::N[i];
					C11 ^= ((C12 >> i) & Mat64::FirstCol) << j;
					C21 ^= ((C22 >> i) & Mat64::FirstCol) << j;

					counterXfound = true;
					break;
				}
	}
}

// Given n×n matrix A and 0<p<n, find idempotents E of rank p in block form
// E=[I0;XI][IY;00][I0;XI] (with the upper left block of size p×p),
// such that (A-E)^3=0. Heights of block X are limited by heights[].
// (E.g. heights[] = {0,0,...,0} means that X = 0 is the only option.)
// If find_one is set, find only one idempotent. Append idempotents to the vector idempotents.
void FindIdempotents(const Matrix& A, int p, int heights[MAX_DIM], bool find_one, Matrices& idempotents) {
	int n = A.Height;
	int q = n - p;

	Matrix B = (A + 1) * (A + 1) * (A + 1);
	Matrix C = A * A * A;

	// The outer loop will run through all matrices P=[I0;XI] with X having given heights.
	Matrix P = Matrix(n, n) + 1;

	// We are looking for Y such that (PAP + [IY;00])^3 = 0.

	// CounterX determines which entry of X is about to change
	Matrix CounterX(n - p, p);

	// Blocks of matrices PAP, PBP, PCP. Will be maintained.
	Matrix A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
	A.Decompose(p, p, A11, A12, A21, A22);
	B.Decompose(p, p, B11, B12, B21, B22);
	C.Decompose(p, p, C11, C12, C21, C22);

	bool counterXfound = true;
	while (counterXfound) {
		if (B21.IsAX(A21) && C21.IsXA(A21)) {
			int rank;
			Matrix U, V;
			A21.Diagonalize(rank, U, V);
			if ((U * C22).IsAX(U * C21, rank) && (B11 * V).IsXA(B21 * V, rank)) {
				Matrix Q(n, n);
				Q.Paste(V.Inv() * A11 * V + 1, 0, 0);
				Q.Paste((V.Inv() * A12 * U.Inv()).SubMatrix(0, 0, rank, rank) + (U * (B21 + A21 * A11) * V).SubMatrix(0, 0, rank, rank), 0, p);
				Q.Paste(Matrix(rank, rank) + 1, p, 0);
				Q.Paste(U * A22 * U.Inv(), p, p);

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
		counterXfound = false;
		for (int j = 0; j < p && !counterXfound; j++)
			for (int i = 0; i < heights[j]; i++)
				if (CounterX.GetEntry(i, j))
					CounterX.SetEntry(i, j, false);
				else {
					CounterX.SetEntry(i, j, true);
					P.Row[i + p] ^= 1 << j;

					A21.AddRow(j, A11, i);
					A22.AddRow(j, A12, i);
					A11.AddCol(i, A12, j);
					A21.AddCol(i, A22, j);

					B21.AddRow(j, B11, i);
					B22.AddRow(j, B12, i);
					B11.AddCol(i, B12, j);
					B21.AddCol(i, B22, j);

					C21.AddRow(j, C11, i);
					C22.AddRow(j, C12, i);
					C11.AddCol(i, C12, j);
					C21.AddCol(i, C22, j);

					counterXfound = true;
					break;
				}
	}
}

// Given n×n matrix A and a subset M of {0,...,n-1}, find idempotents E
// such that (A-E)^3=0 and such that P^T*E*P is of the form [I0;XI][IY;00][I0;XI],
// where P is a permuation matrix associated with M. If find_one is set,
// it finds only one idempotent. Append idempotents to the vector idempotents.
void FindIdempotents(const Matrix& A, const Subset& M, bool find_one, Matrices& idempotents) {
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
		FindIdempotents64(PTAP, p, heights, find_one, new_idempotents);
	else
		FindIdempotents(PTAP, M.Order, heights, find_one, new_idempotents);
	for (Matrix E : new_idempotents)
		idempotents.push_back(P * E * P.Transpose());
}

// Given n×n matrix A, find idempotents E such that (A-E)^3=0.
// Write steps in the procedure. Begin with the subset M0.
// If find_one is set, then find just one E.
// Append idempotents to the vector idempotents.
void FindIdempotents(const Matrix& A, bool find_one, const Subset& M0, Matrices& idempotents) {
	Subset M = M0;
	do {
		std::cout << "Set = {" << M.ToString() << "}\n";
		FindIdempotents(A, M, find_one, idempotents);
		if (find_one && !idempotents.empty()) return;
	} while (M.Inc());
}

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

// Send message.
void Message(const std::string& text) {
	std::cout << text << "\n";
}

// Prompt for yes/no answer.
bool YesNoQuestion(const std::string& question) {
	std::cout << question << " (y/n) ";
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
	std::cout << query << ": ";
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

int main() {
	std::cout << "Nil-clean Matrices Over F_2\n";
	std::cout << "Author: Janez Ster\n";
	while (true) {
		std::cout << "\n";
		std::cout << "1 Find idempotent E such that (A-E)^3=0\n";
		std::cout << "2 Find and save all idempotents E such that (A-E)^3=0\n";
		std::cout << "3 Display matrices from a file\n";
		std::cout << "4 Convert matrices in a file to a user-friendly multi-line form\n";
		std::cout << "5 Exit\n";
		std::cout << "\n";
		std::string command = EnterString("Command");
		if (command == "5")
			return 0;
		else if (command == "4") {
			Matrices matrices;
			if (FileRead(EnterString("Input file name"), matrices) &&
				FileWrite(EnterString("Output file name"), matrices, true))
				Message("Done.");
			else
				Message("Error reading/writing file.");
		}
		else if (command == "3") {
			Matrices matrices;
			if (FileRead(EnterString("File name"), matrices))
				for (int i = 0; i < matrices.size(); i++)
					std::cout << "Matrix no. " << i + 1 << ":\n" << matrices.at(i).ToMultiLineString() << "\n";
			else
				Message("Error opening file.");
		}
		else if (command == "2" || command == "1") {
			std::string stringA = EnterString("Coefficients of companion matrices of A (example: m1111,101)");
			int n = std::count(stringA.begin(), stringA.end(), '0') + std::count(stringA.begin(), stringA.end(), '1');
			if (n == 0 || n > MAX_DIM) {
				Message("Illegal matrix size.");
				continue;
			}
			Matrix A = BlockCompanionMatrix(stringA);
			std::cout << "\n" << A.ToMultiLineString() << "\n";
			Subset M0(A.Height, EnterString("Begin with the subset"));
			Matrices idempotents;
			time_t begin = time(0);
			FindIdempotents(A, command == "1", M0, idempotents);
			time_t end = time(0);
			std::cout << "Finished in " << difftime(end, begin) << " sec.\n";
			if (idempotents.empty()) {
				Message("No idempotents found.");
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
					std::cout << " FAILED!!! Bug in the program!\n";
					return 0;
				}
			}
			Message(" Done.");
			if (idempotents.size() == 1) continue;
			if (YesNoQuestion("Write idempotents into a file?")) {
				if (FileWrite(EnterString("File name"), idempotents))
					Message("Done.");
				else
					Message("Error writing file.");
			}
		}
		else
			Message("Invalid command.");
	}
}