#ifndef POLY_MATRIX_H
#define POLY_MATRIX_H

#include <vector>
#include <cstdint>
#include <string>
#include <utility>
#include "Matrix.h"
#include "Matrix8.h"

// Polynomial over F_2.
class Polynomial {

	// Polynomial coefficients with increasing degree order.
	// Empty if the polynomial is zero.
	std::vector<bool> m_coefficients;

public:

	// Polynomial degree. Degree of 0 is 0.
	std::uint32_t degree() const;

	// Constructors:

	// Zero polynomial.
	Polynomial() = default;

	// Monomial of given degree.
	Polynomial(std::uint32_t degree);

	// Convert string to polynomial.
	// For example, Polynomial("01101") -> x+x^2+x^4,
	// Polynomial("0") -> 0.
	Polynomial(std::string const &inputStr);

	// Conversion to string:

	// Convert polynomial to string.
	// For example, toString(1+x^2) -> "1 + x^2".
	std::string toString() const;

	// Algebraic operations:

	Polynomial operator+(Polynomial const &p) const;

	Polynomial& operator+=(Polynomial const &p);

	Polynomial operator*(Polynomial const &p) const;

	Polynomial& operator*=(Polynomial const &p);

	// Polynomial quotient.
	Polynomial operator/(Polynomial const &p) const;

	// Division remainder.
	Polynomial mod(Polynomial const &p) const;

	// Check if the polynomial is irreducible. Use slow algorithm.
	bool irreducible() const;

	// Relations:

	bool operator==(Polynomial const &p) const;

	// Check if the polynomial is equal to a constant polynomial 0 or 1.
	bool operator==(std::int8_t k) const;

	bool operator!=(Polynomial const &p) const;

	// Check if the polynomial is not equal to a constant polynomial 0 or 1.
	bool operator!=(std::int8_t k) const;

	// Other functions:

	// Replace a polynomial with its successor
	// in the well-ordered set of all polynomials.
	// Polynomials are ordered as binary numbers, for example,
	// 0 < 1 < x^3 < 1 + x^3 < x + x^3 < x^4.
	void inc();

};

// Polynomial matrix of degree at most 1 over
// variables x[i,j] for 0 <= i, j < 8 over the field F_2.
// Coefficients of such a matrix are linear combinations
// x[i1, j1] + ... + x[in, jn] + lambda,
// where lambda is a scalar (0 or 1).
// Entries in the matrix are presented as Matrix8
// (denoting which monomials appear in the entry)
// together with the corresponding entry in another matrix of constants.
class PolyMatrix {

	std::uint8_t m_height;
	
	std::uint8_t m_width;

	// Homogeneous part of matrix entries
	// (homogeneous polynomials of degree 1).
	Matrix8 m_variablePart[MAX_DIM][MAX_DIM];

	// Constant part of matrix entries.
	Matrix m_constantPart;

public:

	std::uint8_t getHeight() const;

	std::uint8_t getWidth() const;

	// Get (i, j)-th matrix entry.
	void getEntry(std::uint8_t i, std::uint8_t j,
		Matrix8 &varPart, bool &constPart) const;

	// Set (i, j)-th matrix entry.
	void setEntry(std::uint8_t i, std::uint8_t j,
		Matrix8 const &varPart, bool constPart);

	// Set a m x n matrix to zero or identity (when m = n).
	PolyMatrix& operator=(std::int8_t value);

	// Constructors:

	// Zero m x n matrix or monomial m x n matrix (x[i, j])_{i, j}
	// (when m, n <= 8).
	PolyMatrix(std::uint8_t m = 0, std::uint8_t n = 0, bool monomial = false);

	// Constant PolyMatrix.
	PolyMatrix(Matrix const &A);

	// Algebraic operations:

	PolyMatrix operator+(PolyMatrix const &A) const;

	// Add a scalar matrix to a n x n matrix.
	PolyMatrix operator+(std::int8_t value) const;

	PolyMatrix& operator+=(PolyMatrix const &A);

	PolyMatrix& operator+=(std::int8_t value);

	// Right multiply by a scalar matrix.
	PolyMatrix operator*(Matrix const &A) const;

	// Relations:

	// Matrix equality. Matrices must be of equal sizes.
	bool operator==(PolyMatrix const &A) const;

	// Matrix inequality. Matrices must be of equal sizes.
	bool operator!=(PolyMatrix const &A) const;

	// Check if a n x n matrix equals a scalar matrix.
	bool operator==(std::int8_t value) const;

	// Check if a n x n matrix equals a scalar matrix.
	bool operator!=(std::int8_t value) const;

	// Linear equation is a pair (Matrix8, bool) where
	// Matrix8 represents all variables in the equation and bool is the
	// constant part of the equation.
	using EquationList = std::vector<std::pair<Matrix8, bool>>;

	// Matrix entries as equation list.
	EquationList equationList() const;

	// In the matrix entries,
	// perform substitutions given in the equations,
	// by replacing the lowest bits with higher bits.
	void substitute(EquationList const &eqnList);

};

// Matrix times a polynomial matrix.
PolyMatrix operator*(Matrix const &A, PolyMatrix const &B);

#endif