#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include "PolyMatrix.h"

std::uint32_t Polynomial::degree() const
{
	return m_coefficients.empty() ? 0 :
		static_cast<std::uint32_t>(m_coefficients.size()) - 1;
}

Polynomial::Polynomial(std::uint32_t degree)
{
	m_coefficients = std::vector<bool>(static_cast<std::size_t>(degree), false);
	m_coefficients.push_back(true);
}

Polynomial::Polynomial(std::string const &inputStr)
{
	if (inputStr == "0") {
		return;
	}
	m_coefficients.reserve(inputStr.length());
	for (char c : inputStr) {
		m_coefficients.push_back(c == '1');
	}
}

std::string Polynomial::toString() const
{
	if (m_coefficients.empty()) {
		return "0";
	}
	std::string result;
	for (std::size_t i = 0; i < m_coefficients.size(); i++) {
		if (m_coefficients[i]) {
			if (result != "")
				result += " + ";
			if (i == 0)
				result += "1";
			else if (i == 1)
				result += "x";
			else {
				result += "x^";
				result += std::to_string(i);
			}
		}
	}
	return result;
}

Polynomial Polynomial::operator+(Polynomial const &p) const
{
	Polynomial result = p;
	for (std::size_t i = 0; i < m_coefficients.size(); i++) {
		if (m_coefficients[i]) {
			if (i < result.m_coefficients.size()) {
				result.m_coefficients[i] = !result.m_coefficients[i];
			}
			else {
				result.m_coefficients.push_back(m_coefficients[i]);
			}
		}
	}
	while (!result.m_coefficients.empty() && !result.m_coefficients.back()) {
		result.m_coefficients.pop_back();
	}
	return result;
}

Polynomial& Polynomial::operator+=(Polynomial const &p)
{
	*this = (*this) + p;
	return *this;
}

Polynomial Polynomial::operator*(Polynomial const &p) const
{
	Polynomial result;
	for (std::size_t i = 0; i < m_coefficients.size(); i++) {
		for (std::size_t j = 0; j < p.m_coefficients.size(); j++) {
			if (m_coefficients[i] && p.m_coefficients[j]) {
				result += Polynomial(static_cast<std::uint32_t>(i + j));
			}
		}
	}
	return result;
}

Polynomial& Polynomial::operator*=(Polynomial const &p)
{
	*this = (*this) * p;
	return *this;
}

Polynomial Polynomial::operator/(Polynomial const &p) const
{
	Polynomial quotient;
	Polynomial remainder = *this;
	while (remainder.m_coefficients.size() >= p.m_coefficients.size()) {
		std::uint32_t n = static_cast<std::uint32_t>(
			remainder.m_coefficients.size() - p.m_coefficients.size());
		quotient = quotient + Polynomial(n);
		remainder = remainder + Polynomial(n) * p;
	}
	return quotient;
}

Polynomial Polynomial::mod(Polynomial const &p) const
{
	return *this + (*this / p) * p;
}

bool Polynomial::irreducible() const
{
	if (m_coefficients.empty()) {
		return false;
	}
	Polynomial divisorCandidate(1);
	while (divisorCandidate.degree() * 2 <= degree()) {
		if (mod(divisorCandidate) == 0) {
			return false;
		}
		divisorCandidate.inc();
	}
	return true;
}

bool Polynomial::operator==(Polynomial const &p) const
{
	if (m_coefficients.size() != p.m_coefficients.size()) {
		return false;
	}
	for (std::size_t i = 0; i < m_coefficients.size(); i++) {
		if (m_coefficients[i] != p.m_coefficients[i]) {
			return false;
		}
	}
	return true;
}

bool Polynomial::operator!=(Polynomial const &p) const
{
	return !(*this == p);
}

bool Polynomial::operator==(std::int8_t k) const
{
	return m_coefficients.size() == k;
}

bool Polynomial::operator!=(std::int8_t k) const
{
	return !(*this == k);
}

void Polynomial::inc()
{
	for (std::size_t i = 0; i < m_coefficients.size(); i++) {
		if (m_coefficients[i]) {
			m_coefficients[i] = false;
		}
		else {
			m_coefficients[i] = true;
			return;
		}
	}
	m_coefficients.push_back(true);
}

std::uint8_t PolyMatrix::getHeight() const
{
	return m_height;
}

std::uint8_t PolyMatrix::getWidth() const
{
	return m_width;
}

void PolyMatrix::getEntry(std::uint8_t i, std::uint8_t j,
	Matrix8 &varPart, bool &constPart) const
{
	varPart = m_variablePart[i][j];
	constPart = m_constantPart.getEntry(i, j);
}

void PolyMatrix::setEntry(std::uint8_t i, std::uint8_t j,
	Matrix8 const &varPart, bool constPart)
{
	m_variablePart[i][j] = varPart;
	m_constantPart.setEntry(i, j, constPart);
}

PolyMatrix& PolyMatrix::operator=(std::int8_t value)
{
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			m_variablePart[i][j] = 0;
		}
	}
	m_constantPart = value;
	return *this;
}

PolyMatrix::PolyMatrix(std::uint8_t m, std::uint8_t n, bool monomial) :
	m_height{m},
	m_width{n},
	m_constantPart{Matrix(m, n)}
{
	if (!monomial) {
		return;
	}
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			m_variablePart[i][j] = Matrix8(i, j);
		}
	}
}

PolyMatrix::PolyMatrix(Matrix const &A) :
	m_height{A.getHeight()},
	m_width{A.getWidth()},
	m_constantPart{A}
{}

PolyMatrix PolyMatrix::operator+(PolyMatrix const &A) const
{
	PolyMatrix sum = A;
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			sum.m_variablePart[i][j] += m_variablePart[i][j];
		}
	}
	sum.m_constantPart += m_constantPart;
	return sum;
}

PolyMatrix PolyMatrix::operator+(std::int8_t value) const
{
	PolyMatrix sum = *this;
	sum.m_constantPart += value;
	return sum;
}

PolyMatrix& PolyMatrix::operator+=(PolyMatrix const &A)
{
	*this = (*this) + A;
	return *this;
}

PolyMatrix& PolyMatrix::operator+=(std::int8_t value)
{
	*this = (*this) + value;
	return *this;
}

PolyMatrix PolyMatrix::operator*(Matrix const &A) const
{
	PolyMatrix result(m_height, A.getWidth());
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < A.getWidth(); j++) {
			for (std::uint8_t k = 0; k < m_width; k++) {
				if (A.getEntry(k, j)) {
					result.m_variablePart[i][j] += m_variablePart[i][k];
				}
			}
		}
	}
	result.m_constantPart *= A;
	return result;
}

bool PolyMatrix::operator==(PolyMatrix const &A) const
{
	if (m_constantPart != A.m_constantPart) {
		return false;
	}
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			if (m_variablePart[i][j] != A.m_variablePart[i][j]) {
				return false;
			}
		}
	}
	return true;
}

bool PolyMatrix::operator!=(PolyMatrix const &A) const
{
	return !(*this == A);
}

bool PolyMatrix::operator==(std::int8_t value) const
{
	if (m_constantPart != value) {
		return false;
	}
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			if (m_variablePart[i][j] != 0) {
				return false;
			}
		}
	}
	return true;
}

bool PolyMatrix::operator!=(std::int8_t value) const
{
	return !(*this == value);
}

PolyMatrix::EquationList PolyMatrix::equationList() const
{
	EquationList result;
	for (std::uint8_t i = 0; i < m_height; i++) {
		for (std::uint8_t j = 0; j < m_width; j++) {
			result.push_back({m_variablePart[i][j], m_constantPart.getEntry(i, j)});
		}
	}
	return result;
}

void PolyMatrix::substitute(EquationList const &eqnList)
{
	for (auto const &equation : eqnList) {
		std::uint64_t pivot = equation.first.getEntries() &
			-static_cast<std::int64_t>(equation.first.getEntries());
		for (std::uint8_t i = 0; i < m_height; i++) {
			for (std::uint8_t j = 0; j < m_width; j++) {
				if (m_variablePart[i][j].getEntries() & pivot) {
					m_variablePart[i][j] += equation.first;
					m_constantPart.setEntry(i, j,
						m_constantPart.getEntry(i, j) ^ equation.second);
				}
			}
		}
	}
}

PolyMatrix operator*(Matrix const &A, PolyMatrix const &B)
{
	PolyMatrix result(A.getHeight(), B.getWidth());
	for (std::uint8_t i = 0; i < A.getHeight(); i++) {
		for (std::uint8_t j = 0; j < B.getWidth(); j++) {
			for (std::uint8_t k = 0; k < A.getWidth(); k++) {
				if (A.getEntry(i, k)) {
					Matrix8 varPart;
					bool constPart;
					B.getEntry(k, j, varPart, constPart);
					PolyMatrix dResult(A.getHeight(), B.getWidth());
					dResult.setEntry(i, j, varPart, constPart);
					result += dResult;
				}
			}
		}
	}
	return result;
}