#include <cstdint>
#include <string>
#include "IntegerSet.h"
#include "Matrix.h"

std::uint8_t IntegerSet::getN() const
{
	return m_n;
}

std::uint8_t IntegerSet::getSize() const
{
	return m_size;
}

std::uint8_t IntegerSet::getElement(std::uint8_t k) const
{
	return m_elements[k];
}

std::uint32_t IntegerSet::complexity() const
{
	std::uint32_t result = 0;
	for (std::uint8_t i = 0; i < m_size; i++) {
		result += m_elements[i];
	}
	std::uint32_t n = m_size;
	return result - n * (n - 1) / 2;
}

IntegerSet::IntegerSet(std::uint8_t n, std::uint8_t size) :
	m_n{n},
	m_size{size}
{
	for (std::uint8_t k = 0; k < size; k++) {
		m_elements[k] = k;
	}
}

IntegerSet::IntegerSet(std::uint8_t n, std::string const &inputStr):
	m_n{n},
	m_size{0}
{
	std::string s = inputStr;
	if (!s.empty()) {
		s += ",";
	}
	while (!s.empty()) {
		std::size_t pos = s.find(',');
		m_elements[m_size++] = std::stoi(s.substr(0, pos));
		s = s.substr(pos + 1);
	}
}

bool IntegerSet::inc()
{
	for (std::uint8_t i = m_size; i > 0; i--) {
		if (m_elements[i - 1] != m_n - m_size + i - 1) {
			m_elements[i - 1]++;
			for (std::uint8_t j = i; j < m_size; j++) {
				m_elements[j] = m_elements[i - 1] + j - i + 1;
			}
			return true;
		}
	}
	if (m_size != m_n) {
		*this = IntegerSet(m_n, m_size + 1);
		return true;
	}
	return false;
}

std::string IntegerSet::toString() const
{
	std::string s;
	for (std::uint8_t i = 0; i < m_size; i++) {
		s += std::to_string(m_elements[i]);
		if (i != m_size - 1) {
			s += ",";
		}
	}
	return s;
}

Matrix IntegerSet::permutationMatrix() const
{
	Matrix P(m_n, m_n);
	for (std::uint8_t j = 0; j < m_size; j++) {
		P.setEntry(m_elements[j], j, true);
	}
	for (std::uint8_t j = m_size; j < m_n; j++) {
		std::uint8_t i = j - m_size;
		for (std::uint8_t k = 0; k < m_size; k++) {
			if (m_elements[k] <= j - m_size + k) {
				i++;
			}
		}
		P.setEntry(i, j, true);
	}
	return P;
}

bool IntegerSet::contains(std::uint8_t value) const
{
	for (std::uint8_t i = 0; i < m_size; i++) {
		if (m_elements[i] == value) {
			return true;
		}
	}
	return false;
}

IntegerSet IntegerSet::complement() const
{
	IntegerSet A(m_n);
	std::uint8_t newSize = 0;
	for (std::uint8_t i = 0; i < m_n; i++) {
		if (!contains(i)) {
			A.m_elements[newSize++] = i;
		}
	}
	A.m_size = newSize;
	return A;
}