#ifndef INTEGER_SET_H
#define INTEGER_SET_H

#include <cstdint>
#include "Matrix.h"

// A subset of {0, ..., n - 1} for 1 <= n <= MAX_DIM.
// Sets are ordered by size,
// and same-sized sets are ordered lexicographically.
// For example, {6} < {1,4} < {1,5} < {2,3}.
class IntegerSet {

	std::uint8_t m_n;

	// Set size.
	std::uint8_t m_size;

	// Set elements.
	std::uint8_t m_elements[MAX_DIM];

public:

	std::uint8_t getN() const;

	std::uint8_t getSize() const;

	std::uint8_t getElement(std::uint8_t k) const;

	// The sum a1 + (a2 - 1) + (a3 - 2) + ... + (an - (n - 1)),
	// where ai are elements of the set.
	std::uint32_t complexity() const;

	// Constructors:

	// The minimal subset of given size.
	IntegerSet(std::uint8_t n, std::uint8_t size = 0);

	// Create a set from a string.
	// Example: IntegerSet("0,3,11") -> {0, 3, 11}.
	IntegerSet(std::uint8_t n, std::string const &inputStr);

	// Increment set by replacing it with its successor.
	// Return false if successor does not exist.
	bool inc();

	// Conversion to string:

	// Convert a set to a string.
	// For example, toStrint({0, 3, 11}) -> "0,3,11".
	std::string toString() const;

	// Permutation n x n matrix associated with a subset of {0, ..., n - 1} of size k.
	// Heights of first k columns are given by the subset in increasing order,
	// and the remaining columns form an increasing sequence as well.
	// Example: subset {1,3} with n=4 is associated with P=[0010;1000;0001;0100]
	// (1st column has height 2, second height 4,
	// and the remaining two columns have heights 1 and 3).
	Matrix permutationMatrix() const;

	// Set operations:

	// Check if the set contains an element.
	bool contains(std::uint8_t value) const;

	// Complement in the universal set {0, ..., N - 1}.
	IntegerSet complement() const;

};

#endif