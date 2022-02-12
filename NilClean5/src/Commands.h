#ifndef COMMANDS_H
#define COMMANDS_H

// A collection of methods that prompt for a user input,
// execute operations, and output the result.
namespace Commands {

	// Find idempotents E such that (A - E) ^ 3 = 0 for a given matrix A.
	void findNilClean();

	// Display matrices from a file.
	void displayMatrices();

	// Convert matrices in a file to multi-line form.
	void convertToMultiline();

	// Select matrices from a file with given entries.
	void selectMatrices();

	// Sort matrices in a file by number of zero entries.
	void sortByZeroCount();

	// Given matrices A, E, and a subspace U,
	// compute dimensions of U, EU, AEU, A^2EU, A^3EU, EA^3EU.
	void computeChainDimensions();

	// Find irreducible polynomials of given dimension.
	void findIrreduciblePolynomials();

	// Bit-wise "and" or "or" operation on a list of matrices in a file.
	void bitwiseAndOrMatrices();

	// Find matrices in lists having maximal number of common entries.
	void countMaximalCommonEntries();

	// Find matrices in lists having maximal number of common columns.
	void countMaximalCommonColumns();

	// Find matrices in lists with common columns with given indices.
	void findMatricesCommonColumns();

	// Find nil-clean decompositions of index at most 3 for a collection of matrices.
	void findNilCleanForAllMatrices();

}

#endif