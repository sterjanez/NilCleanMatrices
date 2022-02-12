// Program: Nil-clean matrices over F_2
// Author: Janez Ster
// Last update: January 2022

// Description: Given a n x n matrix A over F_2 (for small n), the program finds
// idempotent matrices E over F_2 such that (A - E) ^ 3 = 0.

#include <vector>
#include <string>
#include <functional>
#include <cstdint>
#include "UserInterface.h"
#include "Commands.h"

int main() {
	std::vector<std::string> commandStrings = {
		"Find and save idempotents E such that (A - E) ^ 3 = 0",
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
		"Find nil-clean decompositions for a collection of matrices",
		"Exit"
	};
	std::vector<std::function<void(void)>> commandFunctions = {
		Commands::findNilClean,
		Commands::displayMatrices,
		Commands::convertToMultiline,
		Commands::selectMatrices,
		Commands::sortByZeroCount,
		Commands::computeChainDimensions,
		Commands::findIrreduciblePolynomials,
		Commands::bitwiseAndOrMatrices,
		Commands::countMaximalCommonEntries,
		Commands::countMaximalCommonColumns,
		Commands::findMatricesCommonColumns,
		Commands::findNilCleanForAllMatrices
	};
	UserInterface::message("Nil-clean matrices over F_2");
	UserInterface::message("Author: Janez Ster");
	while (true) {
		std::uint8_t command = UserInterface::command(commandStrings);
		if (command == commandStrings.size()) {
			return 0;
		}
		commandFunctions[command - 1]();
	}
}