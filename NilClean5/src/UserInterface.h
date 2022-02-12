#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include <string>
#include <vector>
#include "Matrix.h"

namespace UserInterface {

	// Send message.
	void message(std::string const &text);

	// Prompt for yes/no answer.
	bool yesNoQuestion(std::string const &question);

	// Prompt for a string.
	std::string enterString(std::string const &query);

	// Read a file containing a list of matrices.
	// Return false if the file is not accessible.
	// Append matrices to the list of matrices.
	// Example of input file text format:
	// "10,01\nc1;c11" (returns {[10;01], [1] + [01;11]}).
	bool fileRead(std::string const &fileName, std::vector<Matrix> &matrices);

	// Write matrices into a file. Return false if the file is not accessible.
	bool fileWrite(std::string const &fileName,
		std::vector<Matrix> const &matrices, bool multi_line = false);

	// Prompt for a command. Return the index of the selected command beginning with 1.
	// The number of commands must not exceed 99.
	std::uint8_t command(std::vector<std::string> const &commands);

}

#endif