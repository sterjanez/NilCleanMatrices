#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include "UserInterface.h"

void UserInterface::message(std::string const &text)
{
	std::cout << text << "\n";
}

bool UserInterface::yesNoQuestion(std::string const &question)
{
	std::cout << question << " ";
	std::string s;
	do {
		std::getline(std::cin, s);
		if (s != "y" && s != "n")
			std::cout << "Invalid command. Enter y/n: ";
	}
	while (s != "y" && s != "n");
	return s == "y";
}

std::string UserInterface::enterString(std::string const &query)
{
	std::cout << query;
	std::string s;
	std::getline(std::cin, s);
	return s;
}

bool UserInterface::fileRead(std::string const &fileName, std::vector<Matrix> &matrices)
{
	std::ifstream file(fileName);
	if (!file.is_open()) {
		return false;
	}
	std::string line;
	while (std::getline(file, line)) {
		matrices.push_back(line);
	}
	file.close();
	return true;
}

bool UserInterface::fileWrite(std::string const &fileName,
	std::vector<Matrix> const &matrices, bool multi_line)
{
	std::ofstream file(fileName);
	if (!file.is_open()) {
		return false;
	}
	for (Matrix const &A : matrices) {
		file << (multi_line ? A.toMultiLineString() : A.toString()) << "\n";
	}
	file.close();
	return true;
}

std::uint8_t UserInterface::command(std::vector<std::string> const &commands)
{
	std::cout << "\nCommands:\n";
	for (std::size_t i = 0; i < commands.size(); i++) {
		if (i < 9) {
			std::cout << " ";
		}
		std::cout << i + 1 << " " << commands[i] << "\n";
	}
	while (true) {
		std::string command = enterString("Command: ");
		for (std::size_t i = 0; i < commands.size(); i++) {
			if (std::to_string(i + 1) == command) {
				return static_cast<std::uint8_t>(i + 1);
			}
		}
		message("Invalid command.");
	}
}