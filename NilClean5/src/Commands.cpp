#include <cstddef>
#include <cstdint>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "IntegerSet.h"
#include "NilClean.h"
#include "UserInterface.h"
#include "PolyMatrix.h"
#include "Commands.h"

void Commands::findNilCleanForAllMatrices()
{
	std::string stringA = UserInterface::enterString(
		"Structure of matrices A (example: m1xx0;cxx1;1x,x1): ");
	std::ptrdiff_t xcount = std::count(stringA.begin(), stringA.end(), 'x');
	std::ptrdiff_t n = std::count(stringA.begin(), stringA.end(), '0') +
		std::count(stringA.begin(), stringA.end(), '1') + xcount;
	if (n == 0 || n > MAX_DIM) {
		UserInterface::message("Illegal matrix size.");
		return;
	}
	std::string fileBaseName =
		UserInterface::enterString("Base name of the output files: ");
	std::string startingSet =
		UserInterface::enterString("Starting set: ");
	bool findOne =
		UserInterface::yesNoQuestion("Find one idempotent for each matrix (y/n)?");
	bool random =
		UserInterface::yesNoQuestion("Use random heuristic?");
	std::uint64_t randomCount =
		random ? std::stoi(UserInterface::enterString("Number of tries: ")) : 0;
	std::uint32_t indexMax = static_cast<std::uint32_t>(1) << static_cast<std::uint8_t>(xcount);
	for (std::uint32_t index = 0; index < indexMax; index++) {
		std::string stringAmodified;
		std::uint8_t indexX = 0;
		for (char ch : stringA) {
			if (ch == 'x') {
				stringAmodified +=
					index & (static_cast<std::uint32_t>(1) << (indexX++)) ? '1' : '0';
			}
			else {
				stringAmodified += ch;
			}
		}
		Matrix A(stringAmodified);
		std::cout << "\n" << A.toMultiLineString() << "\n";
		IntegerSet M0(A.getHeight(), startingSet);
		std::vector<Matrix> idempotents;
		NilClean::findIdempotents(A, findOne, M0, idempotents, randomCount, 0);
		std::cout << "Number of idempotents found: " << idempotents.size() << "\n";
		std::cout << "Checking ...";
		for (Matrix const &E : idempotents) {
			Matrix Q = A + E;
			if (((E ^ 2) != E) || ((Q ^ 3) != 0)) {
				std::cout << " FAILED!!! Bug in the program!\n";
				return;
			}
		}
		UserInterface::message(" Done.");
		std::string fileName = fileBaseName + std::to_string(index) + ".txt";
		UserInterface::fileWrite(fileName, idempotents);
		std::cout << "Writing file " << fileName << " complete.\n";
	}
}

std::vector<std::vector<Matrix>> enterMatrixLists()
{
	std::vector<std::vector<Matrix>> lists;
	std::string s;
	do {
		s = UserInterface::enterString(
			"File no. " + std::to_string(lists.size() + 1) + " (Enter to end): ");
		if (!s.empty()) {
			std::vector<Matrix> matrices;
			if (UserInterface::fileRead(s, matrices)) {
				lists.push_back(matrices);
			}
			else {
				UserInterface::message("File could not be read.");
			}
		}
	} while (!s.empty());
	return lists;
}

void Commands::findMatricesCommonColumns()
{
	UserInterface::message("Given nonempty lists l1, ..., lk of m x n matrices");
	UserInterface::message("and column indices, find all matrices A1 in l1 such that");
	UserInterface::message("there exist matrices A2, ..., Ak in l2, ..., lk having");
	UserInterface::message("common columns with A1 with given indices.");
	std::vector<std::vector<Matrix>> lists = enterMatrixLists();
	if (lists.empty() || lists[0].empty()) {
		UserInterface::message("Empty list.");
		return;
	}
	std::uint8_t n = lists[0][0].getWidth();
	if (n == 0) {
		UserInterface::message("Empty matrices.");
		return;
	}
	IntegerSet indicesSet(n, UserInterface::enterString("Columns indices: "));
	std::vector<std::uint8_t> indices;
	for (std::uint8_t i = 0; i < indicesSet.getSize(); i++) {
		indices.push_back(indicesSet.getElement(i));
	}
	std::vector<Matrix> solutions = Matrix::commonColumns(lists, indices);
	std::cout << "Number of solutions found: " << solutions.size();
	if (!solutions.empty() && UserInterface::yesNoQuestion("Write solutions?")) {
		if (UserInterface::fileWrite(
			UserInterface::enterString("File name: "), solutions, false))
		{
			UserInterface::message("Done.");
		}
		else {
			UserInterface::message("Error writing file.");
		}
	}
}

void Commands::countMaximalCommonColumns()
{
	UserInterface::message("Given nonempty lists l1, ..., lk of m x n matrices,");
	UserInterface::message("find maximal integer k such that there exist matrices");
	UserInterface::message("A1, ..., Ak in l1, ..., lk having k common columns.");
	std::vector<std::vector<Matrix>> lists = enterMatrixLists();
	if (lists.empty()) {
		UserInterface::message("Empty list.");
		return;
	}
	for (auto const &list : lists) {
		if (list.empty()) {
			UserInterface::message("Empty list.");
			return;
		}
	}
	std::uint8_t n = lists[0][0].getWidth();
	IntegerSet MissingIndices(n);
	do {
		IntegerSet indicesSet = MissingIndices.complement();
		std::vector<std::uint8_t> indices;
		for (std::uint8_t i = 0; i < indicesSet.getSize(); i++) {
			indices.push_back(indicesSet.getElement(i));
		}
		std::vector<Matrix> solutions = Matrix::commonColumns(lists, indices);
		if (!solutions.empty()) {
			std::cout << "Found matrix:\n" << solutions[0].toMultiLineString() << "\n";
			std::cout << "Indices of common columns: " << indicesSet.toString() << "\n";
			return;
		}
	}
	while (MissingIndices.inc() && MissingIndices.getSize() != n);
	UserInterface::message("No common columns found.");
}

void Commands::countMaximalCommonEntries()
{
	UserInterface::message("Given nonempty lists l1, ..., lk of m x n matrices,");
	UserInterface::message("find matrices Ai in li such that the number of");
	UserInterface::message("common entries in A1, ..., Ak is maximal possible.");
	std::vector<std::vector<Matrix>> lists = enterMatrixLists();
	if (lists.empty()) {
		UserInterface::message("Empty list.");
		return;
	}
	std::uint16_t commonCount;
	std::vector<std::size_t> indices = Matrix::maxCommonEntries(lists, commonCount);
	UserInterface::message("Found indices: ");
	for (std::size_t k = 0; k < indices.size(); k++) {
		if (k != 0) {
			std::cout << ",";
		}
		std::cout << indices[k] + 1;
	}
	std::cout << "\nNumber of common entries: " << commonCount << "\n";
	if (UserInterface::yesNoQuestion("Write matrices into a file?")) {
		std::vector<Matrix> matrices;
		for (std::size_t k = 0; k < indices.size(); k++) {
			matrices.push_back(lists[k][indices[k]]);
		}
		if (UserInterface::fileWrite(UserInterface::enterString("File name: "), matrices)) {
			UserInterface::message("Done.");
		}
		else {
			UserInterface::message("Error writing file.");
		}
	}
}

void Commands::bitwiseAndOrMatrices()
{
	std::vector<Matrix> matrices;
	if (!UserInterface::fileRead(
		UserInterface::enterString("Input file name: "), matrices))
	{
		UserInterface::message("Error reading file.");
		return;
	}
	if (matrices.empty()) {
		UserInterface::message("Empty matrix list.");
		return;
	}
	bool operation = UserInterface::yesNoQuestion("Bitwise operation (Yes=AND, No=OR):");
	Matrix A = matrices[0];
	for (Matrix const &X : matrices) {
		for (std::uint8_t i = 0; i < A.getHeight(); i++) {
			A.setRow(i, operation ? A.getRow(i) & X.getRow(i) : A.getRow(i) | X.getRow(i));
		}
	}
	std::cout << "Output matrix:\n" << A.toMultiLineString() << "\n";
}

void Commands::findIrreduciblePolynomials()
{
	std::uint32_t degree = std::stoi(UserInterface::enterString("Polynomial degree: "));
	UserInterface::message("Irreducible polynomials:");
	Polynomial candidate = Polynomial(degree);
	while (candidate.degree() == degree) {
		if (candidate.irreducible()) {
			std::cout << candidate.toString() << "\n";
		}
		candidate.inc();
	}
	UserInterface::message("Finished.");
}

void Commands::computeChainDimensions()
{
	UserInterface::message("Given matrix A, a set of matrices E, and U = Lin(e1, ..., ek),");
	UserInterface::message("compute dimensions of U, E(U), AE(U), A^2E(U), A^3E(U), EA^3E(U).\n");
	std::string stringA = UserInterface::enterString(
		"Coefficients of companion matrix blocks of A (example: m1111;c101): ");
	std::ptrdiff_t nLong =
		std::count(stringA.begin(), stringA.end(), '0') +
		std::count(stringA.begin(), stringA.end(), '1');
	if (nLong == 0 || nLong > MAX_DIM) {
		UserInterface::message("Illegal matrix size.");
		return;
	}
	std::uint8_t n = static_cast<std::uint8_t>(nLong);
	Matrix A(stringA);
	std::cout << "\n" << A.toMultiLineString() << "\n";
	auto k = std::stoi(UserInterface::enterString("k: "));
	if (k > MAX_DIM) {
		UserInterface::message("Illegal subspace dimension.");
		return;
	}
	Matrix U(n, k);
	U.paste(Matrix(k, k) + 1, 0, 0);
	std::vector<Matrix> matrices;
	if (!UserInterface::fileRead(
		UserInterface::enterString("File containing matrices E: "), matrices))
	{
		UserInterface::message("Error reading file.");
		return;
	}
	std::string fileName = UserInterface::enterString("Output file name: ");
	std::ofstream file(fileName);
	if (!file.is_open()) {
		UserInterface::message("Unable to open output file.");
		return;
	}
	for (Matrix const &E : matrices) {
		std::vector<Matrix> chain;
		chain.push_back(E * U);
		chain.push_back(A * E * U);
		chain.push_back(A * A * E * U);
		chain.push_back(A * A * A * E * U);
		chain.push_back(E * A * A * A * E * U);
		Matrix M = U; // matrix with independent columns; will grow
		std::uint8_t rank = k; // rank of M
		std::string s = std::to_string(rank);
		for (Matrix const &T : chain) {
			for (std::uint8_t j = 0; j < k; j++) { // add column by column to M
				Matrix column = T.submatrix(0, j, n, j + 1);
				if (!column.isAX(M)) {
					Matrix newM(n, rank + 1);
					newM.paste(M, 0, 0);
					newM.paste(column, 0, rank);
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
	UserInterface::message("Done.");
}

void Commands::sortByZeroCount()
{
	std::vector<Matrix> matrices;
	if (!UserInterface::fileRead(UserInterface::enterString("Input file name: "), matrices)) {
		UserInterface::message("Error reading file.");
		return;
	}
	if (matrices.empty()) {
		UserInterface::message("Empty matrix list.");
		return;
	}
	std::cout << "Matrix dimensions: " <<
		matrices[0].getHeight() << " x " << matrices[0].getWidth() << "\n";
	std::cout << "Number of matrices: " << matrices.size() << "\n";
	std::sort(
		matrices.begin(),
		matrices.end(),
		[](Matrix const &A, Matrix const &B) -> bool {
			std::uint16_t zerosA = A.countZeros();
			std::uint16_t zerosB = B.countZeros();
			if (zerosA > zerosB)
				return true;
			if (zerosA < zerosB)
				return false;
			return A < B;
		});
	std::cout << "Maximum number of zeros: " << matrices[0].countZeros() << "\n";
	if (UserInterface::yesNoQuestion("Write sorted list of matrices?")) {
		if (UserInterface::fileWrite(
			UserInterface::enterString("Output file name: "), matrices))
		{
			UserInterface::message("Done.");
		}
		else {
			UserInterface::message("Error reading/writing file.");
		}
	}
}

void Commands::selectMatrices()
{
	std::vector<Matrix> matrices;
	if (!UserInterface::fileRead(UserInterface::enterString("Input file name: "), matrices)) {
		UserInterface::message("Error reading file.");
		return;
	}
	Matrix Ones = Matrix(UserInterface::enterString("Locations of ones: "));
	Matrix Zeros = Matrix(UserInterface::enterString("Locations of zeros: "));
	std::vector<Matrix> selected;
	for (Matrix const &A : matrices) {
		bool success = true;
		for (std::uint8_t i = 0; i < A.getHeight(); i++) {
			if (((A.getRow(i) & Ones.getRow(i)) != Ones.getRow(i)) ||
				((A.getRow(i) & Zeros.getRow(i)) != 0))
			{
				success = false;
				break;
			}
		}
		if (success) {
			selected.push_back(A);
		}
	}
	std::cout << "Number of matrices found: " << selected.size() << "\n";
	if (UserInterface::yesNoQuestion("Write selected matrices?")) {
		if (UserInterface::fileWrite(
			UserInterface::enterString("Output file name: "), selected))
		{
			UserInterface::message("Done.");
		}
		else {
			UserInterface::message("Error writing file.");
		}
	}
}

void Commands::convertToMultiline()
{
	std::vector<Matrix> matrices;
	if (UserInterface::fileRead(UserInterface::enterString("Input file name: "), matrices) &&
		UserInterface::fileWrite(UserInterface::enterString("Output file name: "), matrices, true))
	{
		UserInterface::message("Done.");
	}
	else {
		UserInterface::message("Error reading/writing file.");
	}
}

void Commands::displayMatrices()
{
	std::vector<Matrix> matrices;
	if (UserInterface::fileRead(UserInterface::enterString("File name: "), matrices)) {
		for (std::size_t i = 0; i < matrices.size(); i++) {
			std::cout << "Matrix no. " << i + 1 << ":\n" <<
				matrices[i].toMultiLineString() << "\n";
		}
	}
	else {
		UserInterface::message("Error opening file.");
	}
}

void Commands::findNilClean()
{
	std::string stringA = UserInterface::enterString(
		"Matrix A (example: m111;c101;11,10): ");
	std::ptrdiff_t n =
		std::count(stringA.begin(), stringA.end(), '0') +
		std::count(stringA.begin(), stringA.end(), '1');
	if (n == 0 || n > MAX_DIM) {
		UserInterface::message("Illegal matrix size.");
		return;
	}
	Matrix A(stringA);
	std::cout << "\n" << A.toMultiLineString() << "\n";
	bool findOne = UserInterface::yesNoQuestion("Find one idempotent?");
	bool parallelize = UserInterface::yesNoQuestion("Parallelize?");
	std::uint64_t rand = 0;
	std::uint64_t randSets = 0;
	if (!parallelize && UserInterface::yesNoQuestion("Use random heuristic?")) {
		rand = std::stoi(UserInterface::enterString("Number of random tries: "));
		if (UserInterface::yesNoQuestion("Random set selection?")) {
			randSets = std::stoi(UserInterface::enterString("Number of sets: "));
		}
	}
	IntegerSet M0(A.getHeight());
	if (!parallelize && randSets == 0) {
		M0 = IntegerSet(A.getHeight(), UserInterface::enterString("Begin with the subset: "));
	}
	std::vector<Matrix> idempotents;
	std::time_t t1 = std::time(0);
	if (parallelize) {
		NilClean::findIdempotentsParallel(A, findOne, idempotents);
	}
	else {
		NilClean::findIdempotents(A, findOne, M0, idempotents, rand, randSets);
	}
	std::time_t t2 = std::time(0);
	std::cout << "Finished in " << std::difftime(t2, t1) << " sec.\n";
	if (idempotents.empty()) {
		UserInterface::message("No idempotents found.");
		return;
	}
	if (idempotents.size() == 1) {
		std::cout << "Found idempotent:\n" << idempotents[0].toMultiLineString();
	}
	else {
		std::cout << "Number of idempotents found: " << idempotents.size();
	}
	std::cout << "\nChecking ...";
	for (Matrix const &E : idempotents) {
		Matrix Q = A + E;
		if ((E * E != E) || (Q * Q * Q != 0)) {
			UserInterface::message(" FAILED!!! Bug in the program!");
			return;
		}
	}
	UserInterface::message(" Done.");
	if (UserInterface::yesNoQuestion("Write idempotents into a file?")) {
		if (UserInterface::fileWrite(UserInterface::enterString("File name: "), idempotents)) {
			UserInterface::message("Done.");
		}
		else {
			UserInterface::message("Error writing file.");
		}
	}
}