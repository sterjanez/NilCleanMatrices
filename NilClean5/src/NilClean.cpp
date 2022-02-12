#include "NilClean.h"
#include "Matrix.h"
#include "Matrix8.h"
#include "MatrixInt.h"
#include "IntegerSet.h"
#include <cstdint>
#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <thread>
#include <mutex>
#include <functional>
#include <algorithm>

void NilClean::findIdempotents8(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand)
{
	std::uint8_t n = A.getHeight();
	std::uint8_t q = n - p;

	Matrix B = (A + 1) ^ 3;
	Matrix C = A ^ 3;

	// Allowed nonzero entries of X. Only in use if rand > 0.
	Matrix8 SelectorX;
	for (std::uint8_t j = 0; j < p; j++) {
		for (std::uint8_t i = 0; i < heights[j]; i++) {
			SelectorX += Matrix8(i, j);
		}
	}

	// Random seed, generator and distribution.
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<std::uint64_t> distribution(0, 0xFFFFFFFFFFFFFFFF);

	// The outer loop will run through all matrices P = [I0;XI] with X having given heights.
	// Initializing X in use only if rand = 0.
	Matrix8 X;

	// We are looking for Y such that (PAP + [IY;00]) ^ 3 = 0.

	// CounterX determines which entry of X is about to change. Only in use if rand = 0.
	Matrix8 CounterX;

	// Blocks of matrices A, B, C.
	Matrix8 A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
	A.decompose(p, p, A11, A12, A21, A22);
	B.decompose(p, p, B11, B12, B21, B22);
	C.decompose(p, p, C11, C12, C21, C22);

	// Blocks of matrices PAP, PBP, PCP. Maintained if rand = 0.
	// Initialization relevant only if rand = 0.
	Matrix8 PAP11, PAP12, PAP21, PAP22, PBP11, PBP12, PBP21, PBP22, PCP11, PCP12, PCP21, PCP22;
	A.decompose(p, p, PAP11, PAP12, PAP21, PAP22);
	B.decompose(p, p, PBP11, PBP12, PBP21, PBP22);
	C.decompose(p, p, PCP11, PCP12, PCP21, PCP22);

	// Counter of random attempts completed, irrelevant if rand = 0.
	std::uint64_t attemptsX = 0;

	bool counterXfound = true;
	while (counterXfound) {
		if (rand) {
			// Random choice of X.
			X.setEntries(SelectorX.getEntries() & distribution(generator));

			// Compute blocks of PAP, PBP, PCP.
			PAP11 = A11 + A12 * X;
			PAP21 = A21 + X * PAP11 + A22 * X;
			PAP22 = A22 + X * A12;

			PBP11 = B11 + B12 * X;
			PBP21 = B21 + X * PBP11 + B22 * X;
			PBP22 = B22 + X * B12;

			PCP11 = C11 + C12 * X;
			PCP21 = C21 + X * PCP11 + C22 * X;
			PCP22 = C22 + X * C12;
		}
		if (PBP21.isAX(PAP21) && PCP21.isXA(PAP21)) {
			Matrix8 U, V;
			std::uint8_t rank = PAP21.rank(q, p, U, V);
			if ((U * PCP22).isAX(U * PCP21, rank) && (PBP11 * V).isXA(PBP21 * V, rank)) {
				Matrix8 Uinv = U.inv(q);
				Matrix8 Vinv = V.inv(p);

				// We maintain blocks of Q, R = Q ^ 2 and S = Q ^ 3.
				Matrix8 Q11 = Vinv * PAP11 * V + Matrix8(p);
				Matrix8 Q12 = Vinv * PAP12 * Uinv + U * (PBP21 + PAP21 * PAP11) * V;
				Matrix8 Q21(rank);
				Matrix8 Q22 = U * PAP22 * Uinv;

				Matrix8 R11 = Q11 * Q11 + Q12 * Q21;
				Matrix8 R12 = Q11 * Q12 + Q12 * Q22;
				Matrix8 R21 = Q21 * Q11 + Q22 * Q21;
				Matrix8 R22 = Q21 * Q12 + Q22 * Q22;

				Matrix8 S11 = Q11 * R11 + Q12 * R21;
				Matrix8 S12 = Q11 * R12 + Q12 * R22;
				Matrix8 S21 = Q21 * R11 + Q22 * R21;
				Matrix8 S22 = Q21 * R12 + Q22 * R22;

				// CounterY determines which entry of Q is about to change.
				Matrix8 CounterY;

				bool counterYfound = true;
				while (counterYfound) {
					if (S11.getEntries() == 0 && S12.getEntries() == 0 &&
						S21.getEntries() == 0 && S22.getEntries() == 0)
					{
						Matrix Q(n, n, p, p, Q11, Q12, Q21, Q22);
						Matrix W(n, n, p, p, Vinv, 0, 0, U);
						Matrix P(n, n, p, p, p, 0, X, q);
						idempotents.push_back(A + P * W.inv() * Q * W * P);
						if (find_one) {
							return;
						}
					}
					counterYfound = false;
					for (std::uint8_t i = 0, j0 = rank; i < p && !counterYfound; i++) {
						if (i == rank) {
							j0 = 0;
						}
						for (std::uint8_t j = j0; j < q; j++) {
							Matrix8 Selector(i, j);
							if (CounterY.getEntries() & Selector.getEntries()) {
								CounterY += Selector;
							}
							else {
								CounterY += Selector;
								S11 += Q11.rightMultiply(i, j, Q21);
								S12 += Q11.rightMultiply(i, j, Q22);
								S21 += Q21.rightMultiply(i, j, Q21);
								S22 += Q21.rightMultiply(i, j, Q22);

								S11 += R21.leftMultiply(i, j);
								S12 += R22.leftMultiply(i, j);

								R12 += Q11.rightMultiply(i, j);
								R22 += Q21.rightMultiply(i, j);

								R11 += Q21.leftMultiply(i, j);
								R12 += Q22.leftMultiply(i, j);

								S12 += R11.rightMultiply(i, j);
								S22 += R21.rightMultiply(i, j);

								Q12 += Selector;

								counterYfound = true;
								break;
							}
						}
					}
				}
			}
		}
		if (rand) {
			attemptsX++;
			counterXfound = attemptsX < rand;
		}
		else {
			counterXfound = false;
			for (std::uint8_t j = 0; j < p && !counterXfound; j++) {
				for (std::uint8_t i = 0; i < heights[j]; i++) {
					Matrix8 Selector(i, j);
					if (CounterX.getEntries() & Selector.getEntries()) {
						CounterX += Selector;
					}
					else {
						CounterX += Selector;
						X += Selector;

						PAP21 += PAP11.leftMultiply(i, j);
						PAP22 += PAP12.leftMultiply(i, j);
						PAP11 += PAP12.rightMultiply(i, j);
						PAP21 += PAP22.rightMultiply(i, j);

						PBP21 += PBP11.leftMultiply(i, j);
						PBP22 += PBP12.leftMultiply(i, j);
						PBP11 += PBP12.rightMultiply(i, j);
						PBP21 += PBP22.rightMultiply(i, j);

						PCP21 += PCP11.leftMultiply(i, j);
						PCP22 += PCP12.leftMultiply(i, j);
						PCP11 += PCP12.rightMultiply(i, j);
						PCP21 += PCP22.rightMultiply(i, j);

						counterXfound = true;
						break;
					}
				}
			}
		}
	}
}

void NilClean::findIdempotentsInt(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand)
{
	std::uint8_t n = A.getHeight();
	std::uint8_t q = n - p;

	Matrix B = (A + 1) ^ 3;
	Matrix C = A ^ 3;

	// Allowed nonzero entries of X. Only in use if rand > 0.
	MatrixInt SelectorX = 0;
	for (std::uint8_t j = 0; j < p; j++) {
		for (std::uint8_t i = 0; i < heights[j]; i++) {
			SelectorX ^= MI::create(i, j);
		}
	}

	// Random seed, generator and distribution.
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<std::uint64_t> distribution(0, 0xFFFFFFFFFFFFFFFF);

	// The outer loop will run through all matrices P = [I0;XI] with X having given heights.
	// Initializing X in use only if rand = 0.
	MatrixInt X = 0;

	// We are looking for Y such that (PAP + [IY;00]) ^ 3 = 0.

	// CounterX determines which entry of X is about to change. Only in use if rand = 0.
	MatrixInt CounterX = 0;

	// Blocks of matrices A, B, C.
	MatrixInt A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
	A.decompose(p, p, A11, A12, A21, A22);
	B.decompose(p, p, B11, B12, B21, B22);
	C.decompose(p, p, C11, C12, C21, C22);

	// Blocks of matrices PAP, PBP, PCP. Maintained if rand = 0.
	// Initialization relevant only if rand = 0.
	MatrixInt PAP11, PAP12, PAP21, PAP22, PBP11, PBP12, PBP21, PBP22, PCP11, PCP12, PCP21, PCP22;
	A.decompose(p, p, PAP11, PAP12, PAP21, PAP22);
	B.decompose(p, p, PBP11, PBP12, PBP21, PBP22);
	C.decompose(p, p, PCP11, PCP12, PCP21, PCP22);

	// Counter of random attempts completed, irrelevant if rand = 0.
	std::uint64_t attemptsX = 0;

	bool counterXfound = true;
	while (counterXfound) {
		if (rand) {
			// Random choice of X.
			X = SelectorX & distribution(generator);

			// Compute blocks of PAP, PBP, PCP.
			PAP11 = A11 ^ MI::product(A12, X);
			PAP21 = A21 ^ MI::product(X, PAP11) ^ MI::product(A22, X);
			PAP22 = A22 ^ MI::product(X, A12);

			PBP11 = B11 ^ MI::product(B12, X);
			PBP21 = B21 ^ MI::product(X, PBP11) ^ MI::product(B22, X);
			PBP22 = B22 ^ MI::product(X, B12);

			PCP11 = C11 ^ MI::product(C12, X);
			PCP21 = C21 ^ MI::product(X, PCP11) ^ MI::product(C22, X);
			PCP22 = C22 ^ MI::product(X, C12);
		}
		if (MI::isAX(PBP21, PAP21) && MI::isXA(PCP21, PAP21)) {
			MatrixInt U, V;
			std::uint8_t rank = MI::rank(PAP21, q, p, U, V);
			if (MI::isAX(MI::product(U, PCP22), MI::product(U, PCP21), rank) &&
				MI::isXA(MI::product(PBP11, V), MI::product(PBP21, V), rank))
			{
				MatrixInt Uinv = MI::inv(U, q);
				MatrixInt Vinv = MI::inv(V, p);

				// We maintain blocks of Q, R = Q ^ 2 and S = Q ^ 3.
				MatrixInt Q11 = MI::product(MI::product(Vinv, PAP11), V) ^ MI::create(p);
				MatrixInt Q12 = MI::product(MI::product(Vinv, PAP12), Uinv) ^
					MI::product(MI::product(U, PBP21 ^ MI::product(PAP21, PAP11)), V);
				MatrixInt Q21 = MI::create(rank);
				MatrixInt Q22 = MI::product(MI::product(U, PAP22), Uinv);

				MatrixInt R11 = MI::product(Q11, Q11) ^ MI::product(Q12, Q21);
				MatrixInt R12 = MI::product(Q11, Q12) ^ MI::product(Q12, Q22);
				MatrixInt R21 = MI::product(Q21, Q11) ^ MI::product(Q22, Q21);
				MatrixInt R22 = MI::product(Q21, Q12) ^ MI::product(Q22, Q22);

				MatrixInt S11 = MI::product(Q11, R11) ^ MI::product(Q12, R21);
				MatrixInt S12 = MI::product(Q11, R12) ^ MI::product(Q12, R22);
				MatrixInt S21 = MI::product(Q21, R11) ^ MI::product(Q22, R21);
				MatrixInt S22 = MI::product(Q21, R12) ^ MI::product(Q22, R22);

				// CounterY determines which entry of Q is about to change.
				MatrixInt CounterY = 0;

				bool counterYfound = true;
				while (counterYfound) {
					if (S11 == 0 && S12 == 0 && S21 == 0 && S22 == 0) {
						Matrix Q(n, n, p, p, Q11, Q12, Q21, Q22);
						Matrix W(n, n, p, p, Vinv, 0, 0, U);
						Matrix P(n, n, p, p, MI::create(p), 0, X, MI::create(q));
						idempotents.push_back(A + P * W.inv() * Q * W * P);
						if (find_one) {
							return;
						}
					}
					counterYfound = false;
					for (std::uint8_t i = 0, j0 = rank; i < p && !counterYfound; i++) {
						if (i == rank) {
							j0 = 0;
						}
						for (std::uint8_t j = j0; j < q; j++) {
							MatrixInt Selector = MI::create(i, j);
							if (CounterY & Selector) {
								CounterY ^= Selector;
							}
							else {
								CounterY ^= Selector;
								S11 ^= MI::rightMultiply(Q11, i, j, Q21);
								S12 ^= MI::rightMultiply(Q11, i, j, Q22);
								S21 ^= MI::rightMultiply(Q21, i, j, Q21);
								S22 ^= MI::rightMultiply(Q21, i, j, Q22);

								S11 ^= MI::leftMultiply(R21, i, j);
								S12 ^= MI::leftMultiply(R22, i, j);

								R12 ^= MI::rightMultiply(Q11, i, j);
								R22 ^= MI::rightMultiply(Q21, i, j);

								R11 ^= MI::leftMultiply(Q21, i, j);
								R12 ^= MI::leftMultiply(Q22, i, j);

								S12 ^= MI::rightMultiply(R11, i, j);
								S22 ^= MI::rightMultiply(R21, i, j);

								Q12 ^= Selector;

								counterYfound = true;
								break;
							}
						}
					}
				}
			}
		}
		if (rand) {
			attemptsX++;
			counterXfound = attemptsX < rand;
		}
		else {
			counterXfound = false;
			for (std::uint8_t j = 0; j < p && !counterXfound; j++) {
				for (std::uint8_t i = 0; i < heights[j]; i++) {
					MatrixInt Selector = MI::create(i, j);
					if (CounterX & Selector) {
						CounterX ^= Selector;
					}
					else {
						CounterX ^= Selector;
						X ^= Selector;

						PAP21 ^= MI::leftMultiply(PAP11, i, j);
						PAP22 ^= MI::leftMultiply(PAP12, i, j);
						PAP11 ^= MI::rightMultiply(PAP12, i, j);
						PAP21 ^= MI::rightMultiply(PAP22, i, j);

						PBP21 ^= MI::leftMultiply(PBP11, i, j);
						PBP22 ^= MI::leftMultiply(PBP12, i, j);
						PBP11 ^= MI::rightMultiply(PBP12, i, j);
						PBP21 ^= MI::rightMultiply(PBP22, i, j);

						PCP21 ^= MI::leftMultiply(PCP11, i, j);
						PCP22 ^= MI::leftMultiply(PCP12, i, j);
						PCP11 ^= MI::rightMultiply(PCP12, i, j);
						PCP21 ^= MI::rightMultiply(PCP22, i, j);

						counterXfound = true;
						break;
					}
				}
			}
		}
	}
}

void NilClean::findIdempotents(Matrix const &A, std::uint8_t p,
	std::uint8_t heights[MAX_DIM], bool find_one,
	std::vector<Matrix> &idempotents, std::uint64_t rand)
{
	std::uint8_t n = A.getHeight();
	std::uint8_t q = n - p;

	Matrix B = (A + 1) ^ 3;
	Matrix C = A ^ 3;

	// The outer loop will run through all matrices P = [I0;XI] with X having given heights.
	// Maintained only if rand = 0.
	Matrix P = Matrix(n, n) + 1;

	// We are looking for Y such that (PAP + [IY;00]) ^ 3 = 0.

	// CounterX determines which entry of X is about to change. Relevant if rand = 0.
	Matrix CounterX(n - p, p);

	// Blocks of matrices A, B, C.
	Matrix A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22;
	A.decompose(p, p, A11, A12, A21, A22);
	B.decompose(p, p, B11, B12, B21, B22);
	C.decompose(p, p, C11, C12, C21, C22);

	// Blocks of matrices PAP, PBP, PCP. Maintained only if rand = 0.
	Matrix PAP11, PAP12, PAP21, PAP22, PBP11, PBP12, PBP21, PBP22, PCP11, PCP12, PCP21, PCP22;
	A.decompose(p, p, PAP11, PAP12, PAP21, PAP22);
	B.decompose(p, p, PBP11, PBP12, PBP21, PBP22);
	C.decompose(p, p, PCP11, PCP12, PCP21, PCP22);

	// Allowed nonzero entries of P. Only in use if rand > 0.
	Matrix SelectorP(n, n);
	for (std::uint8_t j = 0; j < p; j++) {
		for (std::uint8_t i = 0; i < heights[j]; i++) {
			SelectorP.setEntry(i + p, j, true);
		}
	}

	// Random seed, generator and distribution.
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<std::uint32_t> distribution(0, 0xFFFFFFFF);

	// Counter of random attempts completed, irrelevant if rand = 0.
	std::uint64_t attemptsX = 0;

	bool counterXfound = true;
	while (counterXfound) {
		if (rand) {
			for (std::uint8_t i = p; i < n; i++) {
				Matrix dP(n, n);
				dP.setRow(i, SelectorP.getRow(i) & distribution(generator));
				P += dP;
			}

			// Compute blocks of PAP, PBP, PCP.
			(P * A * P).decompose(p, p, PAP11, PAP12, PAP21, PAP22);
			(P * B * P).decompose(p, p, PBP11, PBP12, PBP21, PBP22);
			(P * C * P).decompose(p, p, PCP11, PCP12, PCP21, PCP22);
		}
		if (PBP21.isAX(PAP21) && PCP21.isXA(PAP21)) {
			Matrix U, V;
			std::uint8_t rank = PAP21.rank(U, V);
			if ((U * PCP22).isAX(U * PCP21, rank) && (PBP11 * V).isXA(PBP21 * V, rank)) {
				Matrix Q(n, n);
				Q.paste(V.inv() * PAP11 * V + 1, 0, 0);
				Q.paste((V.inv() * PAP12 * U.inv()).submatrix(0, 0, rank, rank) + (U * (PBP21 + PAP21 * PAP11) * V).submatrix(0, 0, rank, rank), 0, p);
				Q.paste(Matrix(rank, rank) + 1, p, 0);
				Q.paste(U * PAP22 * U.inv(), p, p);

				// In the inner loop, we maintain Q and R=Q^2.
				Matrix R = Q * Q;

				// CounterY determines which entry of Q is about to change.
				Matrix CounterY(n, n);

				bool counterYfound = true;
				while (counterYfound) {
					if (Q.zeroProduct(R)) {
						Matrix W = V.inv().blockDiagonal(U);
						idempotents.push_back(A + P * W.inv() * Q * W * P);
						if (find_one) {
							return;
						}
					}
					counterYfound = false;
					for (std::uint8_t i = 0, j0 = p + rank; i < p && !counterYfound; i++) {
						if (i == rank) {
							j0 = p;
						}
						for (std::uint8_t j = j0; j < n; j++) {
							if (CounterY.getEntry(i, j)) {
								CounterY.setEntry(i, j, false);
							}
							else {
								CounterY.setEntry(i, j, true);
								R.addCol(i, Q, j);
								R.addRow(j, Q, i);
								Q.setRow(i, Q.getRow(i) ^
									(static_cast<std::uint32_t>(1) << j));
								counterYfound = true;
								break;
							}
						}
					}
				}
			}
		}
		if (rand) {
			attemptsX++;
			counterXfound = attemptsX < rand;
		}
		else {
			counterXfound = false;
			for (std::uint8_t j = 0; j < p && !counterXfound; j++) {
				for (std::uint8_t i = 0; i < heights[j]; i++) {
					if (CounterX.getEntry(i, j)) {
						CounterX.setEntry(i, j, false);
					}
					else {
						CounterX.setEntry(i, j, true);
						P.setRow(i + p, P.getRow(i + p) ^
							(static_cast<std::uint32_t>(1) << j));

						PAP21.addRow(j, PAP11, i);
						PAP22.addRow(j, PAP12, i);
						PAP11.addCol(i, PAP12, j);
						PAP21.addCol(i, PAP22, j);

						PBP21.addRow(j, PBP11, i);
						PBP22.addRow(j, PBP12, i);
						PBP11.addCol(i, PBP12, j);
						PBP21.addCol(i, PBP22, j);

						PCP21.addRow(j, PCP11, i);
						PCP22.addRow(j, PCP12, i);
						PCP11.addCol(i, PCP12, j);
						PCP21.addCol(i, PCP22, j);

						counterXfound = true;
						break;
					}
				}
			}
		}
	}
}

void NilClean::findIdempotents(Matrix const &A, IntegerSet const &M,
	bool find_one, std::vector<Matrix> &idempotents, std::uint64_t rand)
{
	std::uint8_t n = M.getN();
	std::uint8_t p = M.getSize();
	if (((p & 1) == 1) ^ A.trace()) {
		return;
	}
	Matrix P = M.permutationMatrix();
	std::uint8_t heights[MAX_DIM];
	for (std::uint8_t j = 0; j < p; j++) {
		heights[j] = M.getElement(j) - j;
	}
	std::vector<Matrix> newIdempotents;
	Matrix PTAP = P.transpose() * A * P;
	if (p <= 8 && n - p <= 8) {
//		findIdempotents8(PTAP, p, heights, find_one, newIdempotents, rand);
		findIdempotentsInt(PTAP, p, heights, find_one, newIdempotents, rand);
	}
	else {
		findIdempotents(PTAP, p, heights, find_one, newIdempotents, rand);
	}
	for (Matrix E : newIdempotents) {
		idempotents.push_back(P * E * P.transpose());
	}
}

void NilClean::findIdempotents(Matrix const &A, bool find_one,
	IntegerSet const &M0, std::vector<Matrix> &idempotents,
	std::uint64_t rand, std::uint64_t rand_sets, bool show_steps)
{
	if (rand_sets) {

		// Random seed, generator and distribution.
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::uniform_int_distribution<std::uint32_t>
			distribution(0, (static_cast<std::uint32_t>(1) << A.getHeight()) - 1);

		for (std::uint64_t i = 0; i < rand_sets; i++) {
			std::uint32_t setBinary = distribution(generator);
			std::string setString = "";
			for (std::uint8_t k = 0; k < A.getHeight(); k++) {
				if (setBinary & (static_cast<std::uint32_t>(1) << k)) {
					if (setString != "") {
						setString += ",";
					}
					setString += std::to_string(k);
				}
			}
			IntegerSet M(A.getHeight(), setString);
			if (show_steps) {
				std::cout << "Set = {" << M.toString() << "}\n";
			}
			findIdempotents(A, M, find_one, idempotents, rand);
			if (find_one && !idempotents.empty()) {
				return;
			}
		}
		return;
	}
	IntegerSet M = M0;
	do {
		if (show_steps) {
			std::cout << "Set = {" << M.toString() << "}\n";
		}
		findIdempotents(A, M, find_one, idempotents, rand);
		if (find_one && !idempotents.empty()) {
			return;
		}
	}
	while (M.inc());
}

// Function launched by a thread.
// Given n x n matrix A, find idempotents E such that (A - E) ^ 3 = 0.
// If show_steps is set, write steps in the procedure. Start with the subset with given index.
// Finish if index is the size of the vector of subsets, or if find_one is set and an idempotent
// was found. Append idempotents to the vector idempotents.
// Use mutex to prevent simultaneous reading/writing to index and idempotents,
// and simultaneous use of std::cout.
void findIdempotentsThreadFunction(Matrix A, bool find_one,
	std::vector<IntegerSet> subsets, std::size_t &index,
	std::vector<Matrix> &idempotents, bool show_steps, std::mutex &mutex)
{
	while (true) {
		IntegerSet M(0);
		{
			std::lock_guard<std::mutex> guard(mutex);
			if (index == subsets.size()) {
				return;
			}
			M = subsets[index++];
			if (show_steps) {
				std::cout << "Set = {" << M.toString() << "}; Thread ID: ";
				std::cout << std::this_thread::get_id() << "\n";
			}
		}
		std::vector<Matrix> idempotentsAdded;
		NilClean::findIdempotents(A, M, find_one, idempotentsAdded, 0);
		{
			std::lock_guard<std::mutex> guard(mutex);
			if (!find_one || idempotents.empty()) {
				idempotents.insert(idempotents.end(), idempotentsAdded.begin(), idempotentsAdded.end());
			}
			if (find_one && !idempotents.empty()) {
				index = subsets.size();
			}
		}
	}
}

void NilClean::findIdempotentsParallel(Matrix const &A, bool find_one,
	std::vector<Matrix> &idempotents, bool show_steps)
{
	std::cout << "Number of cores: " << std::thread::hardware_concurrency() << ".\n";
	std::vector<IntegerSet> subsets;
	IntegerSet M(A.getHeight(),"");
	do {
		subsets.push_back(M);
	}
	while (M.inc());
	std::sort(subsets.begin(), subsets.end(),
		[](IntegerSet const &set1, IntegerSet const &set2)
		{ return set1.complexity() > set2.complexity(); });
	std::size_t index = 0;
	std::vector<std::thread> threads;
	std::mutex mutex;
	for (std::uint32_t i = 0; i < std::thread::hardware_concurrency(); i++) {
		threads.push_back(std::thread(findIdempotentsThreadFunction,
			A, find_one, subsets, std::ref(index), std::ref(idempotents),
			show_steps, std::ref(mutex)));
	}
	for (std::size_t i = 0; i < threads.size(); i++) {
		if (threads[i].joinable()) {
			threads[i].join();
		}
	}
}