#include <pari/pari.h>

#include "pari_wrapper/pari.h"

#include <iostream>
#include <vector>

int main() {
    pari_init(1'000'000, 100);

    std::vector<std::vector<int64_t>> matA{ { 0, 2, 3 }, { 0, 0, 3 }, { 0, 0, 0 } };
    std::vector<std::vector<int64_t>> matB{ { 0, 3, 3 }, { 0, 0, 2 }, { 0, 0, 0 } };
    int64_t mA = matA.size(), nA = matA[0].size();
    int64_t mB = matB.size(), nB = matB[0].size();

    auto A = zeromatcopy(mA, nA);
    auto B = zeromatcopy(mB, nA);

    for (int64_t i = 1; i <= mA; ++i) {
        for (int64_t j = 1; j <= nA; ++j) {
            gcoeff(A, i, j) = stoi(matA[i - 1][j - 1]);
        }
    }
    for (int64_t i = 1; i <= mB; ++i) {
        for (int64_t j = 1; j <= nB; ++j) {
            gcoeff(B, i, j) = stoi(matB[i - 1][j - 1]);
        }
    }

    auto E = matid(mA);

    auto S = ZM_sub(pari::kron(E, A), pari::kron(zm_transpose(B), E));

    auto m = nbrows(S);
    auto n = lg(S) - 1;

    std::cout << "System matrix: \n";
    outmat(S);

    auto L = zeromatcopy(m, m);
    auto R = zeromatcopy(n, n);

    auto Ssf = ZM_snfall(S, &L, &R);

    std::cout << "SNF: \n";
    outmat(Ssf);
    std::cout << "L: \n";
    outmat(L);
    std::cout << "R: \n";
    outmat(R);

    int64_t basis_size = 0;
    while (gcoeff(Ssf, basis_size + 1, basis_size + 1) == gen_0) {
        ++basis_size;
    }

    std::cout << "Basis size: " << basis_size << '\n';

    auto basis = zeromatcopy(n, basis_size);
    for (int64_t i = 1; i <= basis_size; ++i) {
        gel(basis, i) = gel(R, i);
    }

    std::cout << "Basis: \n";
    outmat(basis);

    pari_close();
    return 0;
}