#include "pari_wrapper/pari.h"
#include "utils/arg_parse.h"

#include <iostream>

int main(int argc, char* argv[]) {
    StrMatrix str_A, str_B;
    Params params;
    int64_t pari_size;
    int64_t pari_max_prime;

    params.createParam<StrMatrix>("--mat_A", {}, &str_A, "Matrix A");
    params.createParam<StrMatrix>("--mat_B", {}, &str_B, "Matrix B");
    params.createParam<int64_t>("--pari_size", 100'000, &pari_size, "[PARI] Memory size");
    params.createParam<int64_t>("--pari_max_prime", 4, &pari_max_prime, "[PARI] Max prime");

    if (!parseParams(argv, argc, params)) {
        std::cout << "Failed to parse cli args";
        return -1;
    }

    if (str_A.empty()) {
        std::cout << "Invalid matrix A\n";
        return -1;
    }
    if (str_B.empty()) {
        std::cout << "Invalid matrix B\n";
        return -1;
    }

    auto p_init = pari::PariInitializer(pari_size, pari_max_prime);

    auto matA = pari::strzmat_to_zmat(str_A);
    auto matB = pari::strzmat_to_zmat(str_B);

    auto basis = pari::zmat_eq_basis(matA, matB);

    if (basis.empty()) {
        std::cout << "Basis is empty\n";
        return 0;
    }

    for (auto&& x : basis) {
        outmat(x);
    }

    return 0;
}