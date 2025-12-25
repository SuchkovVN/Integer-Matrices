#pragma once

#include <pari/pari.h>

#include <vector>
#include <string>

namespace pari {
class PariInitializer {
public:
    PariInitializer(size_t memory_space, size_t max_prime) {
        pari_init(memory_space, max_prime);
    }
    ~PariInitializer() {
        pari_close();
    }
};

GEN kron(GEN x, GEN y);

GEN strzmat_to_zmat(std::vector<std::vector<std::string>> str_mat);

std::vector<GEN>  zmat_eq_basis(GEN x, GEN y);

}  // namespace pari