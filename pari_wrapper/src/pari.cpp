#include "pari_wrapper/pari.h"

namespace pari {

GEN kron(GEN x, GEN y) {
    if (typ(x) != t_MAT) {
        pari_err(e_TYPE, "kron", x);
    }
    if (typ(y) != t_MAT) {
        pari_err(e_TYPE, "kron", y);
    }

    auto x_num_rows = nbrows(x);
    auto x_num_cols = lg(x) - 1;
    auto y_num_rows = nbrows(y);
    auto y_num_cols = lg(y) - 1;
    auto res_num_rows = x_num_rows * x_num_rows;
    auto res_num_cols = x_num_cols * y_num_cols;

    auto res = zeromatcopy(res_num_rows, res_num_cols);

    for (int64_t i = 0; i < res_num_rows; ++i) {
        for (int64_t j = 0; j < res_num_cols; ++j) {
            gcoeff(res, i + 1, j + 1) = gmul(gcoeff(x, (i / y_num_rows) + 1, (j / y_num_cols) + 1), gcoeff(y, (i % y_num_rows) + 1, (j % y_num_cols) + 1));
        }
    }

    return res;
}

GEN strzmat_to_zmat(std::vector<std::vector<std::string>> str_mat) {
    int64_t m = str_mat.size();
    int64_t n = str_mat[0].size();

    auto res = zeromatcopy(m, n);

    for (int64_t i = 0; i < m; ++i) {
        for (int64_t j = 0; j < n; ++j) {
            if (str_mat[i][j][0] == '-') {
                auto t = str_mat[i][j].substr(1);
                gcoeff(res, i + 1, j + 1) = gneg(strtoi(t.data()));
            } else {
                gcoeff(res, i + 1, j + 1) = strtoi(str_mat[i][j].data());
            }
        }
    }

    return res;
}

std::vector<GEN> zmat_eq_basis(GEN x, GEN y) {
    if (typ(x) != t_MAT) {
        pari_err(e_TYPE, "zmat_eq_basis", x);
    }
    if (typ(y) != t_MAT) {
        pari_err(e_TYPE, "zmat_eq_basis", y);
    }

    auto x_num_rows = nbrows(x);
    auto x_num_cols = lg(x) - 1;
    auto y_num_rows = nbrows(y);
    auto y_num_cols = lg(y) - 1;

    auto E = matid(x_num_rows);

    auto S = ZM_sub(pari::kron(E, x), pari::kron(zm_transpose(y), E));

    auto m = nbrows(S);
    auto n = lg(S) - 1;

    auto L = zeromatcopy(m, m);
    auto R = zeromatcopy(n, n);

    auto Ssf = ZM_snfall(S, &L, &R);

    int64_t basis_size = 0;
    while (gcoeff(Ssf, basis_size + 1, basis_size + 1) == gen_0) {
        ++basis_size;
    }

    std::vector<GEN> res(basis_size);
    for (int64_t i = 0; i < basis_size; ++i) {
        res[i] = zeromatcopy(x_num_cols, y_num_rows);
        for (int64_t col = 1; col <= y_num_rows; ++col) {
            for (int64_t row = 1; row <= x_num_cols; ++row) {
                gcoeff(res[i], row, col) = gcoeff(R, row + (col - 1) * x_num_cols, i + 1);
            }
        }
    }

    return res;
}

}  // namespace pari
