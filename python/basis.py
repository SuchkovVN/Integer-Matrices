import sympy as sp

def kronecker_product(A_input, B_input):
    """
    Computes the Kronecker product of two matrices A and B.
    """
    A = sp.Matrix(A_input)
    B = sp.Matrix(B_input)
    
    return sp.Matrix(
        A.rows * B.rows, 
        A.cols * B.cols, 
        lambda i, j: A[i // B.rows, j // B.cols] * B[i % B.rows, j % B.cols]
    )

def vec_to_matrix(vec, m, n):
    """
    Converts a column vector back to an m x n matrix (undoes vectorization by columns).
    """
    mat = sp.zeros(m, n)
    for col in range(n):
        for row in range(m):
            mat[row, col] = vec[col * m + row]
    return mat

def smith_normal_form_uv(A_input):
    """
    Reduces matrix A to diagonal form S and finds unimodular matrices U and V,
    such that S = U * A * V.
    """
    A = sp.Matrix(A_input)
    m, n = A.shape
    S = A.copy()
    U = sp.eye(m)
    V = sp.eye(n)
    
    for i in range(min(m, n)):
        pivot_row, pivot_col = -1, -1
        for r in range(i, m):
            for c in range(i, n):
                if S[r, c] != 0:
                    pivot_row, pivot_col = r, c
                    break
            if pivot_row != -1: break
                
        if pivot_row == -1:
            break 
            
        if pivot_row != i:
            S.row_swap(i, pivot_row)
            U.row_swap(i, pivot_row)
        if pivot_col != i:
            S.col_swap(i, pivot_col)
            V.col_swap(i, pivot_col)
            
        changed = True
        while changed:
            changed = False
            
            for j in range(i + 1, m):
                if S[j, i] != 0:
                    changed = True
                    q = S[j, i] // S[i, i]
                    for c in range(n): S[j, c] -= q * S[i, c]
                    for c in range(m): U[j, c] -= q * U[i, c]
                    
                    if S[j, i] != 0:
                        S.row_swap(i, j)
                        U.row_swap(i, j)

            for j in range(i + 1, n):
                if S[i, j] != 0:
                    changed = True
                    q = S[i, j] // S[i, i]
                    for r in range(m): S[r, j] -= q * S[r, i]
                    for r in range(n): V[r, j] -= q * V[r, i]
                    
                    if S[i, j] != 0:
                        S.col_swap(i, j)
                        V.col_swap(i, j)
                        
    for i in range(min(m, n)):
        if S[i, i] < 0:
            for c in range(n): S[i, c] = -S[i, c]
            for c in range(m): U[i, c] = -U[i, c]
            
    return S, U, V

def solve_linear_diophantine(A_mat, b_vec):
    """
    Finds integer solutions for the system Ax = b using Smith Normal Form.
    """
    A = sp.Matrix(A_mat)
    b = sp.Matrix(b_vec)
    m, n = A.shape
    
    S, U, V = smith_normal_form_uv(A)
    b_prime = U * b
    
    y_solution = sp.zeros(n, 1)
    free_vars = []
    k_syms = sp.symbols(f'x0:{n}')
    k_idx = 0
    
    for i in range(n):
        if i < m and S[i, i] != 0:
            if b_prime[i] % S[i, i] != 0:
                print(f"No integer solutions: equation {i} requires dividing {b_prime[i]} by {S[i, i]}.")
                return None, None
            else:
                y_solution[i, 0] = b_prime[i] // S[i, i]
        else:
            if i < m and b_prime[i] != 0:
                print(f"No integer solutions: contradiction in row {i} (0 != {b_prime[i]}).")
                return None, None
            
            y_solution[i, 0] = k_syms[k_idx]
            free_vars.append(k_syms[k_idx])
            k_idx += 1
            
    for i in range(n, m):
        if b_prime[i] != 0:
            print(f"No integer solutions: inconsistency in redundant equation {i}.")
            return None, None
            
    x_solution = V * y_solution
    zero_subs = {k: 0 for k in free_vars}
    x_particular = x_solution.subs(zero_subs)
    basis_vectors = {k: x_solution.diff(k) for k in free_vars}
    
    return x_particular, basis_vectors

def solve_sylvester_integer(A_mat, B_mat):
    """
    Solves the matrix equation AX = XB in integers.
    Returns the basis matrices for the null space.
    """
    A = sp.Matrix(A_mat)
    B = sp.Matrix(B_mat)
    m = A.rows
    n = B.rows
    
    I_n = sp.eye(n)
    I_m = sp.eye(m)
    
    term1 = kronecker_product(I_n, A)
    term2 = kronecker_product(B.T, I_m)
    M = term1 - term2
    
    zero_vec = sp.zeros(m * n, 1)
    
    _, basis_vecs = solve_linear_diophantine(M, zero_vec)
    
    basis_matrices = []
    for k, vec in basis_vecs.items():
        mat = vec_to_matrix(vec, m, n)
        basis_matrices.append((k, mat))
        
    return basis_matrices

def analyze_jordan_and_commutator(A_input, B_input):
    """
    Computes Jordan forms for A and B, verifies they share the same structure,
    outputs the transformation matrices, and finds the integer basis for JX = XJ.
    """
    A = sp.Matrix(A_input)
    B = sp.Matrix(B_input)
    
    print("--- Computing Jordan Normal Forms & Transformation Matrices ---")
    P_A, J_A = A.jordan_form()
    P_B, J_B = B.jordan_form()
    
    print("Transformation Matrix P_A (where A = P_A * J_A * P_A^-1):")
    sp.pprint(P_A)
    
    print("\nJordan Form of A (J_A):")
    sp.pprint(J_A)
    
    print("\n" + "-"*30)
    
    print("\nTransformation Matrix P_B (where B = P_B * J_B * P_B^-1):")
    sp.pprint(P_B)
    
    print("\nJordan Form of B (J_B):")
    sp.pprint(J_B)
    
    if J_A != J_B:
        print("\nWarning: J_A and J_B differ (possibly due to Jordan block permutations).")
        print("Using J_A as the target Jordan form J.")
    else:
        print("\nMatrices A and B share the exact same Jordan form J.")
        
    J = J_A
    print("\n--- Solving JX = XJ over integers ---")
    
    basis = solve_sylvester_integer(J, J)
    
    # Returning J, P_A, P_B, and the basis
    return J, P_A, P_B, basis

# === EXAMPLE USAGE ===
if __name__ == "__main__":
    A_test = [[0, 92, -91, -13], 
              [-4, 56, -52,  0 ],
              [-4, 49, -45, 1],
              [0, 25, -25, -3]] 
              
    B_test = [[-4, 1, -2, 10],
              [-26, 6, -8, 42],
              [-25, 4, -6, 41],
              [-6, 1, -2, 12]] 
    
    print("Matrix A:")
    sp.pprint(sp.Matrix(A_test))
    
    print("\nMatrix B:")
    sp.pprint(sp.Matrix(B_test))
    print("\n" + "="*40 + "\n")
    
    J, P_A, P_B, basis = analyze_jordan_and_commutator(A_test, B_test)
    
    if basis:
        print("\nBasis matrices for the Z-module of solutions (X):")
        res_mat = sp.zeros(basis[0][1].shape[0])
        for sym, mat in basis:
            print(f"\nBasis element for {sym}:")
            res_mat += sym * mat
            sp.pprint(mat)

        print("\n Resulting X matrix")
        sp.pprint(res_mat)

        print("\n Determinant")
        print(sp.det(res_mat))
    else:
        print("\nOnly the trivial solution exists (X = 0).")
