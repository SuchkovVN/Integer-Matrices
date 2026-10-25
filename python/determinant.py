import sympy as sp
from gen_matrix_jordan import build_similar_matrices

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
    Reduces matrix A to diagonal form S and finds unimodular matrices U and V.
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
                return None, None
            else:
                y_solution[i, 0] = b_prime[i] // S[i, i]
        else:
            if i < m and b_prime[i] != 0:
                return None, None
            
            y_solution[i, 0] = k_syms[k_idx]
            free_vars.append(k_syms[k_idx])
            k_idx += 1
            
    for i in range(n, m):
        if b_prime[i] != 0:
            return None, None
            
    x_solution = V * y_solution
    basis_vectors = {k: x_solution.diff(k) for k in free_vars}
    
    return basis_vectors

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
    
    basis_vecs = solve_linear_diophantine(M, zero_vec)
    
    basis_matrices = []
    if basis_vecs:
        for k, vec in basis_vecs.items():
            mat = vec_to_matrix(vec, m, n)
            basis_matrices.append((k, mat))
        
    return basis_matrices

def analyze_solution_and_determinant(A_input, B_input):
    """
    Finds the basis for AX = XB, constructs the general solution matrix X,
    and calculates its symbolic determinant.
    """
    A = sp.Matrix(A_input)
    B = sp.Matrix(B_input)
    
    print("--- Solving AX = XB over integers ---")
    basis = solve_sylvester_integer(A, B)
    
    if not basis:
        print("\nOnly the trivial solution exists (X = 0).")
        return
        
    print("\nBasis matrices for the Z-module of solutions (X):")
    
    # Construct the general matrix X dynamically
    m, n = A.rows, B.rows
    X_general = sp.zeros(m, n)
    
    for sym, mat in basis:
        print(f"\nBasis element for parameter {sym}:")
        sp.pprint(mat)
        X_general = X_general + sym * mat
        
    print("\n--- General Form of Matrix X ---")
    sp.pprint(X_general)
    
    # Calculate the symbolic determinant if X is square
    print("\n--- Determinant Analysis ---")
    if m == n:
        print("Calculating the symbolic determinant of the general matrix X...")
        det_X = X_general.det()
        factored_det = sp.factor(det_X)
        print("\nDeterminant of X in terms of basis variables (factored form):")
        sp.pprint(factored_det)
        print("\nLaTeX:")
        print(sp.latex(factored_det))
    else:
        print(f"Matrix X is of size {m}x{n}. Determinant is only defined for square matrices.")

# === EXAMPLE USAGE ===
if __name__ == "__main__":
    J_A = [[1, 1, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]]
    
    A, B, P_A, P_B = build_similar_matrices(J_A)

    # A = [[-5, 6, 0, -3],
    #      [0, 1, 0, 0],
    #      [8, -8, 1, 4],
    #      [12, -12, 0, 7]]
    # B = [[2, 1, -1, 0],
    #      [-1, 0, 1, 0],
    #      [0, 0, 1, 0],
    #      [1, 1, -1, 1]]


    # J_A = sp.Matrix(A).jordan_form(False)
    
    print("Matrix A:")
    sp.pprint(sp.Matrix(A))
    
    print("\nMatrix B:")
    sp.pprint(sp.Matrix(B))
    print("\n" + "="*50 + "\n")

    print("Matrix A jordan form:")
    sp.pprint(sp.Matrix(J_A))
    
    analyze_solution_and_determinant(A, B)