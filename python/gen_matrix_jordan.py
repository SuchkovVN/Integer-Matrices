import sympy as sp
import random

def generate_unimodular_matrix(n, num_operations=15):
    """
    Generates an n x n unimodular matrix (determinant = 1 or -1)
    by applying random elementary row operations to an identity matrix.
    """
    P = sp.eye(n)
    for _ in range(num_operations):
        i = random.randint(0, n - 1)
        j = random.randint(0, n - 1)
        if i != j:
            multiplier = random.randint(-3, 3)
            if multiplier != 0:
                # Add (multiplier * row j) to row i
                E = sp.eye(n)
                E[i, j] = multiplier
                P = E * P
                
        # Occasionally multiply a row by -1 to flip the determinant sign
        if random.random() < 0.2:
            row_idx = random.randint(0, n - 1)
            E_neg = sp.eye(n)
            E_neg[row_idx, row_idx] = -1
            P = E_neg * P
            
    return P

def is_scalar_matrix(M):
    """
    Checks if matrix M is a scalar matrix (c * I).
    """
    n = M.rows
    for i in range(n):
        for j in range(n):
            if i != j and M[i, j] != 0:
                return False
            if i == j and M[i, j] != M[0, 0]:
                return False
    return True

def build_similar_matrices(J_input, max_attempts=100):
    """
    Given a Jordan form J, builds two arbitrary, distinct matrices A and B
    that are similar to J (i.e., they have J as their Jordan form).
    """
    J = sp.Matrix(J_input)
    
    if J.rows != J.cols:
        raise ValueError("Jordan matrix J must be square.")
        
    n = J.rows
    
    if is_scalar_matrix(J):
        print("Warning: The provided matrix J is a scalar matrix.")
        print("It is impossible to generate distinct similar matrices for a scalar matrix.")
        return None, None, None, None

    A = None
    B = None
    P_A = None
    P_B = None
    
    # 1. Generate matrix A
    for _ in range(max_attempts):
        P_temp = generate_unimodular_matrix(n)
        A_temp = P_temp * J * P_temp.inv()
        
        if A_temp != J:
            A = A_temp
            P_A = P_temp
            break
            
    if A is None:
        print("Failed to generate a distinct matrix A within the attempt limit.")
        return None, None, None, None

    # 2. Generate matrix B (must be distinct from J and A)
    for _ in range(max_attempts):
        P_temp = generate_unimodular_matrix(n)
        B_temp = P_temp * J * P_temp.inv()
        
        if B_temp != J and B_temp != A:
            B = B_temp
            P_B = P_temp
            break
            
    if B is None:
        print("Failed to generate a distinct matrix B within the attempt limit.")
        return A, None, P_A, None

    return A, B, P_A, P_B

# === EXAMPLE USAGE ===
if __name__ == "__main__":
    # Define an integer Jordan Normal Form J
    # Eigenvalue 2 with a Jordan block of size 2, and eigenvalue 3
    J_test = [[1, 1, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]]
              
    print("Provided Jordan Matrix J:")
    sp.pprint(sp.Matrix(J_test))
    
    print("\n" + "="*40 + "\n")
    print("Generating matrices A and B...")
    
    A, B, P_A, P_B = build_similar_matrices(J_test)
    
    if A is not None and B is not None:
        print("\n--- MATRIX A ---")
        print("Transformation Matrix P_A:")
        sp.pprint(P_A)
        print("\nMatrix A (A = P_A * J * P_A^-1):")
        sp.pprint(A)
        
        print("\n--- MATRIX B ---")
        print("Transformation Matrix P_B:")
        sp.pprint(P_B)
        print("\nMatrix B (B = P_B * J * P_B^-1):")
        sp.pprint(B)
        
        # Verify Jordan forms back (sympy's jordan_form() returns (P, J))
        _, J_from_A = A.jordan_form()
        _, J_from_B = B.jordan_form()