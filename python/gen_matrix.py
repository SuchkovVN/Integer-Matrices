import sympy as sp
import random

def generate_similar_integer_matrix(A_input, max_attempts=100):
    """
    Generates an integer matrix B similar to A, such that B != A.
    Uses unimodular transformation matrices to ensure integer entries.
    """
    A = sp.Matrix(A_input)
    
    if A.rows != A.cols:
        raise ValueError("Matrix A must be square.")
        
    n = A.rows
    
    # Check if A is a scalar matrix (c * I).
    # A scalar matrix commutes with any matrix, so P^-1 * A * P is always A.
    is_scalar = True
    for i in range(n):
        for j in range(n):
            if i != j and A[i, j] != 0:
                is_scalar = False
            if i == j and A[i, j] != A[0, 0]:
                is_scalar = False
                
    if is_scalar:
        print("Warning: Matrix A is a scalar matrix. No different similar matrix exists.")
        return None, None

    for _ in range(max_attempts):
        # 1. Generate a random unimodular matrix P
        P = sp.eye(n)
        
        # Apply random elementary row operations (determinant remains 1)
        num_operations = random.randint(n, n * 3)
        for _ in range(num_operations):
            i = random.randint(0, n - 1)
            j = random.randint(0, n - 1)
            if i != j:
                multiplier = random.randint(-3, 3)
                if multiplier != 0:
                    # Create elementary matrix
                    E = sp.eye(n)
                    E[i, j] = multiplier
                    # Multiply P by E
                    P = E * P
                    
        # 2. Compute the inverse (guaranteed to be integer since det(P) = 1)
        P_inv = P.inv()
        
        # 3. Calculate B
        B = P_inv * A * P
        
        # 4. Check if B is different from A
        if B != A:
            return B, P
            
    print("Failed to find a different similar matrix within the attempt limit.")
    return None, None

# === EXAMPLE USAGE ===
if __name__ == "__main__":
    # Define an initial integer matrix A
    A_test =[[1,2],
             [0, 1]]
              
    print("Initial Matrix A:")
    sp.pprint(sp.Matrix(A_test))
    
    char_poly_A = sp.Matrix(A_test).charpoly()
    print(f"\nCharacteristic Polynomial of A: {char_poly_A.as_expr()}")
    
    # Generate B
    B_res, P_res = generate_similar_integer_matrix(A_test)
    
    if B_res:
        print("\nGenerated Unimodular Matrix P:")
        sp.pprint(P_res)
        
        print("\nNew Matrix B (B = P^-1 * A * P):")
        sp.pprint(B_res)
        
        char_poly_B = B_res.charpoly()
        print(f"\nCharacteristic Polynomial of B: {char_poly_B.as_expr()}")
        
        if char_poly_A == char_poly_B:
            print("\nVerification Passed: Characteristic polynomials match.")
        else:
            print("\nVerification Failed: Polynomials do not match.")