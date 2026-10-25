import sympy as sp

def build_jordan_matrix(blocks):
    """
    Builds a Jordan Normal Form matrix from a list of blocks.
    Each block is defined by a tuple: (eigenvalue, size).
    """
    jordan_blocks = []
    
    for size, eigenvalue in blocks:
        # Initialize a square zero matrix for the block
        block = sp.zeros(size, size)
        
        # Fill the diagonal and superdiagonal
        for i in range(size):
            block[i, i] =  eigenvalue
            if i < size - 1:
                block[i, i + 1] = 1
                
        jordan_blocks.append(block)
        
    # Combine all blocks into a block diagonal matrix
    return sp.diag(*jordan_blocks)

def get_centralizer_latex(blocks):
    """
    Builds J from block descriptions, calculates the general form 
    of matrix X that commutes with J (XJ = JX), and returns LaTeX.
    """
    # 1. Build Jordan matrix J
    J = build_jordan_matrix(blocks)
    n = J.rows
    
    # 2. Initialize the symbolic matrix X
    X = sp.Matrix(n, n, lambda i, j: sp.Symbol(f'x_{{{i+1},{j+1}}}'))
    
    # 3. Formulate the equation XJ - JX = 0
    commutator = X * J - J * X
    equations = list(commutator)
    
    # 4. Solve the system for all variables in X
    variables = list(X)
    solution = sp.solve(equations, variables)
    
    # 5. Substitute the solution back
    X_general = X.subs(solution)
    
    # 6. Rename free variables to c_1, c_2, ... 
    free_symbols = list(X_general.free_symbols)
    free_symbols.sort(key=lambda sym: str(sym))
    
    c_vars = sp.symbols(f'c_1:{len(free_symbols) + 1}')
    substitution_dict = {old: new for old, new in zip(free_symbols, c_vars)}
    
    X_clean = X_general.subs(substitution_dict)
    
    # 7. Generate LaTeX string
    latex_output = sp.latex(X_clean)
    
    return latex_output, J, X_clean

if __name__ == "__main__":
    # Define Jordan blocks as a list of (eigenvalue, size)
    # Example: one block of size 2 with eigenvalue 2,
    # and one block of size 1 with eigenvalue 3.
    lam = 1
    jordan_structure = [(3, lam), (2, lam), (2, lam), (2, lam), (1, lam), (1, lam)]
    
    print("Jordan Structure Description:")
    for i, (size, eig) in enumerate(jordan_structure):
        print(f"Block {i+1}: Eigenvalue = {eig}, Size = {size}x{size}")
    
    latex_str, J_mat, X_mat = get_centralizer_latex(jordan_structure)
    
    print("\nConstructed Jordan matrix J:")
    sp.pprint(J_mat)
    
    print("\nGeneral form of X (where XJ = JX):")
    sp.pprint(X_mat)
    
    print("\nLaTeX representation:")
    print(latex_str)