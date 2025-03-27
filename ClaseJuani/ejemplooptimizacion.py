import numpy as np

from numba import jit

@jit(nopython=True, parallel=True)
def multiply_matrix_by_gravitational_constant(N):
    """
    Creates an NxN matrix filled with ones and multiplies it by the gravitational constant.
    
    Args:
        N (int): Dimension of the square matrix
        
    Returns:
        np.ndarray: NxN matrix multiplied by gravitational constant
    """
    # Gravitational constant (m³/kg·s²)
    G = 6.67430e-11
    
    # Create NxN matrix filled with ones
    matrix = np.ones((N, N))
    
    # Multiply matrix by gravitational constant
    result = matrix * G
    
    return result

# Example usage
def main():
    # Set matrix dimension
    N = 1000
    
    # Create and multiply matrix
    result_matrix = multiply_matrix_by_gravitational_constant(N)
    
    print(f"Matrix shape: {result_matrix.shape}")
    print(f"All elements equal to G: {np.allclose(result_matrix, 6.67430e-11)}")
    print(f"First element value: {result_matrix[0,0]}")

if __name__ == "__main__":
    main()