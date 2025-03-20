#include <stdio.h>
#include <stdlib.h>

/**
 * Updates a 2D array by multiplying every element by 2.
 * 
 * @param data Pointer to the start of the 2D array
 * @param rows Number of rows in the array
 * @param cols Number of columns in the array
 */
void multiply_matrix_by_two(float *data, int rows, int cols) {
    // Calculate total number of elements
    int total_elements = rows * cols;
    
    // Iterate through all elements and multiply by 2
    for (int i = 0; i < total_elements; i++) {
        data[i] *= 2;
    }
}

int main() {
    // Create sample 2D array
    int rows = 3;
    int cols = 4;
    float *matrix = malloc(rows * cols * sizeof(float));
    
    // Initialize matrix with test values
    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = i + 1.0f;
    }
    
    printf("Original matrix:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.1f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    
    multiply_matrix_by_two(matrix, rows, cols);
    
    printf("\nMatrix after multiplication by 2:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.1f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    
    free(matrix);
    return 0;
}