//Primero: un programa que lea matrices

//Luego las matrices de reescala

//Probamos optimización. O bien ompmD o O1-O2-O3.

//Si lo hacemos en python lo hacemos con numba.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

// Structure to represent a matrix
typedef struct {
    int **data;
    size_t rows;
    size_t cols;
} Matrix;

// Function to read matrix from a file
Matrix* read_matrix(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file: %s\n", filename);
        return NULL;
    }

    // Read dimensions
    size_t rows, cols;
    if (fscanf(file, "%zu %zu", &rows, &cols) != 2) {
        fclose(file);
        fprintf(stderr, "Invalid matrix dimensions in file: %s\n", filename);
        return NULL;
    }

    // Allocate memory for matrix
    int **matrix_data = (int **)calloc(rows, sizeof(int *));
    if (!matrix_data) {
        fclose(file);
        fprintf(stderr, "Memory allocation failed for matrix\n");
        return NULL;
    }

    // Allocate rows
    for (size_t i = 0; i < rows; i++) {
        matrix_data[i] = (int *)calloc(cols, sizeof(int));
        if (!matrix_data[i]) {
            // Free previously allocated rows
            for (size_t j = 0; j < i; j++) {
                free(matrix_data[j]);
            }
            free(matrix_data);
            fclose(file);
            fprintf(stderr, "Memory allocation failed for row %zu\n", i + 1);
            return NULL;
        }
    }

    // Read matrix elements
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            if (fscanf(file, "%d", &matrix_data[i][j]) != 1) {
                // Free all allocated memory
                for (size_t k = 0; k <= i; k++) {
                    free(matrix_data[k]);
                }
                free(matrix_data);
                fclose(file);
                fprintf(stderr, "Invalid matrix element at row %zu, col %zu\n", i + 1, j + 1);
                return NULL;
            }
        }
    }

    fclose(file);

    // Create and populate Matrix structure
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));
    result->data = matrix_data;
    result->rows = rows;
    result->cols = cols;

    return result;
}

// Function to free matrix memory
void free_matrix(Matrix *matrix) {
    if (matrix && matrix->data) {
        for (size_t i = 0; i < matrix->rows; i++) {
            free(matrix->data[i]);
        }
        free(matrix->data);
    }
    free(matrix);
}

// Function to print matrix
void print_matrix(const Matrix *matrix) {
    if (!matrix || !matrix->data) return;

    for (size_t i = 0; i < matrix->rows; i++) {
        for (size_t j = 0; j < matrix->cols; j++) {
            printf("%d ", matrix->data[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    // Check if directory path was provided
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <directory_path>\n", argv[0]);
        return 1;
    }

    DIR *dir;
    struct dirent *ent;
    char filepath[1024];

    // Open directory
    dir = opendir(argv[1]);
    if (!dir) {
        fprintf(stderr, "Could not open directory: %s\n", argv[1]);
        return 1;
    }

    printf("Reading matrices from directory: %s\n", argv[1]);

    while ((ent = readdir(dir)) != NULL) {
        // Skip current directory and parent directory entries
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
            continue;
        }

        // Construct full path to file
        snprintf(filepath, sizeof(filepath), "%s/%s", argv[1], ent->d_name);

        printf("\nReading file: %s\n", ent->d_name);

        Matrix *matrix = read_matrix(filepath);
        if (matrix) {
            printf("Matrix dimensions: %zx x %zx\n", matrix->rows, matrix->cols);
            printf("Matrix contents:\n");
            print_matrix(matrix);
            
            free_matrix(matrix);
        }
    }

    closedir(dir);
    return 0;
}