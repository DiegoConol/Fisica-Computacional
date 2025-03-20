#include <stdio.h>
#include <stdlib.h>

// Función para leer el tamaño de la matriz del archivo
int leerDimension(FILE *archivo) {
    int dimension;
    if (fscanf(archivo, "%d", &dimension) != 1) {
        fprintf(stderr, "Error leyendo la dimensión de la matriz\n");
        return -1;
    }
    return dimension;
}

// Función para leer la matriz del archivo
int **leerMatriz(FILE *archivo, int n) {
    int **matriz = (int **)malloc(n * sizeof(int *));
    if (!matriz) {
        fprintf(stderr, "Error reservando memoria\n");
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        matriz[i] = (int *)malloc(n * sizeof(int));
        if (!matriz[i]) {
            for (int j = 0; j < i; j++) free(matriz[j]);
            free(matriz);
            fprintf(stderr, "Error reservando memoria\n");
            return NULL;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(archivo, "%d", &matriz[i][j]) != 1) {
                fprintf(stderr, "Error leyendo elemento [%d,%d]\n", i, j);
                for (int k = 0; k <= i; k++) free(matriz[k]);
                free(matriz);
                return NULL;
            }
        }
    }

    return matriz;
}

// Función para imprimir la matriz
void imprimirMatriz(int **matriz, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matriz[i][j]);
        }
        printf("\n");
    }
}

// Función para liberar la memoria de la matriz
void liberarMatriz(int **matriz, int n) {
    for (int i = 0; i < n; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

int main() {
    FILE *archivo;
    int **matriz;
    int n;

    // Abrir el archivo
    archivo = fopen("fichero.txt", "r"); // Cambiado de matriz.txt a fichero.txt
    if (!archivo) {
        fprintf(stderr, "Error abriendo el archivo fichero.txt\n");
        return 1;
    }

    // Leer la dimensión
    n = leerDimension(archivo);
    if (n == -1) {
        fclose(archivo);
        return 1;
    }

    // Leer la matriz
    matriz = leerMatriz(archivo, n);
    if (!matriz) {
        fclose(archivo);
        return 1;
    }

    // Imprimir la matriz
    printf("Matriz leída:\n");
    imprimirMatriz(matriz, n);

    // Liberar recursos
    liberarMatriz(matriz, n);
    fclose(archivo);

    return 0;
}