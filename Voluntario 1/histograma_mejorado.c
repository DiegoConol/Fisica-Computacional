/*
Este programa usa los datos de velocidades.txt para hacer un histograma de velocidades
y compararlo con la distribución de Maxwell-Boltzmann. Se analizan los datos
entre los tiempos t=20 y t=50.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

// Constantes del sistema
#define Epsilon 1.0         // Constante de Unidades de potencial
#define Sigma 1.0           // Constante de distancia
#define KB 1.0              // Constante de Boltzmann
#define N 50                // Número de partículas
#define L 10.0              // Longitud de la caja LXL
#define M 1.0               // Masa de las partículas
#define h 0.002             // Paso temporal
#define PI 3.14159265       // Pi
#define T_TOTAL 60          // Tiempo total de simulación

// Tiempos para el histograma
#define T_MIN 20            // Tiempo mínimo para histograma
#define T_MAX 50            // Tiempo máximo para histograma

// Parámetros del histograma
#define NUM_BINS 30         // Número de intervalos para el histograma
#define V_MAX 10.0          // Velocidad máxima esperada para el histograma

int numpasos = (int)(T_TOTAL/h);    // Número de pasos temporales
int t_min_idx = (int)(T_MIN/h);     // Índice correspondiente a T_MIN
int t_max_idx = (int)(T_MAX/h);     // Índice correspondiente a T_MAX

// Función para crear un arreglo dinámico tridimensional
double ***crear_arreglo_dinamico(int numparticulas, int numpasos) 
{
    double ***arreglo = (double ***)malloc(numparticulas * sizeof(double **));
    if (arreglo == NULL) {
        printf("Error de memoria al crear arreglo.\n");
        exit(1);
    }

    for (int i = 0; i < numparticulas; i++) {
        arreglo[i] = (double **)malloc(2 * sizeof(double *));
        if (arreglo[i] == NULL) {
            printf("Error de memoria al crear arreglo[%d].\n", i);
            exit(1);
        }
        for (int j = 0; j < 2; j++) {
            arreglo[i][j] = (double *)malloc(numpasos * sizeof(double));
            if (arreglo[i][j] == NULL) {
                printf("Error de memoria al crear arreglo[%d][%d].\n", i, j);
                exit(1);
            }
        }
    }
    return arreglo;
}

// Función para liberar la memoria de un arreglo dinámico
void liberar_arreglo_dinamico(double ***arreglo, int numparticulas) 
{
    for (int i = 0; i < numparticulas; i++) {
        for (int j = 0; j < 2; j++) {
            free(arreglo[i][j]);
        }
        free(arreglo[i]);
    }
    free(arreglo);
}

// Función para calcular la distribución de Maxwell-Boltzmann en 2D
double maxwell_2D(double v, double T) 
{
    // f(v) = (m*v/kT) * exp(-m*v²/2kT)
    return (M*v/KB/T) * exp(-M*v*v/(2*KB*T));
}

int main(void)
{
    // Abrir archivos
    FILE *veltxt = fopen("velocidades.txt", "r");       // Fichero de velocidades
    FILE *cinetica = fopen("cinetica.txt", "r");        // Fichero de energía cinética
    FILE *histograma = fopen("histograma.txt", "w");    // Fichero del histograma
    FILE *maxwelltxt = fopen("maxwell.txt", "w");       // Fichero de distribución teórica

    if (cinetica == NULL || veltxt == NULL || histograma == NULL || maxwelltxt == NULL) {
        printf("Error al abrir uno o más archivos.\n");
        return 1;
    }

    // Crear arreglos dinámicos
    double ***vel = crear_arreglo_dinamico(N, numpasos);
    double K[numpasos];  // Energía cinética en cada paso

    // Leer energía cinética desde cinetica.txt
    for (int t = 0; t < numpasos; t++) {
        if (fscanf(cinetica, "%lf", &K[t]) != 1) {
            printf("Error al leer el archivo cinetica.txt en el paso %d.\n", t);
            fclose(cinetica);
            fclose(veltxt);
            return 1;
        }
    }

    // Leer velocidades desde velocidades.txt
    for (int t = 0; t < numpasos; t++) {
        for (int i = 0; i < N; i++) {
            if (fscanf(veltxt, "%lf, %lf", &vel[i][0][t], &vel[i][1][t]) != 2) {
                printf("Error al leer el archivo velocidades.txt en el paso %d, partícula %d.\n", t, i);
                fclose(cinetica);
                fclose(veltxt);
                return 1;
            }
        }
        // Saltar la línea en blanco entre pasos temporales
        fgetc(veltxt);
    }

    fclose(cinetica);
    fclose(veltxt);

    // Calcular la temperatura media en el intervalo [T_MIN, T_MAX]
    double suma_energia_cinetica = 0.0;
    int num_muestras = (t_max_idx - t_min_idx) * N;
    
    for (int t = t_min_idx; t < t_max_idx; t++) {
        suma_energia_cinetica += K[t];
    }
    
    double K_media = suma_energia_cinetica / (t_max_idx - t_min_idx);
    double T_media = K_media / (N * KB);
    
    printf("Temperatura media en [%d, %d]: %lf\n", T_MIN, T_MAX, T_media);

    // Crear el histograma
    double bin_width = V_MAX / NUM_BINS;
    int histograma_datos[NUM_BINS];
    double velocidades_totales[N * (t_max_idx - t_min_idx)];
    int contador_velocidades = 0;
    
    // Inicializar el histograma
    for (int i = 0; i < NUM_BINS; i++) {
        histograma_datos[i] = 0;
    }
    
    // Recopilar todas las velocidades en el intervalo y calcular el histograma
    for (int t = t_min_idx; t < t_max_idx; t++) {
        for (int i = 0; i < N; i++) {
            double v_modulo = sqrt(vel[i][0][t] * vel[i][0][t] + vel[i][1][t] * vel[i][1][t]);
            velocidades_totales[contador_velocidades++] = v_modulo;
            
            // Asignar a un bin del histograma
            int bin = (int)(v_modulo / bin_width);
            if (bin >= 0 && bin < NUM_BINS) {
                histograma_datos[bin]++;
            }
        }
    }
    
    // Normalizar el histograma para compararlo con la distribución teórica
    double factor_normalizacion = num_muestras * bin_width;
    
    // Escribir datos del histograma y la curva teórica de Maxwell
    fprintf(histograma, "# Velocidad\tHistograma\tMaxwell-Boltzmann\n");
    for (int i = 0; i < NUM_BINS; i++) {
        double v_central = (i + 0.5) * bin_width;  // Valor central del bin
        double histograma_normalizado = (double)histograma_datos[i] / factor_normalizacion;
        double maxwell_valor = maxwell_2D(v_central, T_media);
        
        fprintf(histograma, "%lf\t%lf\t%lf\n", v_central, histograma_normalizado, maxwell_valor);
        fprintf(maxwelltxt, "%lf\t%lf\n", v_central, maxwell_valor);
    }
    
    // Escribir un archivo separado con todas las velocidades para análisis adicionales
    FILE *todas_velocidades = fopen("todas_velocidades.txt", "w");
    if (todas_velocidades != NULL) {
        for (int i = 0; i < contador_velocidades; i++) {
            fprintf(todas_velocidades, "%lf\n", velocidades_totales[i]);
        }
        fclose(todas_velocidades);
    }
    
    
    
    fclose(histograma);
    fclose(maxwelltxt);
    liberar_arreglo_dinamico(vel, N);
    
    printf("Histograma creado exitosamente entre t=%d y t=%d.\n", T_MIN, T_MAX);
    printf("La temperatura media calculada es: %lf\n", T_media);
    
    return 0;
}