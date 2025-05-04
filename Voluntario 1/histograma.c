/*

Este programa usa los datos de velocidades.txt para hacer un histograma de velocidades.
Copio lo necesario de Voluntario1.c para que funcione.

*/





#include <stdio.h>
#include <string.h>
#include <time.h> //Para la semilla aleatoria
#include <math.h> //Necesario para la potencia
#include <stdlib.h> //Para el uso de rand() y srand()



//Este programa uso el algoritmo de Verlet para simular un gas en dos dimensionas con potencial de Lennard-Jones.

//Ahora defino las constantes:



#define Epsilon 1.0         //Constante de Unidades de potencial
#define Sigma 1.0           //Constante de distancia
#define KB 1.0              //Constante de Boltzmann
#define N 50               //Número de partículas
#define L 10.0              //Longitud de la caja LXL
#define M 1.0               //Masa de las partículas
#define h 0.002              //Paso temporal
#define PI 3.14159265       //Pi
#define T_TOTAL 60        //Tiempo total de simulación

#define a 20 //Valor mínimo del histograma
#define b 50 //Valor máximo del histograma


int numpasos = (int) (T_TOTAL/h) ; //Número de pasos temporales

//ESTE PASO ES NECESARIO PORQUE SINO PETA:
int numparticulas = N; //Número de partículas

double ***crear_arreglo_dinamico(int numparticulas, int numpasos) 
{
    // Asignar memoria para el arreglo dinámico
    double ***arreglo = (double ***)malloc(N * sizeof(double **));

    for (int i = 0; i < numparticulas; i++) {
        arreglo[i] = (double **)malloc(2 * sizeof(double *));
        for (int j = 0; j < 2; j++) {
            arreglo[i][j] = (double *)malloc(numpasos * sizeof(double));
        }
    }

    return arreglo;
}


//Libero la memoria de dichos arrays
void liberar_arreglo_dinamico(double ***arreglo, int numparticulas) {
    // Liberar la memoria asignada al arreglo dinámico
    for (int i = 0; i < numparticulas; i++) {
        for (int j = 0; j < 2; j++) {
            free(arreglo[i][j]);
        }
        free(arreglo[i]);
    }
    free(arreglo);
}


int main (void)
{
    FILE *veltxt = fopen("velocidades.txt", "r"); //Fichero de velocidades
    FILE *cinetica = fopen("cinetica.txt", "r"); //Fichero de energía
    FILE *histograma = fopen("histograma.txt", "w"); // Fichero del histograma
    FILE *maxwelltxt = fopen("maxwell.txt", "w"); // Fichero de la distribución de velocidades de maxwell



    if (cinetica == NULL || veltxt == NULL || histograma == NULL || maxwelltxt == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    double ***vel = crear_arreglo_dinamico(N, numpasos);
    double K[numpasos];
    double T[numpasos];

    // Leer los valores de K desde cinetica.txt
    for (int t = 0; t < numpasos; t++) {
    if (fscanf(cinetica, "%lf", &K[t]) != 1) {
        printf("Error al leer el archivo cinetica.txt en el paso %d.\n", t);
        fclose(cinetica);
        fclose(veltxt);
        return 1;
        }
    }
    // Leer las velocidades desde velocidades.txt
    for (int t = 0; t < numpasos; t++) {
        for (int i = 0; i < N; i++) {
            if (fscanf(veltxt, "%lf, %lf", &vel[i][0][t], &vel[i][1][t]) != 2) {
            printf("Error al leer el archivo velocidades.txt en el paso %d, partícula %d.\n", t, i);
            fclose(cinetica);
            fclose(veltxt);
            return 1;
        }
        }
    // Saltar la línea en blanco que separa los pasos temporales
    fgetc(veltxt);
    }

    fclose(cinetica);
    fclose(veltxt);

    //Vamos ahora con la creación del histograma:
    //Calculo ahora la velocidad media en cada segundo.
    int tiempo = numpasos / T_TOTAL;

    // Controlo los parámetros del histograma:
    int t_min = a * tiempo; // tiempo mínimo
    int t_max = b * tiempo; // tiempo máximo
    double v_media[t_max - t_min]; // velocidad media
    double maxwell[t_max - t_min]; 


    //Calculo la velocidad media y la distribución de velocidades de maxwell:

    for (int t = t_min; t < t_max; t++) {
        if (t >= numpasos) { // Verificar que t no exceda numpasos
            printf("Error: índice t fuera de rango.\n");
            liberar_arreglo_dinamico(vel, N);
            return 1;
        }

        double suma_v = 0.0;
        for (int i = 0; i < N; i++) {
            suma_v += sqrt(vel[i][0][t] * vel[i][0][t] + vel[i][1][t] * vel[i][1][t]);
    }
    v_media[t - t_min] = suma_v / N; // Velocidad media
    maxwell[t - t_min] = sqrt(M/(K[t]))*suma_v/N*exp((-suma_v*suma_v)/(2*K[t]*N*N));

    }

    // Guardo los datos en fichero:
    

    for (int t = t_min; t < t_max; t++) {
        fprintf(histograma, "%lf\n", v_media[t - t_min]);
        fprintf(maxwelltxt, "%lf\n", maxwell[t - t_min]);
    }
    fclose(histograma);
    fclose(maxwelltxt);
    liberar_arreglo_dinamico(vel, N); //Libero la memoria de las velocidades
    
}


