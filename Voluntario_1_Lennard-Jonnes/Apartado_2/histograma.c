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
#define N 20               //Número de partículas
#define L 10.0              //Longitud de la caja LXL
#define M 1.0               //Masa de las partículas
#define h 0.002              //Paso temporal
#define PI 3.14159265       //Pi
#define T_TOTAL 60        //Tiempo total de simulación

#define a 20 //Valor mínimo del histograma
#define b 50 //Valor máximo del histograma
#define NUM_BINS 30
#define V_MAX 20.0 // Velocidad máxima para el histograma


int numpasos = (int) (T_TOTAL/h) ; //Número de pasos temporales
int t_min_pasos = a/h; //Número de pasos temporales mínimo
int t_max_pasos = b/h; //Número de pasos temporales máximo


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


double maxwell_distribucion(double v, double T)
{
    return (M*v/(KB*T))*exp((-M*v*v)/(2*KB*T));
}

int main (void)
{
    FILE *veltxt = fopen("velocidades.txt", "r"); //Fichero de velocidades
    FILE *cinetica = fopen("cinetica.txt", "r"); //Fichero de energía
    FILE *histograma = fopen("histograma.txt", "w"); // Fichero del histograma



    if (cinetica == NULL || veltxt == NULL || histograma == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    double ***vel = crear_arreglo_dinamico(N, numpasos);
    double K[numpasos];

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


    //Calculo la velocidad media para saber la temperatura media: NO SE HACE PARA EL HISTOGRAMA, SINO PARA LA TEMPERATURA MEDIA:
    double sumatotalcuadr=0.0;

    for (int t = t_min; t < t_max; t++) {
        if (t >= numpasos) { // Verificar que t no exceda numpasos
            printf("Error: índice t fuera de rango.\n");
            liberar_arreglo_dinamico(vel, N);
            return 1;
        }
        
        double suma_v = 0.0;
        for (int i = 0; i < N; i++) {
            suma_v += sqrt(vel[i][0][t] * vel[i][0][t] + vel[i][1][t] * vel[i][1][t]);
            sumatotalcuadr += (vel[i][0][t] * vel[i][0][t] + vel[i][1][t] * vel[i][1][t]);
        }
        
        v_media[t - t_min] = suma_v / N; // Velocidad media
    }

    //Para sacar la v_media total cojo sumatotal y lo divido entre el número de pasos y N:
    double v_media_total = sumatotalcuadr /(t_max - t_min);

    //Con esto calculo la temperatura media:
    double T_media = M * v_media_total/ (2*N*KB);
    printf("La temperatura media es: %lf\n", T_media);



    //Ahora calculo la velocidad de verdad para el histograma:

    double v_histograma[N];
    double v_histograma_x[N];
    double v_histograma_y[N];

    for(int t=t_min; t<t_max; t++)
    {
        for (int i=0; i<N; i++)
        {
            v_histograma[i] = sqrt(vel[i][0][t]*vel[i][0][t] + vel[i][1][t]*vel[i][1][t]);
            if(vel[i][0][t]>0.0)
            v_histograma_x[i] = vel[i][0][t];
            else
            v_histograma_x[i] = -vel[i][0][t];

            if(vel[i][1][t]>0.0)
            v_histograma_y[i] = vel[i][1][t];
            else
            v_histograma_y[i] = -vel[i][1][t];

            fprintf(histograma, "%lf %lf %lf\n", v_histograma[i], v_histograma_x[i], v_histograma_y[i]); // Guardar las velocidades en el archivo

            v_histograma[i] = 0.0; // Reiniciar el valor para la siguiente iteración (me ahorro lo de la memoria pero es más ineficiente)
        }
    }


    // Guardo los datos en fichero:
    

    for (int t = 0; t < N; t++) {
        fprintf(histograma, "%lf\n", v_histograma[t]);
    }
    fclose(histograma);
    liberar_arreglo_dinamico(vel, N); //Libero la memoria de las velocidades
    
}


