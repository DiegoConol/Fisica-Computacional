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
#define mod 1 //Modulo de la velocidad







//Primero a qué posición y velocidad inicial tienen las partículas:
//N es las particulas, 2 es para las coordenadas x e y.

double r[N][2]; //Posición de las partículas
double v[N][2]; //Velocidad de las partículas
double a[N][2]; //Aceleración de las partículas

//Ahora hago el vector diferencia de posiciones:

double dr[N][N][2]; //Diferencia de posiciones entre partículas i y j.


//Voy a empezar por la condición de periodicidad. Estamos en una caja de LXL.

void periodicidad(double r[N][2])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<2; j++)
        {
            if (r[i][j] > L)
            {
                r[i][j] -= L;
            }
            else if (r[i][j] < 0)
            {
                r[i][j] += L;
            }
        }
    }
}


//Ahora hago las distancias entre partículas.
/*
    El vector: dr[N][N][2] es la diferencia de posiciones entre partículas i y j.
    donde dr[1][2][0] = x1-x2 y dr[1][2][1] = y1-y2 por ejemplo

*/

void distancia(double r[N][2], double dr[N][N][2])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            if (i != j) //No hacemos de una a si misma y usamos periodicidad
            {
                double a=r[i][0]-r[j][0]; //Diferencia en x
                double b=r[i][1]-r[j][1]; //Diferencia en y

                if (a> L/2)
                dr[i][j][0] = a - L; 
                else if (a < -L/2)
                dr[i][j][0] = a + L; 
                else
                dr[i][j][0] = a; 

                if (b> L/2)
                dr[i][j][1] = b - L; 
                else if (b < -L/2)
                dr[i][j][1] = b + L; 
                else
                dr[i][j][1] = b;

            }
        }
    }
}


//Ahora implemento la aceleración. F=m*a.


void aceleracion(double dr[N][N][2], double a[N][2])
{

    //Inicializo la aceleración a 0.
    for (int i=0; i<N; i++)
    {
        a[i][0] = 0.0; //Aceleración en x
        a[i][1] = 0.0; //Aceleración en y
    }

    for (int i=0; i<N; i++)
    {
    for (int j=0; j<N; j++)
        {
            if (i != j) //No calculamos la aceleración de una partícula sobre sí misma.
            {
                double r = pow(dr[i][j][0]*dr[i][j][0] + dr[i][j][1]*dr[i][j][1], 0.5); //Distancia
                if (r < 1e-10) {
                    continue; // Evita calcular la aceleración si la distancia es muy pequeña
                }
                double acc = 4*Epsilon*(12*pow(Sigma,12)/pow(r,13)-6*pow(Sigma,6)/pow(r,7))/(M*r); //Aceleración

                //Para evitar muchos cálculos, lo guardo en un vector auxiliar.
                a[i][0] += acc*dr[i][j][0]; //Aceleración en x
                a[i][1] += acc*dr[i][j][1]; //Aceleración en y
                a[j][0] -= acc*dr[i][j][0]; //Aceleración en x
                a[j][1] -= acc*dr[i][j][1]; //Aceleración en y

            }
    }
}
}

//Ahora hago el algoritmo de Verlet.

void verlet(double r[N][2], double v[N][2], double a[N][2], double dr[N][N][2], FILE *salir, FILE *pos, FILE *vel, FILE *acel)

{
    double omega[N][2]; //Vector auxiliar para la velocidad

    //Actualizo posicicones (t+h) y calculo omega con velocidad en t
    for (int i=0; i<N; i++)
    {
        r[i][0] += v[i][0]*h + 0.5*a[i][0]*h*h;     //Posición en x
        r[i][1] += v[i][1]*h + 0.5*a[i][1]*h*h;     //Posición en y
        omega[i][0] = v[i][0] + h/2*a[i][0];        //Velocidad en x
        omega[i][1] = v[i][1] + h/2*a[i][1];        //Velocidad en y
    }

    //Actualizo la periodicidad de las posiciones
    periodicidad(r);

    //Actualizo la distancia entre partículas (t+h)
    distancia(r, dr);

    //Actualizo las aceleraciones (t+h)
    aceleracion(dr,a);

    //Actualizo las velocidades (t+h)
    for (int i=0; i<N; i++)
    {
        v[i][0] = omega[i][0] + h/2*a[i][0];        //Velocidad en x
        v[i][1] = omega[i][1] + h/2*a[i][1];        //Velocidad en y
    }

    //Pongo los nuevos parámetros r,v en el fichero.
    /*
    Se guardan así:
    rx1, ry1, vx1, vy1
    rx2, ry2, vx2, vy2
    ...
    rxN, ryN, vxN, vyN
    "Salto de línea para separar los pasos"

    rx1, ry1, vx1, vy1
    rx2, ry2, vx2, vy2
    ...
    rxN, ryN, vxN, vyN
    */

    for (int i=0; i<N; i++)
    {
        fprintf(salir, "%lf, %lf, %lf, %lf, %lf, %lf\n", r[i][0], r[i][1], v[i][0], v[i][1], a[i][0], a[i][1]);
        fprintf(pos, "%lf, %lf\n", r[i][0], r[i][1]); //Posiciones
        fprintf(vel, "%lf, %lf\n", v[i][0], v[i][1]); //Velocidades
        fprintf(acel, "%lf, %lf\n", a[i][0], a[i][1]); //Aceleraciones
    }
    //Salto de línea para separar los pasos
    fprintf(salir, "\n"); 
    fprintf(pos, "\n"); 
    fprintf(vel, "\n"); 
    fprintf(acel, "\n"); 
}



//Ahora hago la energía potencial y cinética. La energía total es la suma de ambas.


double energia(double dr[N][N][2], double v[N][2], FILE *file)
{
    double U = 0.0;          //Energía potencial
    double K = 0.0;          //Energía cinética
    double E = 0.0;         //Energía total

    for (int i=0; i<N; i++)
    {
        K += 0.5*M*(v[i][0]*v[i][0] + v[i][1]*v[i][1]); //Energía cinética
        for (int j=0; j<N; j++)
        {
            if (i != j) //No calculamos la energía potencial de una partícula sobre sí misma.
            {
                double r = pow(dr[i][j][0]*dr[i][j][0] + dr[i][j][1]*dr[i][j][1], 0.5); //Distancia
                if (r < 1e-10) {
                    continue; 
                }
                U += 4*Epsilon*(pow(Sigma,12)/pow(r,12)-pow(Sigma,6)/pow(r,6)); //Energía potencial
            }
        }
    }

    E= K + U; //Energía total


    fprintf(file, "%lf, %lf, %lf\n",K, U, E); //Escribo la energía total en el fichero
    return K; //Devuelvo la energía cinética para luego ussarla en K = kbT y calcular T.
}



//Ahora hago la función de guardar las velocidades.

/* ################### VOY A USAR MEMORIA DINÁMICA PARA GUARDAR LAS VELOCIDADES Y POSICIONES EN CADA PASO ###################


Es una función muy grande.
Crea memoria para un array tridimensional de tamaño N x 2 x numpasos.
El primer índice es el número de partículas, el segundo índice es 2 (x e y) y el tercer índice es el número de pasos temporales.
*/

//Ahora hago la función de guardar las velocidades y posiciones con memoria dinámica

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





//Ahora hago la función principal.


int main(void)
{

    FILE *salida = fopen("SALIDA.txt", "w");        //Fichero de salida
    FILE *postxt = fopen("posiciones.txt", "w"); //Fichero de posiciones
    FILE *veltxt = fopen("velocidades.txt", "w"); //Fichero de velocidades
    FILE *aceltxt = fopen("aceleraciones.txt", "w"); //Fichero de aceleraciones
    FILE *energiatxt = fopen("energia.txt", "w"); //Fichero de energía
    FILE *cinetica = fopen("cinetica.txt", "w"); //Fichero de energía cinética


    if (salida == NULL || energiatxt == NULL || postxt == NULL || veltxt == NULL || aceltxt == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    //Las posiciones iniciales y velocidad son aleatorias. Velocidades con módulo 1 y direcciones aleatorias.
    
    //Pongo la semilla aleatoria para que cada vez que ejecute el programa me de un resultado diferente.
    srand(time(NULL));

    //Las posiciones tienen que estar en la caja LXL, así que las inicializo aleatoriamente.
    
    for (int i=0; i<N; i++)
    {
        //Voy a poner las posiciones como si estuvieran en una cuadricula, así estarán distribuidas uniformemente.
        //El poner el N+1 es para que no ocupe los bordes de la caja.

        r[i][0] = (i % (int)sqrt(N) + 1) * (L / ((int)sqrt(N) + 1)) + 0.1*(rand()/(double)RAND_MAX);
        r[i][1] = (i / (int)sqrt(N) + 1) * (L / ((int)sqrt(N) + 1)) + 0.1*(rand()/(double)RAND_MAX);

        //Para la dirección cojo un ángulo aleatorio entre 0 y 2pi y lo paso a coordenadas cartesianas.

        //PONGO UNA VELOCIDAD DE MODULO 4 por ejemplo

        double theta = ((double) rand() / (double) RAND_MAX); //Dirección aleatoria
        v[i][0] = theta*mod;   //Velocidad en x
        v[i][1] = 0.0;   //Velocidad en y

        
        
    }
    //Inicializo la distancia entre partículas y la aceleración.
    //También calculo la energía inicial.

    double ***vel = crear_arreglo_dinamico(N, numpasos);
    double ***pos = crear_arreglo_dinamico(N, numpasos);

    double K[numpasos];
    double T[numpasos];  

    distancia(r, dr); 

    aceleracion(dr,a);

    K[0]= energia(dr, v, energiatxt);
    T[0]= K[0]/(KB);


    for(int i=0; i<N; i++)
    {
        //Imprimo las posiciones y velocidades iniciales.
        fprintf(salida, "%lf, %lf, %lf, %lf, %lf, %lf\n", r[i][0], r[i][1], v[i][0], v[i][1], a[i][0], a[i][1]);
        fprintf(postxt, "%lf, %lf\n", r[i][0], r[i][1]); //Posiciones
        fprintf(veltxt, "%lf, %lf\n", v[i][0], v[i][1]); //Velocidades
        fprintf(aceltxt, "%lf, %lf\n", a[i][0], a[i][1]); //Aceleraciones
        fprintf(cinetica, "%lf\n", K[0]); //Energía cinética inicial
    }
    //Salto de línea para separar los pasos
    fprintf(salida, "\n");
    fprintf(postxt, "\n");
    fprintf(veltxt, "\n");
    fprintf(aceltxt, "\n");
    fprintf(energiatxt, "\n"); 


    //Ahora hago el bucle de la simulación. El tiempo total es T_TOTAL y el paso temporal es h.
    for (int t=0; t<numpasos; t++)
    {
        if (t % (numpasos / 100) == 0) {
            printf("Progreso: %d%% completado\n", (t * 100) / numpasos);
            fflush(stdout); // Con esto no se raya, lo suelta por pantalla siempre.
        }
        
        verlet(r, v, a, dr, salida, postxt, veltxt, aceltxt);
        K[t+1]=energia(dr, v, energiatxt);
        T[t+1]= K[t+1]/(KB);

        //Actualizo las posiciones y velocidades en el array tridimensional.
        for (int i=0; i<N; i++)
        {
            vel[i][0][t] = v[i][0]; //Velocidad en x
            vel[i][1][t] = v[i][1]; //Velocidad en y
            pos[i][0][t] = r[i][0]; //Posición en x
            pos[i][1][t] = r[i][1]; //Posición en y
        }
        fprintf(cinetica, "%lf\n", K[t+1]); //Energía cinética
    }



    //Cierro los ficheros.
    fclose(salida); 
    fclose(postxt); 
    fclose(veltxt); 
    fclose(aceltxt); 
    fclose(energiatxt);
    fclose(cinetica);
    
    liberar_arreglo_dinamico(vel, N); //Libero la memoria de las velocidades
    liberar_arreglo_dinamico(pos, N); //Libero la memoria de las posiciones
    return 0;

}
