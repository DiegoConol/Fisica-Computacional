#include <stdio.h>
#include <string.h>
#include <math.h> //Necesario para la potencia
#include <stdlib.h> //Para el uso de rand() y srand()



//Este programa uso el algoritmo de Verlet para simular un gas en dos dimensionas con potencial de Lennard-Jones.

//Ahora defino las constantes:



#define Epsilon 1.0      //Constante de Unidades de potencial
#define Sigma 1.0        //Constante de distancia
#define KB 1.0           //Constante de Boltzmann
#define N 100            //Número de partículas
#define L 10.0           //Longitud de la caja LXL
#define T 1.0            //Temperatura
#define M 1.0            //Masa de las partículas
#define h 0.01           //Paso temporal



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




void aceleracion(double r[N][2], double a[N][2])
{

    for (int i=0; i<N; i++)
    {
    for (int j=0; j<N; j++)
        {
            if (i != j) //No calculamos la aceleración de una partícula sobre sí misma.
            {
                double r = pow(dr[i][j][0]*dr[i][j][0] + dr[i][j][1]*dr[i][j][1], 0.5); //Distancia al cuadrado
                double acc = 4*Epsilon*(12*pow(Sigma,12)/pow(r,13)-6*pow(Sigma,6)/pow(r,7))/(M*r); //Aceleración

                //Para evitar muchos cálculos, lo guardamos en un vector auxiliar.
                a[i][0] += acc*dr[i][j][0]; //Aceleración en x
                a[i][1] += acc*dr[i][j][1]; //Aceleración en y

            }
    }
}

