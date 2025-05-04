#include <stdio.h>
#include <string.h>
#include <math.h> //Necesario para la potencia
#include <stdlib.h> //Para el uso de rand() y srand()



//Este programa uso el algoritmo de Verlet para simular un gas en dos dimensionas con potencial de Lennard-Jones.

//Ahora defino las constantes:



#define Epsilon 1.0         //Constante de Unidades de potencial
#define Sigma 1.0           //Constante de distancia
#define KB 1.0              //Constante de Boltzmann
#define N 20               //Número de partículas
#define L 10.0              //Longitud de la caja LXL
#define T 1.0               //Temperatura
#define M 1.0               //Masa de las partículas
#define h 0.01              //Paso temporal
#define PI 3.14159265       //Pi
#define T_TOTAL 10.0        //Tiempo total de simulación




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

                //Para evitar muchos cálculos, lo guardo en un vector auxiliar.
                a[i][0] += acc*dr[i][j][0]; //Aceleración en x
                a[i][1] += acc*dr[i][j][1]; //Aceleración en y

            }
    }
}
}

//Ahora hago el algoritmo de Verlet.

void verlet(double r[N][2], double v[N][2], double a[N][2], FILE *file)
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

    //Actualizo las aceleraciones (t+h)
    aceleracion(r,a);

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
        fprintf(file, "%lf, %lf, %lf, %lf\n", r[i][0], r[i][1], v[i][0], v[i][1]);
    }
    fprintf(file, "\n"); //Salto de línea para separar los pasos
}


//Ahora hago la función principal.

int main(void)
{

    FILE *salida = fopen("SALIDA.txt", "w");    //Fichero de salida
    if (salida == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    //Las posiciones iniciales y velocidad son aleatorias. Velocidades con módulo 1 y direcciones aleatorias.
    srand(0); //Semilla para la aleatoriedad

    //Las posiciones tienen que estar en la caja LXL, así que las inicializo aleatoriamente.
    
    for (int i=0; i<N; i++)
    {
        r[i][0] = ((double) rand() / (double) RAND_MAX)*L;  //Posición en x
        r[i][1] = ((double) rand() / (double) RAND_MAX)*L;  //Posición en y

        //Para la dirección cojo un ángulo aleatorio entre 0 y 2pi y lo paso a coordenadas cartesianas.

        double theta = ((double) rand() / (double) RAND_MAX)*2*PI; //Dirección aleatoria
        v[i][0] = cos(theta);   //Velocidad en x
        v[i][1] = sin(theta);   //Velocidad en y

        //Inicializo la aceleración en 0.0.
        a[i][0] = 0.0;  //Aceleración en x
        a[i][1] = 0.0;  //Aceleración en y
    }



    //Ahora hago el bucle de la simulación. El tiempo total es T_TOTAL y el paso temporal es h.
    for (double t=0; t<T_TOTAL; t+=h)
    {
        //Calculo la distancia entre partículas.
        distancia(r, dr);

        //Calculo la aceleración de cada partícula.
        aceleracion(r, a);

        //Hago el algoritmo de Verlet.
        verlet(r, v, a, salida);
    }

}
