
/*

Este es un programa que hace un pendulo doble y exporta los datos de los angulos, angulo's y hamiltoniano a un archivo.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// DEFINIMOS LAS CONSTANTES

#define g 9.81
#define PI 3.14159265
#define T_TOTAL 60 //tiempo total.
#define h 0.01 //paso temporal

//Los parámetros del pendulo (según el voluntario son =1 para simplifcar el problema)

#define m1 1.0 // masa del primer pendulo
#define m2 1.0 // masa del segundo pendulo
#define l1 1.0 // longitud del primer pendulo
#define l2 1.0 // longitud del segundo pendulo




//Creo la función que calcula el hamiltoniano.

double hamiltoniano(double theta, double psi, double dtheta, double dpsi, FILE *file)
{

    double T = 1/2*dpsi*dpsi + dtheta*dtheta +dtheta*dpsi*cos(psi-theta);
    double V = 2*g*(1-cos(theta))+g*(1-cos(psi));
    double H = T+V;

    fprintf(file, "%lf\n", H);

    return H;

}

//Como es un sistema aislado y no hay ligaduras reónomas la energía se conserva, 
//y como T es función cuadrática de las velocidades y V no depende de estas entonces H=E.


//Resolvemos las ecuaciones diferenciales para que nos den soluciones separadas para theta y psi.




//Ahora hago el algoritmo de Runge-Kutta de cuarto orden.

void rungekutta (double theta, double psi, double dtheta, double dpsi, FILE *file)
{
    
}





int main(void)
{

    FILE *posiciones=("posiciones.txt");
    FILE *velocidades=("velocidades.txt");
    FILE *hamiltoniano=("hamiltoniano.txt");
    FILE *pendulo=("pendulo.txt");

    if (posiciones == NULL || velocidades == NULL || hamiltoniano == NULL || pendulo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }



    fclose(posiciones);
    fclose(velocidades);
    fclose(hamiltoniano);
    fclose(pendulo);


}