
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
#define E 10.0 //Energía total del sistema.

//Creo el vector que tendrá las coordenadas: [theta, phi, momento de theta, momento de phi]

double y[4];


//Creo la función que calcula el hamiltoniano.

void hamiltoniano(double vector[4], FILE *file)
{
    double H=0.0;
    double aux1 = 2-cos(y[0]-y[1])*cos(y[0]-y[1]);
    double auxt = y[2]-y[3]*cos(y[0]-y[1]);
    double auxp = 2*y[3]-y[2]*cos(y[0]-y[1]);

    H = 1.0/aux1/aux1*(auxt*auxt+1/2*auxp*auxp+auxt*auxp*cos(y[0]-y[1]))+g*(3-cos(y[0])-cos(y[1]));

    fprintf(file, "%lf\n", H);

}

//Como es un sistema aislado y no hay ligaduras reónomas la energía se conserva, 
//y como T es función cuadrática de las velocidades y V no depende de estas entonces H=E.


//A continuación resuelvo el Hamiltoniano, obteniendo las ecuaciones de las derivadas de los momentos y las velocidades.


// ####################### ECUACIONES DE MOVIMIENTO ######################


double dtheta (double theta, double phi, double mtheta, double mphi)
{
    double aux=0.0;
    aux = 1/(2-cos(theta - phi)*cos(theta-phi))*(mtheta-mphi);
    return aux;
}

double dphi(double theta, double phi, double mtheta, double mphi)
{

    double aux=0.0;
    aux = 1/(2-cos(theta - phi)*cos(theta-phi))*(2*mphi-mtheta*cos(theta-phi));
    return aux;
}

double dmtheta(double theta, double phi, double mtheta, double mphi)
{
    double aux=0.0;
    double c = cos(theta-phi);
    aux=sin(theta-phi)*(1/(c*c))*(mtheta*mphi*(2+c*c)-(mtheta*mtheta*c+2*mphi*mphi*c)) - 2*g*sin(theta);
    return aux;

}

double dmphi(double theta, double phi, double mtheta, double mphi)
{

    double aux=0.0;
    double c = cos(theta-phi);
    aux=sin(phi-theta)*(1/(c*c))*(mtheta*mphi*(2+c*c)-(mtheta*mtheta*c+2*mphi*mphi*c)) - g*sin(theta);
    return aux;

}


/* ###################### RUNGE KUTTA ###################

//Resolvemos las ecuaciones diferenciales para que nos den soluciones separadas para theta y psi.


//Ahora hago el algoritmo de Runge-Kutta de cuarto orden.
//Resolviendo el hamiltoniano tenemos las funciones p'(theta, psi) y las velocidades. 4 ecuaciones de 4 incógnitas
*/

void rungekutta (double vector[4])
{

    double k[4][4];

    //Primero evaluamos k1, que es el paso*funcion.

    k[0][0]=h*dtheta(vector[0], vector[1], vector[2], vector[3]);
    k[0][1]=h*dphi(vector[0], vector[1], vector[2], vector[3]);
    k[0][2]=h*dmtheta(vector[0], vector[1], vector[2], vector[3]);
    k[0][3]=h*dmphi(vector[0], vector[1], vector[2], vector[3]);

    //Ahora calculo k2.

    k[1][0]=h*dtheta(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][1]=h*dphi(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][2]=h*dmtheta(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][3]=h*dmphi(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);

    //Ahora calculo k3.


    k[2][0]=h*dtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][1]=h*dphi(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][2]=h*dmtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][3]=h*dmphi(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);

    //Por último, k4.

    k[3][0]=h*dtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][1]=h*dphi(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][2]=h*dmtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][3]=h*dmphi(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);

    //actualizamos el vector inicial

    for(int i=0; i<4; ++i)
    {
        vector[i]+= 1.0/6.0*(k[i][0]+2*k[i][1]+2*k[i][2]+k[i][3]);

    }
}





int main(void)
{

    FILE *angulostxt = fopen("angulos.txt", "w"); 
    FILE *momentostxt = fopen("momentos.txt", "w");
    FILE *hamiltonianotxt = fopen("hamiltoniano.txt", "w");
    FILE *pendulotxt = fopen("pendulo.txt", "w");

    if (angulostxt == NULL || momentostxt == NULL || hamiltonianotxt == NULL || pendulotxt == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    //Establecemos las condiciones iniciales. Tenemos 4 parámetros libres pero lo reducimos a 2: los ángulos iniciales.
    // Para ello asumimos phi'=0, por tanto calculamos su momento -> obtenemos mphi.
    //También tenemos que la energía tiene un valor constante, que usamos para sacar el otro valor -> obtenemos mtheta.

    y[0]=0.3; //Ángulo theta.
    y[1]=0.4; //Ángulo phi.
    y[2]=2*sqrt(E-2*g*(1-cos(y[0])-(1-cos(y[1])))); //Momento theta.
    y[3]=y[2]/2*cos(y[0]-y[1]); //Momento phi.

    //Hago el bucle principal.

    int numpasos = T_TOTAL/h;

    for(int i=0; i<numpasos; i++)
    {
        rungekutta(y);
        hamiltoniano(y, hamiltonianotxt);
        fprintf(angulostxt, "%lf, %lf\n", y[0], y[1]);
        fprintf(momentostxt, "%lf, %lf\n", y[2], y[3]);

    }



    fclose(angulostxt);
    fclose(momentostxt);
    fclose(hamiltonianotxt);
    fclose(pendulotxt);


}