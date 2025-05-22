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
#define T_TOTAL 30 //tiempo total. (No recomendable poner más de 30, sino se raya. Poner 30 o bajar el paso.)
#define h 0.001 //paso temporal

//Los parámetros del pendulo (según el voluntario son =1 para simplifcar el problema)

#define m1 1.0 // masa del primer pendulo
#define m2 1.0 // masa del segundo pendulo
#define l1 1.0 // longitud del primer pendulo
#define l2 1.0 // longitud del segundo pendulo
double E = 0.0;

#define Energiamax 15 //Energía máxima que alcanza el sistema en el bucle.

//Creo el vector que tendrá las coordenadas: [theta, phi, momento de theta, momento de phi]

double y[4];


//Creo la función que calcula el hamiltoniano.

void hamiltoniano(double y[4], FILE *file)
{

    double H=0.0;
    double c = cos(y[0]-y[1]);
    double dtheta1 = 1/(2-c*c)*(y[2]-y[3]*c);
    double dphi1 = 1/(2-c*c)*(2*y[3]-y[2]*c);

    H= dtheta1*dtheta1 +0.5*dphi1*dphi1 + dtheta1*dphi1*c+2*g*(1-cos(y[0]))+g*(1-cos(y[1]));
    fprintf(file, "%lf\n", H);

}

//Como es un sistema aislado y no hay ligaduras reónomas la energía se conserva, 
//y como T es función cuadrática de las velocidades y V no depende de estas entonces H=E.


//A continuación resuelvo el Hamiltoniano, obteniendo las ecuaciones de las derivadas de los momentos y las velocidades.


// ####################### ECUACIONES DE MOVIMIENTO ######################

double dtheta (double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta-phi);
    double d = 2.0 - c*c; 
    double aux= 1.0/d*(mtheta-mphi*c);
    return aux;
}

double dphi(double theta, double phi, double mtheta, double mphi)
{

    double c = cos(theta-phi);
    double d = 2.0 - c*c;
    double aux= 1.0/d*(2.0*mphi-mtheta*c);
    return aux;
}

double dmtheta(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    double aux=s/(d*d)*(2.0*mtheta*mphi + c*c*mtheta*mphi - c*(mtheta*mtheta + 2.0*mphi*mphi))-2.0*g*sin(theta);
    return aux;

}

double dmphi(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    double aux = -s/(d*d)*(2.0*mtheta*mphi + c*c*mtheta*mphi - c*(mtheta*mtheta + 2.0*mphi*mphi))-g*sin(phi);
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
    //El primer indice indica la ecuación y el segundo el k.

    
    //Calculamos k1
    k[0][0] = h*dtheta(vector[0], vector[1], vector[2], vector[3]);
    k[1][0] = h*dphi(vector[0], vector[1], vector[2], vector[3]);
    k[2][0] = h*dmtheta(vector[0], vector[1], vector[2], vector[3]);
    k[3][0] = h*dmphi(vector[0], vector[1], vector[2], vector[3]);

    //Calculamos k2
    k[0][1] = h*dtheta(vector[0]+k[0][0]*0.5, vector[1]+k[1][0]*0.5, vector[2]+k[2][0]*0.5, vector[3]+k[3][0]*0.5);
    k[1][1] = h*dphi(vector[0]+k[0][0]*0.5, vector[1]+k[1][0]*0.5, vector[2]+k[2][0]*0.5, vector[3]+k[3][0]*0.5);
    k[2][1] = h*dmtheta(vector[0]+k[0][0]*0.5, vector[1]+k[1][0]*0.5, vector[2]+k[2][0]*0.5, vector[3]+k[3][0]*0.5);
    k[3][1] = h*dmphi(vector[0]+k[0][0]*0.5, vector[1]+k[1][0]*0.5, vector[2]+k[2][0]*0.5, vector[3]+k[3][0]*0.5);

    //Calculamos k3
    k[0][2] = h*dtheta(vector[0]+k[0][1]*0.5, vector[1]+k[1][1]*0.5, vector[2]+k[2][1]*0.5, vector[3]+k[3][1]*0.5);
    k[1][2] = h*dphi(vector[0]+k[0][1]*0.5, vector[1]+k[1][1]*0.5, vector[2]+k[2][1]*0.5, vector[3]+k[3][1]*0.5);
    k[2][2] = h*dmtheta(vector[0]+k[0][1]*0.5, vector[1]+k[1][1]*0.5, vector[2]+k[2][1]*0.5, vector[3]+k[3][1]*0.5);
    k[3][2] = h*dmphi(vector[0]+k[0][1]*0.5, vector[1]+k[1][1]*0.5, vector[2]+k[2][1]*0.5, vector[3]+k[3][1]*0.5);

    //Calculamos k4
    k[0][3] = h*dtheta(vector[0]+k[0][2], vector[1]+k[1][2], vector[2]+k[2][2], vector[3]+k[3][2]);
    k[1][3] = h*dphi(vector[0]+k[0][2], vector[1]+k[1][2], vector[2]+k[2][2], vector[3]+k[3][2]);
    k[2][3] = h*dmtheta(vector[0]+k[0][2], vector[1]+k[1][2], vector[2]+k[2][2], vector[3]+k[3][2]);
    k[3][3] = h*dmphi(vector[0]+k[0][2], vector[1]+k[1][2], vector[2]+k[2][2], vector[3]+k[3][2]);

    for(int i=0; i<4; ++i)
    {
        vector[i]+= 1.0/6.0*(k[i][0]+2*k[i][1]+2*k[i][2]+k[i][3]);

    }
}





int main(void)
{
    FILE *pendulotxt = fopen("pendulo.txt", "w");
    FILE *angulosenergiastxt = fopen("angulosenergias.txt", "w");
    FILE *energiasdistintastxt = fopen("energiasdistintas.txt", "w");

    //El vector pos[4] indica las posiciones en el eje X,Y de la particula 1 y la particula 2. siendo [x1, y1, x2, y2].
    double pos[4];
    int numpasos = T_TOTAL/h;

    // Condiciones iniciales
    double energias[] = {1.0, 3.0, 5.0, 10.0, 15.0};
    int num_energias = sizeof(energias)/sizeof(energias[0]);

    for (int e = 0; e < num_energias; ++e) {
        double E = energias[e];
        char fname_angulos[64], fname_momentos[64], fname_hamiltoniano[64], fname_posiciones[64];

        snprintf(fname_angulos, sizeof(fname_angulos), "angulos_%.1f.txt", E);
        snprintf(fname_momentos, sizeof(fname_momentos), "momentos_%.1f.txt", E);
        snprintf(fname_hamiltoniano, sizeof(fname_hamiltoniano), "hamiltoniano_%.1f.txt", E);
        snprintf(fname_posiciones, sizeof(fname_posiciones), "posiciones_%.1f.txt", E);

        FILE *angulostxt = fopen(fname_angulos, "w");
        FILE *momentostxt = fopen(fname_momentos, "w");
        FILE *hamiltonianotxt = fopen(fname_hamiltoniano, "w");
        FILE *posicionestxt = fopen(fname_posiciones, "w");

        // Condiciones iniciales
        y[0]=0.1; // Ángulo theta
        y[1]=0.1;  // Ángulo phi

        double arg = E - 2*g*(1-cos(y[0])) - g*(1-cos(y[1]));
        if (arg < 0) {
            printf("El argumento es negativo para E=%.1f, no podemos hacer eso\n", E);
            fclose(angulostxt);
            fclose(momentostxt);
            fclose(hamiltonianotxt);
            fclose(posicionestxt);
            continue;
        }
        y[2]=2*sqrt(arg); // Momento theta
        y[3]=y[2]*0.5*cos(y[0]-y[1]); // Momento phi

        for(int i=0; i<numpasos; i++) {
            rungekutta(y);
            hamiltoniano(y, hamiltonianotxt);
            fprintf(angulostxt, "%lf %lf\n", y[0], y[1]);
            fprintf(momentostxt, "%lf %lf\n", y[2], y[3]);

            // Calculo las posiciones para la animación
            pos[0]= 3 + sin(y[0]);
            pos[1]= 3 - cos(y[0]);
            pos[2]= pos[0] + sin(y[1]);
            pos[3]= pos[1] - cos(y[1]);
            fprintf(posicionestxt, "%lf %lf %lf %lf\n", pos[0], pos[1], pos[2], pos[3]);
        }

        fclose(angulostxt);
        fclose(momentostxt);
        fclose(hamiltonianotxt);
        fclose(posicionestxt);
    }

    fclose(pendulotxt);
    fclose(angulosenergiastxt);
    fclose(energiasdistintastxt);

    return 0;
}