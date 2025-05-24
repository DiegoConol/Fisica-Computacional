

/*

Este es el segundo intento para el pendulo doble. Contiene lo básico.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Defino las constantes:

#define g 9.81
#define PI 3.14159265
#define h 0.001         //paso de tiempo
#define T_Total 60      //tiempo total de simulacion

//Condiciones iniciales:

double thetaini = 0.1;
double phiini = 0.2;


//El vector que usaré para guardar las coordenadas del péndulo será [theta, phi, momento theta, momento phi], con theta el primer ángulo.

double y[4];

double pos[4]; //Tendrá las coordenadas x e y de ambos pendulos.
double vel[2]; //Contiene las theta' y phi'.


// ################### ECUACIONES DE MOVIMIENTO ######################

double dtheta(double theta, double phi, double mtheta, double mphi)
{
    double a = 1/(2-pow(cos(theta-phi),2));
    return a*(mtheta-mphi*cos(theta-phi));

}

double dphi(double theta, double phi, double mtheta, double mphi)
{
    double a = 1/(2-pow(cos(theta-phi),2));
    return a*(2*mphi-mtheta*cos(theta-phi));

}

double dmtheta(double theta, double phi, double mtheta, double mphi)
{
    double a = 1/(2-pow(cos(theta-phi),2));
    double derivadatheta = a*(mtheta-mphi*cos(theta-phi));
    double derivadaphi = a*(2*mphi-mtheta*cos(theta-phi));

    return derivadaphi*derivadatheta*sin(theta-phi)-2*g*sin(theta);

}

double dmphi(double theta, double phi, double mtheta, double mphi)
{
    double a = 1/(2-pow(cos(theta-phi),2));
    double derivadatheta = a*(mtheta-mphi*cos(theta-phi));
    double derivadaphi = a*(2*mphi-mtheta*cos(theta-phi));

    return derivadaphi*derivadatheta*sin(phi-theta)-g*sin(phi);

}



// ################### ECUACIONES DE RUNGE-KUTTA #####################


void rungekutta(double a[4])

{
    //Recuerdo: a[0]=theta, a[1] = phi, a[2] = momento theta, a[3] = momento phi.
    double k[4][4];

    //Calculo k1:

    k[0][0] = h*dtheta  (a[0], a[1], a[2], a[3]);
    k[0][1] = h*dphi    (a[0], a[1], a[2], a[3]);
    k[0][2] = h*dmtheta (a[0], a[1], a[2], a[3]);
    k[0][3] = h*dmphi   (a[0], a[1], a[2], a[3]);

    //Calculo k2:

    k[1][0] = h*dtheta  (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);
    k[1][1] = h*dphi  (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);
    k[1][2] = h*dmtheta  (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);
    k[1][3] = h*dmphi  (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);

    //Calculo k3:

    k[2][0] = h*dtheta  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[2][3]*0.5);
    k[2][1] = h*dphi  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[2][3]*0.5);
    k[2][2] = h*dmtheta  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[2][3]*0.5);
    k[2][3] = h*dmphi  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[2][3]*0.5);

    //Calculo k4:

    k[3][0] = h*dtheta  (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][1] = h*dphi  (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][2] = h*dmtheta  (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][3] = h*dmphi  (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);

    //Actualizo los resultados de a.

    for(int i=0; i<4; i++)
    {
        a[i]+=1.0/6.0*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i]);
    }




}


//Creo ahora el main con los archivos necesarios.

int main(void)
{

    //Primero creo los ficheros, para ello hago un array con las diferentes energías.

    double energias[] = {1.0, 3.0, 5.0, 10.0, 15.0};
    int num_energias = sizeof(energias)/sizeof(energias[0]);

    for (int e=0; e<num_energias; e++)
    {

        
        double E=energias[e];

        //Creo los ficheros para cada energía.
        
        char fname_posiciones[64], fname_angulos[64], fname_momentos[64], fname_fasetheta[64], fname_fasephi[64];

        snprintf(fname_angulos,     sizeof(fname_angulos),      "posiciones/angulos_%.1f.txt", E);
        snprintf(fname_posiciones,  sizeof(fname_posiciones),   "posiciones/posiciones_%.1f.txt",E);
        snprintf(fname_fasetheta,   sizeof(fname_fasetheta),    "espaciofase/fasetheta_%.1f.txt", E);
        snprintf(fname_fasephi,     sizeof(fname_fasephi),      "espaciofase/fasephi_%.1f.txt", E);
        snprintf(fname_momentos,    sizeof(fname_momentos),     "espaciofase/momentos_%.1f.txt", E);


        FILE *angulostxt =      fopen(fname_angulos, "w");
        FILE *posicionestxt =   fopen(fname_posiciones, "w");
        FILE *fasetheta =       fopen(fname_fasetheta, "w");
        FILE *fasephi =         fopen(fname_fasephi, "w");
        FILE *momentostxt =     fopen(fname_momentos, "w");

        //Comprobamos que lo abre:

        if (!angulostxt || !posicionestxt || !fasetheta || !fasephi || !momentostxt) {
        printf("Error al abrir uno de los archivos para E=%.1f\n", E);
        // Cierra los que sí se abrieron
        if (angulostxt) fclose(angulostxt);
        if (posicionestxt) fclose(posicionestxt);
        if (fasetheta) fclose(fasetheta);
        if (fasephi) fclose(fasephi);
        if (momentostxt) fclose(momentostxt);
        continue;
        }
        

        //Condiciones iniciales:

        y[0] = thetaini;
        y[1] = phiini;

        //Comprobamos que se puede para esa energía que phi' = 0.

        double arg = E -2*g*(1-cos(y[0]))-g*(1-cos(y[1]));
        if(arg < 0)
        { 
        printf("El argumento es negativo, no es físicamente posible");
        fclose(angulostxt);
        fclose(posicionestxt);
        fclose(fasetheta);
        fclose(fasephi);
        fclose(momentostxt);
        continue;
        }
        
        y[2] = 2*sqrt(arg);
        y[3] = sqrt(arg)*cos(y[0]-y[1]);

        //Voy a poner al péndulo en una caja 3x3 para que nunca choque con ella.

        pos[0]= 3 + sin(y[0]);
        pos[1]= 3 - cos(y[0]);
        pos[2]= pos[0] + sin(y[1]);
        pos[3]= pos[1] - cos(y[1]);

        //Calculo las velocidades iniciales.

        vel[0] = 1/(2-pow(cos(y[0]-y[1]),2))*(y[2]-y[3]*cos(y[0]-y[1]));
        vel[1] = 1/(2-pow(cos(y[0]-y[1]),2))*(2*y[3]-y[2]*cos(y[0]-y[1]));
       

        //Ponemos las condiciones iniciales en cada fichero:

        fprintf(angulostxt, "%lf %lf\n", y[0], y[1]);
        fprintf(posicionestxt, "%lf %lf %lf %lf\n", pos[0], pos[1], pos[2], pos[3]);
        fprintf(momentostxt, "%lf %lf\n", y[2], y[3]);
        fprintf(fasetheta, "%lf %lf\n", y[0], vel[0]);
        fprintf(fasephi, "%lf %lf\n", y[1], vel[1]);
        


        //Condiciones iniciales listas, vamos al bucle del programa.
        double t=0.0;
        while(t<T_Total)
        {
            rungekutta(y);

            vel[0] = 1/(2-pow(cos(y[0]-y[1]),2))*(y[2]-y[3]*cos(y[0]-y[1]));
            vel[1] = 1/(2-pow(cos(y[0]-y[1]),2))*(2*y[3]-y[2]*cos(y[0]-y[1]));

            pos[0]= 3 + sin(y[0]);
            pos[1]= 3 - cos(y[0]);
            pos[2]= pos[0] + sin(y[1]);
            pos[3]= pos[1] - cos(y[1]);

            fprintf(angulostxt, "%lf %lf\n", y[0], y[1]);
            fprintf(posicionestxt, "%lf %lf %lf %lf\n", pos[0], pos[1], pos[2], pos[3]);
            fprintf(momentostxt, "%lf %lf\n", y[2], y[3]);
            fprintf(fasetheta, "%lf %lf\n", y[0], vel[0]);
            fprintf(fasephi, "%lf %lf\n", y[1], vel[1]);

            t=t+h;
        }
        fclose(angulostxt);
        fclose(posicionestxt);
        fclose(fasetheta);
        fclose(fasephi);
        fclose(momentostxt);

        printf("He terminado una energía\n");
        fflush(stdout);


    }

}