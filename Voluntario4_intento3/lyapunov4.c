
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Defino las constantes:

#define g 9.81
#define PI 3.14159265
#define h 0.001         //paso de tiempo
#define T_Total 10      //tiempo total de simulacion

#define incremento 10   //Incremento de tiempo
#define Tmaximo  200  //Tiempo máximo para lyapunov.

#define E 10.0 //Energía total de cada péndulo.

//Condiciones iniciales:

double thetaini = 0.1;
double phiini = 0.3;

//Diferencias en los ángulos:

double thetadiff = 0.0;
double phidiff = 0.00005;


//El vector que usaré para guardar las coordenadas del péndulo será [theta, phi, momento theta, momento phi], con theta el primer ángulo.

double y[4];
double l[4]; //Coordenadas y momentos del segundo péndulo.

double pos[4]; //Tendrá las coordenadas x e y de ambos pendulos.
double vel[2]; //Contiene las theta' y phi'.
double vell[2]; //Velocidades del segundo péndulo.


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

    k[2][0] = h*dtheta  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][1] = h*dphi  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][2] = h*dmtheta  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][3] = h*dmphi  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);

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


int main(void)
{

    FILE* lyapunovtxt = fopen("lyapunov.txt", "w");

    for (int tiempo = T_Total; tiempo < Tmaximo; tiempo = tiempo + incremento)
    {
        double Lyapunov = 0.0;
        double diferencia_inicial = 0.0;
        double diferencia_total = 0.0;

        // Condiciones iniciales.

        y[0] = thetaini;
        y[1] = phiini;

        l[0] = thetaini - thetadiff;
        l[1] = phiini - phidiff;

        double arg1 = E - 2 * g * (1 - cos(y[0])) - g * (1 - cos(y[1]));
        double arg2 = E - 2 * g * (1 - cos(l[0])) - g * (1 - cos(l[1]));

        y[2] = 2 * sqrt(arg1);
        y[3] = sqrt(arg1) * cos(y[0] - y[1]);

        l[2] = 2 * sqrt(arg2);
        l[3] = sqrt(arg2) * cos(l[0] - l[1]);

        // Calculo las velocidades 

        vel[0] = dtheta(y[0], y[1], y[2], y[3]);
        vel[1] = dphi(y[0], y[1], y[2], y[3]);

        vell[0] = dtheta(l[0], l[1], l[2], l[3]);
        vell[1] = dphi(l[0], l[1], l[2], l[3]);

        diferencia_inicial = sqrt(pow(y[0] - l[0], 2) + pow(y[1] - l[1], 2) +
                                  pow(vel[0] - vell[0], 2) + pow(vel[1] - vell[1], 2));

        double t = 0.0;
        while (t < tiempo)
        {
            rungekutta(y);
            rungekutta(l);

            // Actualizo velocidades:
            vel[0] = dtheta(y[0], y[1], y[2], y[3]);
            vel[1] = dphi(y[0], y[1], y[2], y[3]);

            vell[0] = dtheta(l[0], l[1], l[2], l[3]);
            vell[1] = dphi(l[0], l[1], l[2], l[3]);

            t += h;
        }

        // Diferencia final:
        diferencia_total = sqrt(pow(y[0] - l[0], 2) + pow(y[1] - l[1], 2) +
                                pow(vel[0] - vell[0], 2) + pow(vel[1] - vell[1], 2));

        Lyapunov = log(diferencia_total / diferencia_inicial) / tiempo;

        fprintf(lyapunovtxt, "%lf %d\n", Lyapunov, tiempo);

        printf("He terminado un tiempo: %d\n", tiempo);
        fflush(stdout);
    }

    fclose(lyapunovtxt);

    return 0;
}
