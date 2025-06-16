
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

#define E 1.0 //Energía total de cada péndulo.

//Condiciones iniciales:

double thetaini = 0.1;
double phiini = 0.3;

//Diferencias en los ángulos:

double thetadiff = 0.0;
double phidiff = 0.0000005;


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

    double delta_t = 1.0;
    int pasos = (int)(delta_t / h);
    int N;

    for (int tiempo = T_Total; tiempo < Tmaximo; tiempo += incremento)
    {
        double suma_logs = 0.0;
        N = 0;

        // Condiciones iniciales
        y[0] = thetaini;
        y[1] = phiini;
        l[0] = thetaini - thetadiff;
        l[1] = phiini - phidiff;

        // Cálculo de momentos a partir de energía (deben ser positivos los argumentos)
        double arg1 = E - 2 * g * (1 - cos(y[0])) - g * (1 - cos(y[1]));
        double arg2 = E - 2 * g * (1 - cos(l[0])) - g * (1 - cos(l[1]));

        if (arg1 <= 0 || arg2 <= 0)
        {
            printf("Condiciones iniciales inconsistentes con energía. arg1 = %lf, arg2 = %lf\n", arg1, arg2);
            continue;
        }

        y[2] = 2 * sqrt(arg1);
        y[3] = sqrt(arg1) * cos(y[0] - y[1]);

        l[2] = 2 * sqrt(arg2);
        l[3] = sqrt(arg2) * cos(l[0] - l[1]);

        // Velocidades iniciales
        vel[0] = dtheta(y[0], y[1], y[2], y[3]);
        vel[1] = dphi(y[0], y[1], y[2], y[3]);
        vell[0] = dtheta(l[0], l[1], l[2], l[3]);
        vell[1] = dphi(l[0], l[1], l[2], l[3]);

        double dx0[4] = {
            l[0] - y[0],
            l[1] - y[1],
            vell[0] - vel[0],
            vell[1] - vel[1]
        };

        double d0 = sqrt(dx0[0]*dx0[0] + dx0[1]*dx0[1] + dx0[2]*dx0[2] + dx0[3]*dx0[3]);

        if (d0 == 0.0 || isnan(d0))
        {
            printf("Separación inicial cero o indefinida.\n");
            continue;
        }

        double t = 0.0;
        while (t < tiempo)
        {
            for (int i = 0; i < pasos && t < tiempo; i++)
            {
                rungekutta(y);
                rungekutta(l);
                t += h;
            }

            // Calcular velocidades actualizadas
            vel[0] = dtheta(y[0], y[1], y[2], y[3]);
            vel[1] = dphi(y[0], y[1], y[2], y[3]);
            vell[0] = dtheta(l[0], l[1], l[2], l[3]);
            vell[1] = dphi(l[0], l[1], l[2], l[3]);

            // Calcular nueva separación
            double dx[4] = {
                l[0] - y[0],
                l[1] - y[1],
                vell[0] - vel[0],
                vell[1] - vel[1]
            };

            double d = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3]);

            if (d <= 1e-15 || isnan(d))  // evitar log(0) o NaN
                continue;

            suma_logs += log(d / d0);
            N++;

            // Reescalar diferencia a longitud d0
            double factor = d0 / d;
            for (int j = 0; j < 2; j++) {
                l[j] = y[j] + factor * dx[j];
            }

            // Recalcular momentos desde energía para mantener consistencia
            double arg = E - 2 * g * (1 - cos(l[0])) - g * (1 - cos(l[1]));
            if (arg <= 0)
                break;

            l[2] = 2 * sqrt(arg);
            l[3] = sqrt(arg) * cos(l[0] - l[1]);
        }

        if (N > 0)
        {
            double Lyapunov = suma_logs / (N * delta_t);
            fprintf(lyapunovtxt, "%lf %d\n", Lyapunov, tiempo);
            printf("Tiempo %d terminado. λ ≈ %.5lf\n", tiempo, Lyapunov);
        }
        else
        {
            printf("No se pudo calcular Lyapunov para tiempo = %d\n", tiempo);
        }

        fflush(stdout);
    }

    fclose(lyapunovtxt);
    return 0;
}
