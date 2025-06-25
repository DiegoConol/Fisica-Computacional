/*
Este es el tercer intento para el pendulo doble. Calculamos lyapunov
CORREGIDO: Algoritmo de Lyapunov con renormalización y Runge-Kutta arreglado
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Defino las constantes:
#define g 9.81
#define PI 3.14159265
#define h 0.001         //paso de tiempo
#define T_Total 0     //tiempo total de simulacion

#define incremento 10   //Incremento de tiempo
#define Tmaximo  100  //Tiempo máximo para lyapunov.

#define E 10.0 //Energía total de cada péndulo.

//Condiciones iniciales:
double thetaini = 0.1;
double phiini = 0.3;

//Diferencias en los ángulos:
double thetadiff = 0.0;
double phidiff = 0.05;

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
    k[1][1] = h*dphi    (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);
    k[1][2] = h*dmtheta (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);
    k[1][3] = h*dmphi   (a[0]+k[0][0]*0.5, a[1]+k[0][1]*0.5, a[2]+k[0][2]*0.5, a[3]+k[0][3]*0.5);

    //Calculo k3: CORREGIDO - era k[2][3] en el último argumento
    k[2][0] = h*dtheta  (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][1] = h*dphi    (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][2] = h*dmtheta (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);
    k[2][3] = h*dmphi   (a[0]+k[1][0]*0.5, a[1]+k[1][1]*0.5, a[2]+k[1][2]*0.5, a[3]+k[1][3]*0.5);

    //Calculo k4:
    k[3][0] = h*dtheta  (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][1] = h*dphi    (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][2] = h*dmtheta (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);
    k[3][3] = h*dmphi   (a[0]+k[2][0], a[1]+k[2][1], a[2]+k[2][2], a[3]+k[2][3]);

    //Actualizo los resultados de a.
    for(int i=0; i<4; i++)
    {
        a[i]+=1.0/6.0*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i]);
    }
}

// ################### FUNCIONES PARA LYAPUNOV #####################

double calcular_separacion(double y1[4], double y2[4], double vel1[2], double vel2[2])
{
    // Calculo la separación en el espacio de fases (4D)
    return sqrt(pow(y1[0]-y2[0], 2) + pow(y1[1]-y2[1], 2) + 
                pow(vel1[0]-vel2[0], 2) + pow(vel1[1]-vel2[1], 2));
}

void renormalizar_trayectoria(double y_ref[4], double y_pert[4], double epsilon)
{
    // Renormaliza la trayectoria perturbada para mantener separación pequeña
    double separacion_actual = sqrt(pow(y_ref[0]-y_pert[0], 2) + pow(y_ref[1]-y_pert[1], 2) + 
                                   pow(y_ref[2]-y_pert[2], 2) + pow(y_ref[3]-y_pert[3], 2));
    
    if(separacion_actual > 0) {
        double factor = epsilon / separacion_actual;
        for(int i = 0; i < 4; i++) {
            y_pert[i] = y_ref[i] + (y_pert[i] - y_ref[i]) * factor;
        }
    }
}

int main(void)
{
    FILE* lyapunovtxt = fopen("lyapunov.txt", "w");
    
    printf("Calculando exponente de Lyapunov para el péndulo doble...\n");
    printf("Rango de tiempos: %d a %d segundos\n", T_Total, Tmaximo);
    
    for (int tiempo = T_Total; tiempo < Tmaximo; tiempo = tiempo + incremento)
    {
        double suma_logaritmos = 0.0;
        double epsilon = 1e-12;  // Separación inicial muy pequeña
        int num_renormalizaciones = 0;
        
        printf("Calculando para tiempo = %d s\n", tiempo);
        
        //Condiciones iniciales del primer péndulo (referencia)
        y[0] = thetaini;
        y[1] = phiini;
        
        //Condiciones iniciales del segundo péndulo (perturbado)
        l[0] = thetaini + epsilon;  // Perturbación muy pequeña
        l[1] = phiini;
        
        // Inicialización correcta de momentos usando la energía
        double arg1 = E - 2*g*(1-cos(y[0])) - g*(1-cos(y[1]));
        double arg2 = E - 2*g*(1-cos(l[0])) - g*(1-cos(l[1]));
        
        if(arg1 <= 0 || arg2 <= 0) {
            printf("Error: Energía insuficiente para las condiciones iniciales\n");
            continue;
        }
        
        y[2] = 2*sqrt(arg1);
        y[3] = sqrt(arg1)*cos(y[0]-y[1]);
        
        l[2] = 2*sqrt(arg2);
        l[3] = sqrt(arg2)*cos(l[0]-l[1]);
        
        double t = 0.0;
        int paso = 0;
        int pasos_renorm = 10;  // Renormalizar cada 10 pasos
        
        while(t < tiempo)
        {
            // Evolucionar ambas trayectorias
            rungekutta(y);
            rungekutta(l);
            
            // Renormalizar periódicamente
            if(paso % pasos_renorm == 0 && paso > 0) {
                // Calcular velocidades actuales
                vel[0] = dtheta(y[0], y[1], y[2], y[3]);
                vel[1] = dphi(y[0], y[1], y[2], y[3]);
                vell[0] = dtheta(l[0], l[1], l[2], l[3]);
                vell[1] = dphi(l[0], l[1], l[2], l[3]);
                
                // Calcular separación en espacio de fases
                double separacion = calcular_separacion(y, l, vel, vell);
                
                if(separacion > 0) {
                    suma_logaritmos += log(separacion / epsilon);
                    num_renormalizaciones++;
                    
                    // Renormalizar la trayectoria perturbada
                    renormalizar_trayectoria(y, l, epsilon);
                }
            }
            
            t += h;
            paso++;
        }
        
        // Calcular el exponente de Lyapunov
        double Lyapunov = suma_logaritmos / (double)tiempo;
        
        fprintf(lyapunovtxt, "%d %.8f %d\n", tiempo, Lyapunov, num_renormalizaciones);
        printf("Tiempo: %d s, Lyapunov: %.6f s^-1, Renormalizaciones: %d\n", 
               tiempo, Lyapunov, num_renormalizaciones);
        fflush(stdout);
    }
    
    fclose(lyapunovtxt);
    
    printf("\nCálculo completado. Resultados guardados en 'lyapunov.txt'\n");
    printf("Formato: [tiempo] [exponente_lyapunov] [num_renormalizaciones]\n");
    
    return 0;
}