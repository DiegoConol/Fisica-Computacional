/*

Este es un programa que hace un pendulo doble y exporta los datos de los angulos, angulo's y hamiltoniano a un archivo.
Versión mejorada con cálculo correcto de diferencias para exponente de Lyapunov.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// DEFINIMOS LAS CONSTANTES

#define g 9.81
#define PI 3.14159265
#define T_TOTAL 600 //tiempo total. (No recomendable poner más de 30, sino se raya. Poner 30 o bajar el paso.)
#define h 0.01 //paso temporal

//Los parámetros del pendulo (según el voluntario son =1 para simplifcar el problema)

#define m1 1.0 // masa del primer pendulo
#define m2 1.0 // masa del segundo pendulo
#define l1 1.0 // longitud del primer pendulo
#define l2 1.0 // longitud del segundo pendulo
double E = 0.0;


// ##########   DEFINE LAS CONDICIONES INICIALES PARA CADA ENERGÍA: ###############

double thetaini = PI/16;
double phiini = PI/16;

//Diferencias en las condiciones iniciales para el segundo péndulo.
double difftheta = 0.0;
double diffphi = 0.0001; 

//Creo el vector que tendrá las coordenadas: [theta, phi, momento de theta, momento de phi]

double y[4];

//El vector para el segundo pendulo:
double l[4];


//Creo la función que calcula el hamiltoniano.

double hamiltoniano(double y[4])
{

    double H=0.0;
    double c = cos(y[0]-y[1]);
    double dtheta1 = 1/(2-c*c)*(y[2]-y[3]*c);
    double dphi1 = 1/(2-c*c)*(2*y[3]-y[2]*c);

    H= dtheta1*dtheta1 +0.5*dphi1*dphi1 + dtheta1*dphi1*c+2*g*(1-cos(y[0]))+g*(1-cos(y[1]));
    return H;

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

// Función para normalizar ángulos entre -π y π
double normalize_angle(double angle) {
    while (angle > PI) angle -= 2*PI;
    while (angle < -PI) angle += 2*PI;
    return angle;
}

// Función para calcular la diferencia angular considerando periodicidad
double angle_difference(double a1, double a2) {
    double diff = a1 - a2;
    return normalize_angle(diff);
}


int main(void)
{
    FILE *pendulotxt = fopen("pendulo.txt", "w");
    FILE *angulosenergiastxt = fopen("angulosenergias.txt", "w");
    FILE *energiasdistintastxt = fopen("energiasdistintas.txt", "w");

    //El vector pos[4] indica las posiciones en el eje X,Y de la particula 1 y la particula 2. siendo [x1, y1, x2, y2].
    double pos[4];
    double posl[4];
    double diffv[4];

    //El vector vel[2] indica las velocidades de los ángulos [theta', phi']
    double vel[2];

    //El valor auxiliar de diferencia total para los angulos:
    double diferenciatotal = 0.0;

    int numpasos = T_TOTAL/h;

    // Condiciones iniciales
    double energias[] = {1.0, 3.0, 5.0, 10.0, 15.0};
    int num_energias = sizeof(energias)/sizeof(energias[0]);

    for (int e = 0; e < num_energias; ++e) {
        double E = energias[e];
        char fname_angulos[64], fname_lyapunov[64], fname_momentos[64], fname_hamiltoniano[64], fname_posiciones[64], fname_diffangulos[64], fname_posiciones2[64];

        snprintf(fname_angulos, sizeof(fname_angulos), "angulos_%.1f.txt", E);
        snprintf(fname_momentos, sizeof(fname_momentos), "momentos_%.1f.txt", E);
        snprintf(fname_hamiltoniano, sizeof(fname_hamiltoniano), "hamiltoniano_%.1f.txt", E);
        snprintf(fname_posiciones, sizeof(fname_posiciones), "posiciones_%.1f.txt", E);
        snprintf(fname_diffangulos, sizeof(fname_diffangulos), "diferenciaangulos_%.1f.txt", E);
        snprintf(fname_posiciones2, sizeof(fname_posiciones2), "posiciones2_%.1f.txt", E);
        snprintf(fname_lyapunov, sizeof(fname_lyapunov), "lyapunov_%.1f.txt", E);



        FILE *angulostxt = fopen(fname_angulos, "w");
        FILE *momentostxt = fopen(fname_momentos, "w");
        FILE *hamiltonianotxt = fopen(fname_hamiltoniano, "w");
        FILE *posicionestxt = fopen(fname_posiciones, "w");
        FILE *diffangulos = fopen(fname_diffangulos, "w");
        FILE *posiciones2txt = fopen(fname_posiciones2, "w");
        FILE *lyapunovtxt = fopen(fname_lyapunov, "w");

        // Condiciones iniciales
        y[0]=thetaini; // Ángulo theta
        y[1]=phiini;  // Ángulo phi

        l[0]=thetaini + difftheta; // Corregido: suma en lugar de resta
        l[1]=phiini + diffphi;     // Corregido: suma en lugar de resta

        //Asignamos el coeficiente de lyapunov para cada energía y lo guardamos en cada fichero.
        double lyapunov_sum = 0.0;

        double arg = E - 2*g*(1-cos(y[0])) - g*(1-cos(y[1]));
        double arg2 = E - 2*g*(1-cos(l[0])) - g*(1-cos(l[1]));
        
        if (arg < 0 || arg2 < 0) {
            printf("El argumento es negativo para E=%.1f, no podemos hacer eso\n", E);
            fclose(angulostxt);
            fclose(momentostxt);
            fclose(hamiltonianotxt);
            fclose(posicionestxt);
            fclose(diffangulos);
            fclose(posiciones2txt);
            fclose(lyapunovtxt);
            continue;
        }
        
        y[2]=2*sqrt(arg); // Momento theta
        y[3]=y[2]*0.5*cos(y[0]-y[1]); // Momento phi

        l[2]=2*sqrt(arg2);
        l[3]=l[2]*0.5*cos(l[0]-l[1]);

        // Cálculo inicial de diferencias usando diferencias angulares correctas
        diffv[0] = fabs(angle_difference(y[0], l[0]));
        diffv[1] = fabs(angle_difference(y[1], l[1]));
        
        double ypunto[2]; // velocidad en [thetapunto, phipunto];
        double lpunto[2];

        ypunto[0] = 1/(2-cos(y[0]-y[1])*cos(y[0]-y[1]))*(y[2]-y[3]*cos(y[0]-y[1]));
        ypunto[1] = 1/(2-cos(y[0]-y[1])*cos(y[0]-y[1]))*(y[3]*2-y[2]*cos(y[0]-y[1]));

        lpunto[0] = 1/(2-cos(l[0]-l[1])*cos(l[0]-l[1]))*(l[2]-l[3]*cos(l[0]-l[1]));
        lpunto[1] = 1/(2-cos(l[0]-l[1])*cos(l[0]-l[1]))*(l[3]*2-l[2]*cos(l[0]-l[1]));

        diffv[2] = fabs(ypunto[0]-lpunto[0]);
        diffv[3] = fabs(ypunto[1]-lpunto[1]);

        //Para lyapunov necesitamos delta(t)=diferenciatotal y delta(0) = diferenciainicial. Luego sacaremos el lambda.
        double diferenciainicial = sqrt(diffv[0]*diffv[0]+diffv[1]*diffv[1]+diffv[2]*diffv[2]+diffv[3]*diffv[3]);
        diferenciatotal = diferenciainicial;

        fprintf(diffangulos, "%lf\n", diferenciatotal);

        for(int i=0; i<numpasos; i++) {

            //Hacemos el algoritmo de rungekutta y la energía.
            rungekutta(y);
            rungekutta(l);

            double Hamil = 0.0;
            Hamil = hamiltoniano(y);
            
            fprintf(angulostxt, "%lf %lf\n", y[0], y[1]);
            fprintf(momentostxt, "%lf %lf\n", y[2], y[3]);
            fprintf(hamiltonianotxt, "%lf\n", Hamil);

            //Calculo la diferencia de ángulos usando diferencias angulares correctas
            diffv[0] = fabs(angle_difference(y[0], l[0]));
            diffv[1] = fabs(angle_difference(y[1], l[1]));
            
            ypunto[0] = 1/(2-cos(y[0]-y[1])*cos(y[0]-y[1]))*(y[2]-y[3]*cos(y[0]-y[1]));
            ypunto[1] = 1/(2-cos(y[0]-y[1])*cos(y[0]-y[1]))*(y[3]*2-y[2]*cos(y[0]-y[1]));

            lpunto[0] = 1/(2-cos(l[0]-l[1])*cos(l[0]-l[1]))*(l[2]-l[3]*cos(l[0]-l[1]));
            lpunto[1] = 1/(2-cos(l[0]-l[1])*cos(l[0]-l[1]))*(l[3]*2-l[2]*cos(l[0]-l[1]));

            diffv[2] = fabs(ypunto[0]-lpunto[0]);
            diffv[3] = fabs(ypunto[1]-lpunto[1]);
            
            diferenciatotal = sqrt(diffv[0]*diffv[0]+diffv[1]*diffv[1]+diffv[2]*diffv[2]+diffv[3]*diffv[3]);

            // Cálculo correcto del exponente de Lyapunov
            if (diferenciatotal > 0 && diferenciainicial > 0) {
                lyapunov_sum += log(diferenciatotal/diferenciainicial);
            }
            
            double lyapunov_actual = lyapunov_sum / ((i+1) * h);
            fprintf(lyapunovtxt, "%lf\n", lyapunov_actual);
            fprintf(diffangulos, "%lf\n", diferenciatotal);

            // Calculo las posiciones para la animación de ambos pendulos
            pos[0]= 3 + sin(y[0]);
            pos[1]= 3 - cos(y[0]);
            pos[2]= pos[0] + sin(y[1]);
            pos[3]= pos[1] - cos(y[1]);

            posl[0]= 3 + sin(l[0]);
            posl[1]= 3 - cos(l[0]);
            posl[2]= posl[0] + sin(l[1]);
            posl[3]= posl[1] - cos(l[1]);

            fprintf(posicionestxt, "%lf %lf %lf %lf\n", pos[0], pos[1], pos[2], pos[3]);
            fprintf(posiciones2txt, "%lf %lf %lf %lf\n%lf %lf %lf %lf\n\n", pos[0], pos[1], pos[2], pos[3], posl[0], posl[1], posl[2], posl[3]);
        }

        fclose(angulostxt);
        fclose(momentostxt);
        fclose(hamiltonianotxt);
        fclose(posicionestxt);
        fclose(diffangulos);
        fclose(posiciones2txt);
        fclose(lyapunovtxt);

        printf("Una energía más hecha :D (E=%.1f)\n", E);
        fflush(stdout);
    }

    fclose(pendulotxt);
    fclose(angulosenergiastxt);
    fclose(energiasdistintastxt);

    return 0;
}