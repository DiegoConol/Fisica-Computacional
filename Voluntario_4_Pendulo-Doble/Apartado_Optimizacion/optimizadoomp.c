/*
Este es un programa que hace un pendulo doble REALMENTE optimizado
Solo paraleliza donde tiene sentido y mejora el rendimiento de verdad.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// DEFINIMOS LAS CONSTANTES
#define g 9.81
#define PI 3.14159265
#define T_TOTAL 10000.0  // tiempo total como double
#define h 0.01         // paso temporal

// Los parámetros del pendulo
#define m1 1.0  // masa del primer pendulo
#define m2 1.0  // masa del segundo pendulo
#define l1 1.0  // longitud del primer pendulo
#define l2 1.0  // longitud del segundo pendulo
#define E 20.0  // Energía total del sistema

// CONDICIONES INICIALES
double thetaini = 1.0;  // Ángulo inicial en theta
double phiini = 1.2;    // Ángulo inicial en phi

// Función inline para el hamiltoniano (más rápida)
static inline double hamiltoniano(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta - phi);
    double d = 2.0 - c*c;
    double dtheta1 = (mtheta - mphi*c) / d;
    double dphi1 = (2.0*mphi - mtheta*c) / d;

    return dtheta1*dtheta1 + 0.5*dphi1*dphi1 + dtheta1*dphi1*c + 
           2.0*g*(1.0 - cos(theta)) + g*(1.0 - cos(phi));
}

// ####################### ECUACIONES DE MOVIMIENTO ######################
// Funciones inline para mejor rendimiento

static inline double dtheta(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta - phi);
    return (mtheta - mphi*c) / (2.0 - c*c);
}

static inline double dphi(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta - phi);
    return (2.0*mphi - mtheta*c) / (2.0 - c*c);
}

static inline double dmtheta(double theta, double phi, double mtheta, double mphi)
{
    double c = cos(theta - phi);
    double s = sin(theta - phi);
    double d = 2.0 - c*c;
    
    double term1 = s * c * (mtheta - mphi*c) * (2.0*mphi - mtheta*c) / (d*d);
    double term2 = -2.0*g*sin(theta);
    
    return term1 + term2;
}

static inline double dmphi(double theta, double phi, double mtheta, double mphi)
{   
    double c = cos(theta - phi);
    double s = sin(theta - phi);
    double d = 2.0 - c*c;
    
    double term1 = -s * c * (mtheta - mphi*c) * (2.0*mphi - mtheta*c) / (d*d);
    double term2 = -g*sin(phi);
    
    return term1 + term2;
}

/* ###################### RUNGE KUTTA SECUENCIAL ################### */
// El Runge-Kutta DEBE ser secuencial, no se puede paralelizar

void rungekutta(double vector[4])
{
    double k1[4], k2[4], k3[4], k4[4];
    double temp[4];
    
    // k1 - estado actual
    k1[0] = h * dtheta(vector[0], vector[1], vector[2], vector[3]);
    k1[1] = h * dphi(vector[0], vector[1], vector[2], vector[3]);
    k1[2] = h * dmtheta(vector[0], vector[1], vector[2], vector[3]);
    k1[3] = h * dmphi(vector[0], vector[1], vector[2], vector[3]);
    
    // k2 - punto medio con k1
    temp[0] = vector[0] + k1[0] * 0.5;
    temp[1] = vector[1] + k1[1] * 0.5;
    temp[2] = vector[2] + k1[2] * 0.5;
    temp[3] = vector[3] + k1[3] * 0.5;
    
    k2[0] = h * dtheta(temp[0], temp[1], temp[2], temp[3]);
    k2[1] = h * dphi(temp[0], temp[1], temp[2], temp[3]);
    k2[2] = h * dmtheta(temp[0], temp[1], temp[2], temp[3]);
    k2[3] = h * dmphi(temp[0], temp[1], temp[2], temp[3]);
    
    // k3 - punto medio con k2
    temp[0] = vector[0] + k2[0] * 0.5;
    temp[1] = vector[1] + k2[1] * 0.5;
    temp[2] = vector[2] + k2[2] * 0.5;
    temp[3] = vector[3] + k2[3] * 0.5;
    
    k3[0] = h * dtheta(temp[0], temp[1], temp[2], temp[3]);
    k3[1] = h * dphi(temp[0], temp[1], temp[2], temp[3]);
    k3[2] = h * dmtheta(temp[0], temp[1], temp[2], temp[3]);
    k3[3] = h * dmphi(temp[0], temp[1], temp[2], temp[3]);
    
    // k4 - final con k3 completo
    temp[0] = vector[0] + k3[0];
    temp[1] = vector[1] + k3[1];
    temp[2] = vector[2] + k3[2];
    temp[3] = vector[3] + k3[3];
    
    k4[0] = h * dtheta(temp[0], temp[1], temp[2], temp[3]);
    k4[1] = h * dphi(temp[0], temp[1], temp[2], temp[3]);
    k4[2] = h * dmtheta(temp[0], temp[1], temp[2], temp[3]);
    k4[3] = h * dmphi(temp[0], temp[1], temp[2], temp[3]);
    
    // Actualización final usando la fórmula de RK4
    const double sixth = 1.0/6.0;
    vector[0] += sixth * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
    vector[1] += sixth * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
    vector[2] += sixth * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]);
    vector[3] += sixth * (k1[3] + 2.0*k2[3] + 2.0*k3[3] + k4[3]);
}

int main(void)
{
    printf("Usando OpenMP con %d threads máximo\n", omp_get_max_threads());
    
    FILE *angulostxt = fopen("angulos.txt", "w"); 
    FILE *momentostxt = fopen("momentos.txt", "w");
    FILE *hamiltonianotxt = fopen("hamiltoniano.txt", "w");
    FILE *posicionestxt = fopen("posiciones.txt", "w");

    if (!angulostxt || !momentostxt || !hamiltonianotxt || !posicionestxt) {
        printf("Error al abrir archivos.\n");
        if (angulostxt) fclose(angulostxt);
        if (momentostxt) fclose(momentostxt);
        if (hamiltonianotxt) fclose(hamiltonianotxt);
        if (posicionestxt) fclose(posicionestxt);
        return 1;
    }

    double y[4];
    
    // Condiciones iniciales
    y[0] = thetaini;
    y[1] = phiini;
    double arg = E - 2.0*g*(1.0 - cos(y[0])) - g*(1.0 - cos(y[1]));
    if (arg < 0.0) {
        printf("El argumento es negativo, no podemos hacer eso\n");
        fclose(angulostxt);
        fclose(momentostxt);
        fclose(hamiltonianotxt);
        fclose(posicionestxt);
        return 1;
    }
    y[2] = 2.0 * sqrt(arg);
    y[3] = y[2] * 0.5 * cos(y[0] - y[1]);

    int numpasos = (int)(T_TOTAL / h);
    
    // Pre-reservar memoria para I/O en lotes (única optimización real)
    const int BATCH_SIZE = 1000;
    double *batch_angulos = malloc(BATCH_SIZE * 2 * sizeof(double));
    double *batch_momentos = malloc(BATCH_SIZE * 2 * sizeof(double));
    double *batch_hamiltonianos = malloc(BATCH_SIZE * sizeof(double));
    double *batch_posiciones = malloc(BATCH_SIZE * 4 * sizeof(double));
    
    if (!batch_angulos || !batch_momentos || !batch_hamiltonianos || !batch_posiciones) {
        printf("Error de memoria.\n");
        if (batch_angulos) free(batch_angulos);
        if (batch_momentos) free(batch_momentos);
        if (batch_hamiltonianos) free(batch_hamiltonianos);
        if (batch_posiciones) free(batch_posiciones);
        fclose(angulostxt);
        fclose(momentostxt);
        fclose(hamiltonianotxt);
        fclose(posicionestxt);
        return 1;
    }

    double start_time = omp_get_wtime();
    
    // Bucle principal - DEBE ser secuencial
    for(int i = 0; i < numpasos; i++) {
        rungekutta(y);
        
        int batch_idx = i % BATCH_SIZE;
        
        // Guardar en batch
        batch_angulos[batch_idx * 2] = y[0];
        batch_angulos[batch_idx * 2 + 1] = y[1];
        batch_momentos[batch_idx * 2] = y[2];
        batch_momentos[batch_idx * 2 + 1] = y[3];
        batch_hamiltonianos[batch_idx] = hamiltoniano(y[0], y[1], y[2], y[3]);
        
        // Calcular posiciones
        batch_posiciones[batch_idx * 4] = 3.0 + sin(y[0]);
        batch_posiciones[batch_idx * 4 + 1] = 3.0 - cos(y[0]);
        batch_posiciones[batch_idx * 4 + 2] = batch_posiciones[batch_idx * 4] + sin(y[1]);
        batch_posiciones[batch_idx * 4 + 3] = batch_posiciones[batch_idx * 4 + 1] - cos(y[1]);
        
        // Escribir batch cuando esté lleno o al final
        if (batch_idx == BATCH_SIZE - 1 || i == numpasos - 1) {
            int write_count = (i == numpasos - 1) ? (batch_idx + 1) : BATCH_SIZE;
            
            // AQUÍ sí podemos paralelizar la escritura
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    for(int j = 0; j < write_count; j++) {
                        fprintf(angulostxt, "%.6f %.6f\n", 
                               batch_angulos[j * 2], batch_angulos[j * 2 + 1]);
                    }
                }
                #pragma omp section
                {
                    for(int j = 0; j < write_count; j++) {
                        fprintf(momentostxt, "%.6f %.6f\n", 
                               batch_momentos[j * 2], batch_momentos[j * 2 + 1]);
                    }
                }
                #pragma omp section
                {
                    for(int j = 0; j < write_count; j++) {
                        fprintf(hamiltonianotxt, "%.6f\n", batch_hamiltonianos[j]);
                    }
                }
                #pragma omp section
                {
                    for(int j = 0; j < write_count; j++) {
                        fprintf(posicionestxt, "%.6f %.6f %.6f %.6f\n", 
                               batch_posiciones[j * 4], batch_posiciones[j * 4 + 1],
                               batch_posiciones[j * 4 + 2], batch_posiciones[j * 4 + 3]);
                    }
                }
            }
        }
    }
    
    double end_time = omp_get_wtime();
    printf("Tiempo de ejecución: %f segundos\n", end_time - start_time);

    // Limpiar memoria
    free(batch_angulos);
    free(batch_momentos);
    free(batch_hamiltonianos);
    free(batch_posiciones);

    fclose(angulostxt);
    fclose(momentostxt);
    fclose(hamiltonianotxt);
    fclose(posicionestxt);
    printf("Terminé el programa.\n");

    return 0;
}