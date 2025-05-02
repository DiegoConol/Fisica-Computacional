#include <stdio.h>
#include <string.h>
#include <math.h> //Necesario para la potencia

typedef struct {
    double rx, ry;      // Posición reescalada: r' = r / C
    double vx, vy;      // Velocidad reescalada
    double ax, ay;      // Aceleración reescalada
    double m;           // Masa reescalada: m' = m / MS
    double e;          // Excentricidad
    double t;          // Período reescalado:
} Body;

#define G 6.67e-11        // Constante de gravitación en N·m²/kg² (no usada directamente aquí)
#define MS 1.99e30        // Masa del Sol en kg
#define C 1.496e11        // Unidad de distancia: distancia Tierra-Sol en metros

#define planetas 10              // Número de cuerpos (por ejemplo, Sol y Tierra)
#define T_TOTAL 100.0    // Tiempo total de simulación (en unidades reescaladas)
#define h 0.1           // Paso temporal

int N = planetas; //Número de planetas, por si acaso.
void aceleracion(Body cuerpos[]) 
{
    double distcubo;
    double dx, dy;
    //Inicializamos la aceleración a 0 SKIPEABLE.
    for (int i=0; i<planetas; i++) 
    {
        cuerpos[i].ax = 0.0;
        cuerpos[i].ay = 0.0; 
    }

    //Calculamos la aceleración de cada planeta.

    for (int i=0; i<planetas; i++)
    {
        for (int j=0; j<planetas; j++)
        {
            if (i != j) //No calculamos la aceleración de un planeta sobre sí mismo.
            {
                dx = cuerpos[i].rx-cuerpos[j].rx; //Diferencia en x
                dy = cuerpos[i].ry-cuerpos[j].ry; //Diferencia en y
                distcubo = pow(dx*dx + dy*dy, 1.5); //Distancia al cubo
                cuerpos[i].ax -= cuerpos[j].m*dx/distcubo; //Aceleración en x
                cuerpos[i].ay -= cuerpos[j].m*dy/distcubo; //Aceleración en y
            }
        }
    }
    return;
    
}


void verlet(Body cuerpos[], FILE *file)
{
    double acelvieja[planetas][2]; //Aceleración vieja de cada planeta.
    double omega[planetas][2]; //Vector auxiliar
    //Actualizamos posiciones y calculamos vector auxiliar omega
    for (int i=0; i<planetas; i++)
    {
        cuerpos[i].rx += cuerpos[i].vx * h + 0.5 * cuerpos[i].ax * h * h;
        cuerpos[i].ry += cuerpos[i].vy * h + 0.5 * cuerpos[i].ay * h * h;
        omega[i][0]=cuerpos[i].vx + h/2*cuerpos[i].ax;
        omega[i][1]=cuerpos[i].vy + h/2*cuerpos[i].ay; 
    }

    //Hemos obtenido r(t+h) y omega(t+h).

    //Actualizamos las aceleraciones (a(t+h)).
    aceleracion(cuerpos);

    //Por último actualizamos las velocidades.
    for (int i=0; i<planetas; i++)
    {
        cuerpos[i].vx = omega[i][0] + h/2*cuerpos[i].ax;
        cuerpos[i].vy = omega[i][1] + h/2*cuerpos[i].ay;
    }
    //Con esto ya tenemos r(t+h), v(t+h) y a(t+h). Dichos parámetros se han actualizado en el cuerpo.
    for (int i=0; i<planetas; i++)
    {
        fprintf(file, "%lf, %lf\n", cuerpos[i].rx, cuerpos[i].ry);
    }
    fprintf(file, "\n"); //Salto de línea para separar los pasos
}
/*
void inicializarCuerpos(Body cuerpos[]) {
    // Cuerpo 0: Sol (centrado en el origen)
    cuerpos[0].rx = 0.0;
    cuerpos[0].ry = 0.0;
    cuerpos[0].vx = 0.0;
    cuerpos[0].vy = 0.0;
    cuerpos[0].m  = 1.0;   // m' = 1 (masa solar)

    // Cuerpo 1: Planeta 
    cuerpos[1].rx = 0.387;   
    cuerpos[1].ry = 0.0;
    cuerpos[1].vx = 0.0;
    cuerpos[1].vy = sqrt(1.0 / fabs(cuerpos[1].rx)); 
    cuerpos[1].m  = 1.658e-7;  
    
     // Cuerpo 2: Venus 
     cuerpos[2].rx = 0.723;   
     cuerpos[2].ry = 0.0;
     cuerpos[2].vx = 0.0;
     cuerpos[2].vy = sqrt(1.0 / fabs(cuerpos[2].rx)); 
     cuerpos[2].m  = 2.447e-6;

      // Cuerpo 3: Tierra 
    cuerpos[3].rx = 1.0;   
    cuerpos[3].ry = 0.0;
    cuerpos[3].vx = 0.0;
    cuerpos[3].vy = sqrt(1.0 / fabs(cuerpos[3].rx)); 
    cuerpos[3].m  = 3.0e-6;

     // Cuerpo 4: Marte 
     cuerpos[4].rx = 1.524;   
     cuerpos[4].ry = 0.0;
     cuerpos[4].vx = 0.0;
     cuerpos[4].vy = sqrt(1.0 / fabs(cuerpos[4].rx)); 
     cuerpos[4].m  = 3.227e-7;

      // Cuerpo 5: Júpiter 
    cuerpos[5].rx = 5.203;   
    cuerpos[5].ry = 0.0;
    cuerpos[5].vx = 0.0;
    cuerpos[5].vy = sqrt(1.0 / fabs(cuerpos[5].rx)); 
    cuerpos[5].m  = 9.547e-4;

     // Cuerpo 6: Saturno 
     cuerpos[6].rx = 9.582;   
     cuerpos[6].ry = 0.0;
     cuerpos[6].vx = 0.0;
     cuerpos[6].vy = sqrt(1.0 / fabs(cuerpos[6].rx)); 
     cuerpos[6].m  = 2.857e-4;

      // Cuerpo 7: Urano 
    cuerpos[7].rx = 19.22;   
    cuerpos[7].ry = 0.0;
    cuerpos[7].vx = 0.0;
    cuerpos[7].vy = sqrt(1.0 / fabs(cuerpos[7].rx)); 
    cuerpos[7].m  = 4.366e-5;

     // Cuerpo 8: Neptuno 
     cuerpos[8].rx = 30.070;   
     cuerpos[8].ry = 0.0;
     cuerpos[8].vx = 0.0;
     cuerpos[8].vy = sqrt(1.0 / fabs(cuerpos[8].rx)); 
     cuerpos[8].m  = 5.151e-5;
}

*/

void inicializarCuerpos(Body cuerpos[], int N, FILE *archivo) 
{
    int i; 
    
    for(i=0; i<N; i++){
        fscanf(archivo, "%lf %lf %lf", 
            &cuerpos[i].m, &cuerpos[i].rx, &cuerpos[i].ry);
        
        cuerpos[i].vx = 0.0;
        cuerpos[i].vy = sqrt(1.0 / fabs(cuerpos[i].rx)) ;
    }
    cuerpos[0].vx = 0.0;
    cuerpos[0].vy = 0.0;

    cuerpos[9].vy=0.31;
    cuerpos[9].vx=-0.275;
    /*
    cuerpos[10].vx = cuerpos[5].vx;
    cuerpos[10].vy = cuerpos[5].vy;
    cuerpos[11].vx = cuerpos[5].vx;
    cuerpos[11].vy = cuerpos[5].vy;
    */
}

//Ahora para calcular la energía

void Energia(Body cuerpos[], int N, FILE *archivo)
{
    double E[N];
    double l;
    for (int i=0; i<N; i++)
    {
        l=cuerpos[i].m*(cuerpos[i].rx*cuerpos[i].vy-cuerpos[i].ry*cuerpos[i].vx);
        E[i]=((cuerpos[i].e*cuerpos[i].e)-1)*(cuerpos[i].m*cuerpos[i].m*cuerpos[i].m)/(2*l*l);
        fprintf(archivo, "%e\n", E[i]);
    }
    fprintf(archivo, "\n"); //Salto de línea para separar los pasos
}


//Esta energía es mucho mejor
void computeTotalEnergy(Body bodies[], int N_BODIES, FILE *archivo) {
    int i, j;
    double kinetic = 0.0, potential = 0.0;
    double energy=0.0;
    // Energía cinética: (1/2) m v²
    for (i = 0; i < N_BODIES; i++) {
        double v2 = bodies[i].vx * bodies[i].vx + bodies[i].vy * bodies[i].vy;
        kinetic += 0.5 * bodies[i].m * v2;
    }

    // Energía potencial: -G * m_i * m_j / r_ij, con G=1 en unidades reescaladas
    for (i = 0; i < N_BODIES; i++) {
        for (j = i + 1; j < N_BODIES; j++) {
            double dx = bodies[i].rx - bodies[j].rx;
            double dy = bodies[i].ry - bodies[j].ry;
            double r = sqrt(dx * dx + dy * dy);
            if (r > 1e-10) { // Evitar división por cero
                potential -= bodies[i].m * bodies[j].m / r;
            }
        }
    }

    energy = kinetic + potential;
    fprintf(archivo, "%e\n", energy);
}














int main (void)
{
    FILE *archivo = fopen("inicialreescalado.txt", "r");
    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    FILE *file = fopen("SALIDA.txt", "w");
    if (file == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    FILE *fileenergia = fopen("SuperEnergia.txt", "w");
    if (fileenergia == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }


    Body cuerpos[N]; // Arreglo de cuerpos
    inicializarCuerpos(cuerpos, N, archivo);
    aceleracion(cuerpos);
    for(float i=0; i<T_TOTAL; i=i+h)
    {   
        verlet(cuerpos, file);
        computeTotalEnergy(cuerpos, N, fileenergia);
    }
    fclose(file);
    fclose(fileenergia);
    fclose(archivo);

    return 0;
}
