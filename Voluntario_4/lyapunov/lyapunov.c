
/*

Este programa hace el pendulo doble 2 VECES. Y se calcula la separación para ser usada en el coeficiente de Lyapunov.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// DEFINIMOS LAS CONSTANTES

#define g 9.81
#define PI 3.14159265
#define T_TOTAL 10 //tiempo total.
#define h 0.01 //paso temporal

//Los parámetros del pendulo (según el voluntario son =1 para simplifcar el problema)

#define m1 1.0 // masa del primer pendulo
#define m2 1.0 // masa del segundo pendulo
#define l1 1.0 // longitud del primer pendulo
#define l2 1.0 // longitud del segundo pendulo
#define E 1.0 //Energía total del sistema.

//CONDICIONES INCIALES

double thetaini = PI/16;  //Ángulo inicial en theta
double phiini = PI/16;    //Ángulo inicial en phi

double difftheta = 0.001; //Diferencia de ángulo inicial en theta.
double diffphi = 0.001; //Diferencia de ángulo inicial en phi.



//Creo el vector que tendrá las coordenadas: [theta, phi, momento de theta, momento de phi]

double y[4];

//También creo que el segundo vector

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
    /*
    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    double aux=s/(d*d)*(2.0*mtheta*mphi + c*c*mtheta*mphi - c*(mtheta*mtheta + 2.0*mphi*mphi))-2.0*g*sin(theta);
    return aux;
    */

    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    
    // Término que viene de ∂H/∂θ
    double term1 = s * c / (d*d) * (mtheta - mphi*c) * (2*mphi - mtheta*c);
    double term2 = -2.0*g*sin(theta);
    
    return term1 + term2;

}

double dmphi(double theta, double phi, double mtheta, double mphi)
{   /*
    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    double aux = -s/(d*d)*(2.0*mtheta*mphi + c*c*mtheta*mphi - c*(mtheta*mtheta + 2.0*mphi*mphi))-g*sin(phi);
    return aux;
    */
    double c = cos(theta-phi);
    double s = sin(theta-phi);
    double d = 2.0 - c*c;
    
    // Término que viene de ∂H/∂ψ  
    double term1 = -s * c / (d*d) * (mtheta - mphi*c) * (2*mphi - mtheta*c);
    double term2 = -g*sin(phi);
    
    return term1 + term2;
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

    
    /*
    //k1
    k[0][0] = h*dtheta(vector[0], vector[1], vector[2], vector[3]);
    k[0][1] = h*dphi(vector[0], vector[1], vector[2], vector[3]);
    k[0][2] = h*dmtheta(vector[0], vector[1], vector[2], vector[3]);
    k[0][3] = h*dmphi(vector[0], vector[1], vector[2], vector[3]);

    //k2
    k[1][0]= h*dtheta(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][1]= h*dphi(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][2]= h*dmtheta(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);
    k[1][3]= h*dmphi(vector[0]+k[0][0]/2, vector[1]+k[0][1]/2, vector[2]+k[0][2]/2, vector[3]+k[0][3]/2);

    //k3
    k[2][0]= h*dtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][1]= h*dtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][2]= h*dtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);
    k[2][3]= h*dtheta(vector[0]+k[1][0]/2, vector[1]+k[1][1]/2, vector[2]+k[1][2]/2, vector[3]+k[1][3]/2);

    //k4
    k[3][0]= h*dtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][1]= h*dtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][2]= h*dtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);
    k[3][3]= h*dtheta(vector[0]+k[2][0], vector[1]+k[2][1], vector[2]+k[2][2], vector[3]+k[2][3]);

    */


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
    FILE *posicionestxt = fopen("posiciones.txt", "w");

    FILE *diffangulos = fopen("diferenciaangulos.txt", "w");
    FILE *posiciones2txt = fopen("posiciones2.txt", "w");

    if (angulostxt == NULL || momentostxt == NULL || hamiltonianotxt == NULL ||  posicionestxt == NULL || diffangulos == NULL || posiciones2txt == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }



    //El vector pos[4] indica las posiciones en el eje X,Y de la particula 1 y la particula 2. siendo [x1, y1, x2, y2].
    //El vector posly es del otro vector, y el de diffv es la diferencia entre ángulos [theta1 - theta2, phi1 - phi2]

    double pos[4];
    double posl[4];

    double diffv[2];

    //Establecemos las condiciones iniciales. Tenemos 4 parámetros libres pero lo reducimos a 2: los ángulos iniciales.
    // Para ello asumimos phi'=0, por tanto calculamos su momento -> obtenemos mphi.
    //También tenemos que la energía tiene un valor constante, que usamos para sacar el otro valor -> obtenemos mtheta.

    y[0]=thetaini; //Ángulo theta.
    y[1]=phiini; //Ángulo phi.

    l[0]=thetaini - difftheta;
    l[1]=phiini - diffphi;

    //Evaluamos la diferencia inicial.

    diffv[0] = y[0]-l[0];
    diffv[1] = y[1]-l[1];

    fprintf(diffangulos, "%lf %lf\n", diffv[0], diffv[1]);


    double arg=E-2*g*(1-cos(y[0]))-g*(1-cos(y[1]));

    double arg2 = E-2*g*(1-cos(l[0]))-g*(1-cos(l[1]));
    if (arg<0)
    {
        printf("El argumento es negativo, no podemos hacer eso");
        return 1;
    }
    y[2]=2*sqrt(arg); //Momento theta.
    y[3]=y[2]*0.5*cos(y[0]-y[1]); //Momento phi.

    l[2]=2*sqrt(arg);
    l[3]=l[2]*0.5*cos(l[0]-l[1]); 



    //Hago el bucle principal.

    int numpasos = T_TOTAL/h;

    //Condiciones iniciales;

    //y[0]=0.25; //Ángulo theta.
    //y[1]=0.3; //Ángulo phi.


    /* SI QUEREMOS SOLO EJECUTAR EL PROGRAMA 1 VEZ*/
    
    for(int i=0; i<numpasos; i++)
    {
        rungekutta(y);
        rungekutta(l);

        //sumo las energías para
        
        double Hamil = 0.0;
        Hamil = hamiltoniano(y);
        Hamil += hamiltoniano(l);

        fprintf(hamiltonianotxt, "%lf\n", Hamil);
        fprintf(angulostxt, "%lf %lf\n", y[0], y[1]);
        fprintf(momentostxt, "%lf %lf\n", y[2], y[3]);
        fprintf(diffangulos, "%lf %lf\n", diffv[0], diffv[1]);

        //Calculo las posiciones para la animación. Para ello asumo una caja de 6x6, en cuyo centro se encuentra el pendulo.

        pos[0]= 3 + sin(y[0]);
        pos[1]= 3 - cos(y[0]);
        pos[2]= pos[0] + sin(y[1]);
        pos[3]= pos[1] - cos(y[1]);

        posl[0]= 3 + sin(l[0]);
        posl[1]= 3 - cos(l[0]);
        posl[2]= posl[0] + sin(l[1]);
        posl[3]= posl[1] - cos(l[1]);

        fprintf(posicionestxt, "%lf %lf %lf %lf\n", pos[0], pos[1], pos[2], pos[3]);
        fprintf(posiciones2txt, "%lf %lf %lf %lf\n %lf %lf %lf %lf\n \n", pos[0], pos[1], pos[2], pos[3], posl[0], posl[1], posl[2], posl[3]);


    }
    


    fclose(angulostxt);
    fclose(momentostxt);
    fclose(hamiltonianotxt);
    fclose(diffangulos);
    fclose(posiciones2txt);
    fclose(posicionestxt);


}