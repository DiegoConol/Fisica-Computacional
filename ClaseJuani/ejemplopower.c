#include <stdio.h>

int power(int base, int n);
int power(int base, int n)
{

    int i, p;
    p = 1;
    for (i=1; i<=n; ++i) {
        p=p*base;
    }
    return p;
}
int main()
//prueba la función power
{
    int i;
    for (i=0; i<10; ++i){
        printf("%i\t%i\t%i\n", i, power(2,i), power(-3,i));
    }
    return 0;

}

