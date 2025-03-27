#include <omp.h>
#include <stdio.h>

int main() {
    #pragma omp parallel
    {
        printf("OpenMP is working!\n");
    }
    return 0;
}