#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double dot_product(double v[], double u[], int n){
    int i;
    double result = 0;
    for (i = 0; i < n; i++)
    {
    result += v[i]*u[i];
    }
    return result;
}

int main(int argc, char*argv[]){
    //FILE *f;
    unsigned long cpu_time_used;
    
    double *ar1 = malloc(10000000 * sizeof(double));
    double *ar2 = malloc(10000000 * sizeof(double));
    if (ar1 == NULL) {
        fprintf(stderr, "malloc failed for ar1\n");
        return 1;
    }
    if (ar2 == NULL) {
        fprintf(stderr, "malloc failed for ar2\n");
        return 1;
    }
    /*f = fopen("D:\\User\\Downloads\\test\\test.txt", "w");
    if (f == NULL)
    {
        fprintf(stderr, "Error opening file!\n");
        exit(1);
    }*/

    srand(time(NULL));
    for (int i = 0; i < 10000000; i++){
        ar1[i] = (double)rand()/RAND_MAX;
        ar2[i] = (double)rand()/RAND_MAX;
    }

    //start = clock();
    struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long startTime = 1000000 * tv.tv_sec + tv.tv_usec;
    double s = dot_product(ar1, ar2, 10000000);
    //end = clock();
    gettimeofday(&tv,NULL);
    unsigned long endTime = 1000000 * tv.tv_sec + tv.tv_usec;
    cpu_time_used = endTime - startTime;
    printf("scalar: %0.3f; elapsed [µs]: %d\n", s, cpu_time_used);
    //fclose(f); 
    free(ar1);
    free(ar2);
    exit(0);
}