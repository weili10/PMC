#include <stdio.h>
#include <omp.h>

#define MAX_THREAD 8
#define NUM_STEPS 100000;


int main()
{
    int i, t;
    double step, st, et, x, pi = 0.0;
    double sum[MAX_THREAD] = {0.0};

    omp_set_num_threads(MAX_THREAD);

    step = 1.0/(double) NUM_STEPS;
    st = omp_get_wtime();

    #pragma omp parallel for
    for (int i = 0; i < NUM_STEPS; ++i)
    {
        t = omp_get_thread_num();
        x = (i+0.5) * step;
        sum[t] = sum[t] + 4.0/(1.0 + x*x);
    }

    et = omp_get_wtime();

    for (int i = 0; i < MAX_THREAD; ++i)
    {
        pi = pi + step * sum[i];
        printf("sum[%d] = %lf\n", i, sum[i]);
    }

    printf("pi = %lf\n", pi);
    printf("run time = %lf s\n", et-st);

    return 0;
}