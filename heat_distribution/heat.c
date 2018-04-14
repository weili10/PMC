// author: 
// wei li      wli10@student.unimelb.edu.au
// side lu     sidel@student.unimelb.edu.au

// date: 21/10/2017   

// This program is used to compute the heat distribution.
// It used the Jacobi methods to find the steady state, and used MPI and OpenMP
// to parallelize the computation.

#include<math.h>
#include<stdio.h>
#include"mpi.h"
#include<omp.h>
#define N 1000
#define EPSILON 0.0001
#define THREADS 16

double global_u[N][N];
double local_u[N][N];

void main(int argc, char *argv[])
{
    int i,j;
    double mean,global_mean;  //mean temperature of the boundaries

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);


    int size,rank;

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    omp_set_num_threads(THREADS);

    /** initialize the matrix in this node **/
    int my_rows = (int)N/size;
    if(rank == 0 || rank == size-1){
        my_rows += 1;
    }else{
        my_rows += 2;
    }
    
    double u[my_rows][N];
    double w[my_rows][N];
    
    // initialize the boundaries
    mean = 0.0;
    if(rank == size-1){
        for(j = 0; j<N;j++){
            u[my_rows-1][j] = 0.0;
        }
    }else if(rank == 0){
        for(j = 0; j<N;j++){
            u[0][j] = 100.0;
            mean += 100.0;
        }
    }
    
    for(i = 1; i<my_rows-1;i++){
        u[i][0]  = 100.0;
        u[i][N-1] = 100.0;
        mean+=u[i][0]+u[i][N-1];
    }
    
    MPI_Allreduce (&mean, &global_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
    // initialize the matrix except boundaries
    global_mean/=(4.0*N);
    for(i=1;i<my_rows-1;i++)
        for(j=1;j<N-1;j++)
            u[i][j]=global_mean;


    /** calculate heat distribution, use MPI and openMP **/
    double diff;            //maximum difference in this node
    double tdiff;           //maximum difference in one thread
    double global_diff;     //global maximum difference of all the nodes
    
    MPI_Status status;

    // double startTime, endTime;
    // startTime = MPI_Wtime();

    for (;;) {
        //  MPI exchange ghost row's data
        if (rank > 0){
            MPI_Send (u[1], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if (rank < size-1) {
            MPI_Send (u[my_rows-2], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            MPI_Recv (u[my_rows-1], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD,&status);
        }
        if (rank > 0){
            MPI_Recv (u[0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
        }
        diff = 0.0;

        //OpenMP update heat distribution in parallel
        #pragma omp parallel private (i, j, tdiff)
        {
            tdiff = 0.0;
            #pragma omp for
            for (i = 1; i < my_rows-1; i++){
                for (j = 1; j < N-1; j++) {
                    w[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1])/4.0;
                    if (fabs(w[i][j] - u[i][j]) > tdiff){
                        tdiff = fabs(w[i][j] - u[i][j]);
                    }
                    
                }
            }
            #pragma omp for nowait
            for (i = 1; i < my_rows-1; i++){
                for (j = 1; j < N-1; j++){
                    u[i][j] = w[i][j];
                }
            }

            #pragma omp critical
            if (tdiff > diff) diff = tdiff;
        } 
        MPI_Allreduce (&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (global_diff <= EPSILON) break;
    }
    
    // endTime = MPI_Wtime();

    // printf("node %d runs %f s\n",rank,endTime - startTime);
    
    /** output the result **/

    // MPI gather the data in node 0
    int rowSize = (int)N/size;
    if(rank == 0){
        for(i=0;i<rowSize;i++){
            for(j = 0; j< N; j++){
                local_u[i][j] = u[i][j];
            }
        }
    }else{
        for(i=0;i<rowSize;i++){
            for(j = 0; j< N; j++){
                local_u[i][j] = u[i+1][j];
            }
        }
    }

    MPI_Gather(local_u,rowSize*N,MPI_DOUBLE,global_u,rowSize*N,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // node 0 print the result
    if(rank == 0){
        for(i = 0; i <N; i++){
            for(j = 0; j< N; j++){
                printf("%8.5f\t",global_u[i][j]);
            }
            putchar('\n');
        }
    }

    MPI_Finalize();
    
}