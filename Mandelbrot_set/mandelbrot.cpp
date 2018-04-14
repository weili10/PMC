#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>

// return 1 if in set, 0 otherwise
int inset(double real, double img, int maxiter){
    double z_real = real;
    double z_img = img;
    for(int iters = 0; iters < maxiter; iters++){
        double z2_real = z_real*z_real-z_img*z_img;
        double z2_img = 2.0*z_real*z_img;
        z_real = z2_real + real;
        z_img = z2_img + img;
        if(z_real*z_real + z_img*z_img > 4.0) return 0;
    }
    return 1;
}

// using escape time algorithm to count the Mandelbrot Set in the region
// 
// [real_lower+real_num_lower*real_step,
//  real_lower+real_num_upper*real_step,
//  img_lower,
//  img_lower+img_num*img_step]
// 
int mandelbrotSetCount(double real_lower, int real_num_lower, int real_num_upper, double real_step, 
    double img_lower, int img_num, double img_step, int maxiter){
    
    int count=0;
    
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(32) reduction(+: count)
    for(int real=real_num_lower; real<=real_num_upper; real++){  
        for(int img=0; img<=img_num; img++){
            double x = real_lower+real*real_step;
            double y = img_lower+img*img_step;
            double q = (x-0.25)*(x-0.25)+y*y;
            if(q*(q+x-0.25) < 0.25*y*y){
                // (x,y) is in Cardioid
                count++;
            }
            else if(((x+1)*(x+1)+y*y) < 0.0625){
                // (x,y) is in the period-2 bulb
                count++;
            }
            else{
                count+=inset(x,y,maxiter);
            }
        }
    }
    return count;
}


// count the duplicated Mandelbrot Set on line y = img in the x range of 
// 
//  [real_lower+real_num_lower*real_step,
//   real_lower+real_num_upper*real_step]
// 
int duplicateCount(double real_lower,int real_num_low, int real_num_up, double real_step, 
                    double img, int maxiter){

    int count = 0;

    #pragma omp parallel for schedule(dynamic) num_threads(32) reduction(+: count)
    for(int real = real_num_low; real<=real_num_up;real++){
        // count += inset(real_lower+real*real_step, img, maxiter);
        double x = real_lower+real*real_step;
        double y = img;
        // double y = img_lower+img*img_step;
        double q = (x-0.25)*(x-0.25)+y*y;
        if(q*(q+x-0.25) < 0.25*y*y){
            // (x,y) is in Cardioid
            count++;
        }
        else if(((x+1)*(x+1)+y*y) < 0.0625){
            // (x,y) is in the period-2 bulb
            count++;
        }
        else{
            count+=inset(x,y,maxiter);
        }
    }
    return count;
}

// main
int main(int argc, char *argv[]){

    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double real_lower;
    double real_upper;
    double img_lower;
    double img_upper;
    int num;
    int maxiter;
    int num_regions = (argc-1)/6;
    int chunk_num = 200;   // split data set into 100 chunks along the real axis
    
    for(int region=0;region<num_regions;region++){
        // double st,et,tt;
        // st = omp_get_wtime();
        
        // scan the arguments
        sscanf(argv[region*6+1],"%lf",&real_lower);
        sscanf(argv[region*6+2],"%lf",&real_upper);
        sscanf(argv[region*6+3],"%lf",&img_lower);
        sscanf(argv[region*6+4],"%lf",&img_upper);
        sscanf(argv[region*6+5],"%i",&num);
        sscanf(argv[region*6+6],"%i",&maxiter);

        if(num < 200){
            chunk_num = num;
        }

        int localCount = 0;
        int globalCount = 0;
        double real_step = (real_upper-real_lower)/num;
        double img_step = (img_upper-img_lower)/num;

        int numgap = num/chunk_num;
        int chunks_on_node = chunk_num/world_size;
        int real_num_low[chunks_on_node];
        int real_num_up[chunks_on_node];
        // split data set into chunk_num chunks along the real axis
        #pragma omp parallel for num_threads(32)
        for(int i = 0; i < chunks_on_node; i++){
            real_num_low[i]=(world_rank+i*4)*numgap;
            real_num_up[i] = (world_rank+i*4+1)*numgap-1;
        }
        // amend the last boundary, in case that num is not divisible by chunk_num
        if(world_rank == world_size){
            real_num_up[chunks_on_node-1] = num;
        }

        if(img_upper*img_lower<0){
            // apply symmetry property
            if(img_upper>(-img_lower)){
                int img_num = num*(-img_lower)/(img_upper-img_lower);
                for(int i = 0; i < chunks_on_node; i++){
                    localCount += 2*mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,img_lower,img_num,img_step,maxiter)
                                    +mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,-img_lower,num-2*img_num,img_step,maxiter)
                                    -duplicateCount(real_lower,real_num_low[i],real_num_up[i],real_step,0,maxiter)
                                    -duplicateCount(real_lower,real_num_low[i],real_num_up[i],real_step,img_lower,maxiter);
                }   
            }
            else if(img_upper<(-img_lower)){
                int img_num = num*(img_upper)/(img_upper-img_lower);
                for(int i = 0; i < chunks_on_node; i++){
                    localCount += 2*mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,-img_upper,img_num,img_step,maxiter)
                                    +mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,img_lower,num-2*img_num,img_step,maxiter)
                                    -duplicateCount(real_lower,real_num_low[i],real_num_up[i],real_step,0,maxiter)
                                    -duplicateCount(real_lower,real_num_low[i],real_num_up[i],real_step,-img_lower,maxiter);
                }
            }
            else{
                for(int i = 0; i < chunks_on_node; i++){
                    localCount += 2*mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,img_lower,num/2,img_step,maxiter)
                                    -duplicateCount(real_lower,real_num_low[i],real_num_up[i],real_step,0,maxiter);
                }  
            }
        }
        else{
            for(int i = 0; i < chunks_on_node; i++){
                localCount += mandelbrotSetCount(real_lower,real_num_low[i],real_num_up[i],real_step,img_lower,num,img_step,maxiter);
            }
        } 

        // et = omp_get_wtime();
        // printf("proccessor %d runs %.6f s\n",world_rank,et-st);

        //reduce the local count to global count
        // here 1 represents the count, 0 represents the root
        MPI_Reduce(&localCount,&globalCount,1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        // print result in root
        if(world_rank == 0){
            printf("%d\n",globalCount);
            // tt = omp_get_wtime();
            // printf("total runtime %.6f s\n",tt-st);
        }
    }
    MPI_Finalize();

    return EXIT_SUCCESS;
}
