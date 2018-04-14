author: 
wei li      wli10@student.unimelb.edu.au
side lu     sidel@student.unimelb.edu.au

date: 21/10/2017    

This program intends to calculate the heat distribution of a plane.

The program consists of 4 pieces of code, namely
    heat.c      the main program using MPI and OMP to do the computations
    compile.sh  the script used to compile heat.c
    script.sh   the script used to submit the job in slurm system
    plot.py     a python program used to plot the result

To run the program, you need to go through following steps:
    1.  use command 
            $ bash compile.sh
        to compile the heat.c

    2.  use command
            $ sbatch script.sh
        to require resources and submit the job in slurm system. It
        would give you a job ID, and after finishing the job, it would 
        create an output file in the name of < slurm-jobID.out >

    3.  use command
            $ python plot.py < slurm-jobID.out >
        to plot the result.
        It should be noticed that plot.py needs to be run in the environment
        of python3, and it requires the libraries of "numpy" and "matplotlib"
        
