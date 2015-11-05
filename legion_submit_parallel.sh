#!/bin/bash -l

# Batch script to run an MPI parallel job on Legion with the upgraded software
# stack under SGE with OpenMPI.

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=03:00:0

# 3. Request 1 gigabyte of RAM 
#$ -l mem=1G

# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=10G

# 5. Set the name of the job.
#$ -N Test

# 6. Select the QLogic parallel environment (qlc) and 16 process.
#$ -pe qlc 64

# 7. Select the project that this job will run under.
# Find <your_project_id> by running the command "groups"
#$ -P OPTIMETPROJECT

# 8. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/uceealj/Scratch/output/0_scaling/2_weak_scaling/2_Np2/64_Np64

# Set up teh desired module environment 
# compilers
module unload compilers
module load compilers/intel/13.0/028
module load compilers/intel/13.0/028_cxx11
# MPI
module unload mpi
module load mpi/intel/4.1.0/024/testing
# MKL
module unload mkl
module load mkl/11.0/001
# GSL
module load gsl/1.15/intel
# HDF5
module load hdf/5-1.8.7/intel

# 9. Run the application.
gerun $HOME/1_AJ/3_OPTIMET/1_OPTIMET_Parallel_test/5_serial_to_parallel/2_INTEL/0_test_jobs/2_parallel/2_weak_scaling/2_Np2/64_Np64/Optimet3DP