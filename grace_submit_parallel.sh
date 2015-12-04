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
#$ -N OPTIMET3DP

# 6. Select the QLogic parallel environment (qlc) and 16 process.
#$ -pe qlc 64

# 7. Select the project that this job will run under.
# Find <your_project_id> by running the command "groups"
# -P OPTIMETPROJECT

# 8. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/ucakgim/Scratch

# Set up teh desired module environment
module purge
module load rcps-core
# compilers/MKL
module load gcc-libs compilers/intel/2015/update2
# MPI
module load mpi/intel/2015/update3/intel
# GSL
module load gsl
# HDF5
module load hdf/5-1.8.15-p1-impi
# F2C
module load personal-modules
module load f2c

# 9. Run the application.
gerun $HOME/git/OPTIMET/build/Optimet3DP
