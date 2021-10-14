#!/bin/bash -l

#Batch script to run an MPI parallel job with the upgraded software
# stack under SGE with Intel MPI.

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=23:10:0

# 3. Request 1 gigabyte of RAM per process (must be an integer)
#$ -l mem=7G

module unload default-modules/2018

module load gerun

module remove compilers mpi

module load compilers/gnu/4.9.2 mpi/openmpi/3.1.4/gnu-4.9.2 hdf/5-1.10.5/gnu-4.9.2 gsl/1.16/gnu-4.9.2
module load openblas/0.3.7-serial/gnu-4.9.2 scalapack/2.0.2/gnu-4.9.2/openblas-0.3.7 optimet/1.0.1 f2c

module list

cd /home/uceeise/Scratch/OPTIMET/build_belos
make

# 5. Set the name of the job.
#$ -N ManySpheres4

# 6. Select the MPI parallel environment and number of processes.
#$ -pe mpi 120

# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.

#$ -wd /home/uceeise/Scratch/OPTIMET/examples/output

# 8. Run our MPI job.  GERun is a wrapper that launches MPI jobs on our clusters.
cd /home/uceeise/Scratch/OPTIMET/examples
gerun /home/uceeise/Scratch/OPTIMET/build_belos/Optimet3D ManySpheres4.xml | grep "e-"

