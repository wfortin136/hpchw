#!/bin/bash
#SBATCH --job-name=insitu

#SBATCH --nodes=32
#SBATCH --ntasks-per-node=1
#SBATCH --output=midinsitu.out
# this job requests exclusive access to the nodes it is given
# this mean it will be the only job running on the node

# --constraint=ib must be give to guarantee a job is allocated 
# nodes with Infiniband
#SBATCH --constraint=ib

# load the openmpi module
module load mvapich2

# Run the process with mpirun. Notice -n is not required. mpirun will
# automatically figure out how many processes to run from the slurm option
mpirun -N 32 ./noblock_parallel inputfile.in
mpirun -N 16 ./noblock_parallel inputfile.in
mpirun -N 8 ./noblock_parallel inputfile.in
mpirun -N 4 ./noblock_parallel inputfile.in
mpirun -N 1 ./noblock_parallel inputfile.in
mpirun -N 32 ./block_parallel inputfile.in
mpirun -N 16 ./noblock_parallel inputfile.in
mpirun -N 16 ./block_parallel inputfile.in
mpirun -N 8 ./block_parallel inputfile.in
mpirun -N 4 ./block_parallel inputfile.in
mpirun -N 1 ./block_parallel inputfile.in
