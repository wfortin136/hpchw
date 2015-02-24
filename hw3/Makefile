MPIF =	mpifort
FORTRAN = gfortran

all: serial block_parallel noblock_parallel openmp

serial:serial_laxwendroff.f95
	$(FORTRAN) -o serial serial_laxwendroff.f95
openmp: openmp_laxwendroff.f95
	$(FORTRAN) -o openmp -fopenmp openmp_laxwendroff.f95

block_parallel:parallel_laxwendroff.f95
	$(MPIF) -o block_parallel parallel_laxwendroff.f95	

noblock_parallel:noblock_parallel_laxwendroff.f95
	$(MPIF) -o noblock_parallel noblock_parallel_laxwendroff.f95
clean:
	rm -f *.a *.o a.out core* serial block_parallel noblock_parallel openmp
