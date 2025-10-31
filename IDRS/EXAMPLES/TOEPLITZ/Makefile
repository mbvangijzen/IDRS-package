#Standard gfortran, single processor. Should always work:
F90 = gfortran
FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=single

#gfortran, with co-arrays. Make sure that the communication library is installed and linked correctly.
#F90 = caf
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi

#ifort compiler with single node support for co-arrays:
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = precision_module.o dense_la_module.o matrix_module.o idrs_module.o

SRC = ../../MODULES
IDRS = $(SRC)/IDRS

all:	toeplitz

clean:  
	rm -f *.mod *.o toeplitz

toeplitz:toeplitz.o $(OBJ)
	$(F90) toeplitz.o $(OBJ) -o toeplitz

matrix_module.o:matrix_module.f90
	$(F90) -c  matrix_module.f90

precision_module.o:     $(IDRS)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/precision_module.f90

dense_la_module.o:      $(IDRS)/dense_la_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/dense_la_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs.inc $(IDRS)/msidrs.inc $(IDRS)/qmridr.inc $(IDRS)/msqmridr.inc matrix_module.o precision_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

toeplitz.o:toeplitz.f90 $(OBJ)
	$(F90) -c  toeplitz.f90
