F90=gfortran 
FFLAGS = -O3 -fcoarray=single -fcheck=all -std=f2018
%FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
OBJ = interface_module.o precision_module.o eigen_module.o user_module.o matrix_module.o idrs_module.o ritz_module.o
IDRS = ../IDRS_FORTRAN
MATRIX = ../MATRIX_MODULE
INTERFACE = ../INTERFACE

all:	integer_random real_random complex_random

clean:
	rm -f *.mod *.o integer_random real_random complex_random

integer_random:integer_random.o $(OBJ)
	$(F90) $(FFLAGS) integer_random.o $(OBJ) -o integer_random

real_random:real_random.o $(OBJ)
	$(F90) $(FFLAGS) real_random.o $(OBJ) -o real_random

complex_random:complex_random.o $(OBJ)
	$(F90) $(FFLAGS) complex_random.o $(OBJ) -o complex_random

interface_module.o:  $(INTERFACE)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INTERFACE)/interface_module.f90

eigen_module.o: $(IDRS)/eigen_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/eigen_module.f90

user_module.o:  $(MATRIX)/user_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/user_module.f90

matrix_module.o:$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/matrix_module.f90

precision_module.o:     $(MATRIX)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/precision_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msqmridr_body.f90 matrix_module.o precision_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

ritz_module.o:  $(IDRS)/ritz_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/ritz_module.f90

integer_random.o:integer_random.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  integer_random.f90

real_random.o:real_random.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  real_random.f90

complex_random.o:complex_random.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  complex_random.f90
