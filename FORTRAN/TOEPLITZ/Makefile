F90=gfortran -O3 -fcoarray=single -fcheck=all -std=f2018
OBJ = precision_module.o matrix_module.o idrs_module.o
IDRS = ../IDRS_FORTRAN
MATRIX = ../MATRIX_MODULE

toeplitz:toeplitz.o $(OBJ)
	$(F90) toeplitz.o $(OBJ) -o toeplitz

matrix_module.o:matrix_module.f90
	$(F90) -c  matrix_module.f90

precision_module.o:     $(MATRIX)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/precision_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msqmridr_body.f90 matrix_module.o precision_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

toeplitz.o:toeplitz.f90 $(OBJ)
	$(F90) -c  toeplitz.f90
