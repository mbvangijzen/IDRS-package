F90=gfortran 
#FFLAGS = -O3 -fcoarray=single -fcheck=all -std=f2018
FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = interface_module.o precision_module.o eigen_module.o user_module.o matrix_module.o idrs_module.o ritz_module.o
IDRS = ../IDRS_FORTRAN
MATRIX = ../MATRIX_MODULE
INTERFACE = ../INTERFACE

all:	pagerank ms_pagerank

clean:
	rm -f *.o *.mod pagerank ms_pagerank

pagerank:pagerank.o $(OBJ)
	$(F90) $(FFLAGS) pagerank.o $(OBJ) -o pagerank

ms_pagerank:ms_pagerank.o $(OBJ)
	$(F90) $(FFLAGS) ms_pagerank.o $(OBJ) -o ms_pagerank

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

pagerank.o:pagerank.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  pagerank.f90

ms_pagerank.o:ms_pagerank.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  ms_pagerank.f90
