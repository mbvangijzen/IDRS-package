#Standard gfortran, single processor. Should always work:
#F90 = gfortran
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=single

#gfortran, with co-arrays. Make sure that the communication library is installed and linked correctly.
F90 = caf
FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi

#ifort compiler with single node support for co-arrays:
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = interface_module.o precision_module.o dense_la_module.o user_module.o matrix_module.o idrs_module.o ritz_module.o pp_idrs_module.o
SRC = ../../MODULES
IDRS = $(SRC)/IDRS
MATRIX = $(SRC)/PP_MATRIX
INTERFACE = $(SRC)/INTERFACE
PP_IDRS = $(SRC)/PP_IDRS

all:	pagerank

clean:
	rm -f *.o *.mod pagerank *.dat *.plt

pagerank:pagerank.o $(OBJ)
	$(F90) $(FFLAGS) pagerank.o $(OBJ) -o pagerank

interface_module.o:  $(INTERFACE)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INTERFACE)/interface_module.f90

dense_la_module.o: $(IDRS)/dense_la_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/dense_la_module.f90

user_module.o:  $(MATRIX)/user_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/user_module.f90

matrix_module.o:$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/matrix_module.f90

precision_module.o:     $(IDRS)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/precision_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs.inc $(IDRS)/msidrs.inc $(IDRS)/qmridr.inc $(IDRS)/msqmridr.inc matrix_module.o precision_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

ritz_module.o:  $(IDRS)/ritz_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/ritz_module.f90

pp_idrs_module.o:	$(PP_IDRS)/pp_idrs_module.f90
	$(F90) $(FFLAGS) -c  $(PP_IDRS)/pp_idrs_module.f90

pagerank.o:pagerank.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  pagerank.f90
