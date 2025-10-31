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

OBJ = interface_module.o precision_module.o dense_la_module.o user_module.o matrix_module.o idrs_module.o ritz_module.o pp_idrs_module.o
SRC = ../../MODULES
IDRS = $(SRC)/IDRS
MATRIX = $(SRC)/PP_MATRIX
INPUT = $(SRC)/INTERFACE
PP_IDRS = $(SRC)/PP_IDRS

all:    top_opt

clean:
	rm -f *.o *.mod top_opt *.dat *.plt *.gif *.jpeg

top_opt:top_opt.o $(OBJ)
	$(F90) $(FFLAGS) top_opt.o $(OBJ) -o top_opt

precision_module.o:     $(IDRS)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/precision_module.f90

dense_la_module.o: $(IDRS)/dense_la_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/dense_la_module.f90

user_module.o:user_module.f90
	$(F90) $(FFLAGS) -c  user_module.f90

matrix_module.o:	$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/matrix_module.f90

interface_module.o:$(INPUT)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INPUT)/interface_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs.inc $(IDRS)/msidrs.inc $(IDRS)/qmridr.inc $(IDRS)/msqmridr.inc 
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

ritz_module.o:	$(IDRS)/ritz_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/ritz_module.f90

pp_idrs_module.o:	$(PP_IDRS)/pp_idrs_module.f90
	$(F90) $(FFLAGS) -c  $(PP_IDRS)/pp_idrs_module.f90

top_opt.o:top_opt.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  top_opt.f90
