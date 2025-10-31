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

all:	dense_integer dense_real dense_complex band_integer band_real band_complex

clean:
	rm -f *.mod *.o dense_integer dense_real dense_complex band_integer band_real band_complex *.dat *.plt

dense_integer:dense_integer.o $(OBJ)
	$(F90) $(FFLAGS) dense_integer.o $(OBJ) -o dense_integer

dense_real:dense_real.o $(OBJ)
	$(F90) $(FFLAGS) dense_real.o $(OBJ) -o dense_real

dense_complex:dense_complex.o $(OBJ)
	$(F90) $(FFLAGS) dense_complex.o $(OBJ) -o dense_complex

band_integer:band_integer.o $(OBJ)
	$(F90) $(FFLAGS) band_integer.o $(OBJ) -o band_integer

band_real:band_real.o $(OBJ)
	$(F90) $(FFLAGS) band_real.o $(OBJ) -o band_real

band_complex:band_complex.o $(OBJ)
	$(F90) $(FFLAGS) band_complex.o $(OBJ) -o band_complex

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

dense_integer.o:dense_integer.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  dense_integer.f90

dense_real.o:dense_real.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  dense_real.f90

dense_complex.o:dense_complex.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  dense_complex.f90

band_integer.o:band_integer.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  band_integer.f90

band_real.o:band_real.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  band_real.f90

band_complex.o:band_complex.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  band_complex.f90
