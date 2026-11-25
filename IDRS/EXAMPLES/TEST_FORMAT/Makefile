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

OBJ = interface_module.o precision_module.o dense_la_module.o user_module.o matrix_module.o idrs_module.o recycle_module.o idrs_driver_module.o

SRC = ../../MODULES
IDRS = $(SRC)/IDRS_ALGORITHMS
MATRIX = $(SRC)/MATRIX_OPERATIONS
INTERFACE = $(SRC)/INTERFACE
RECYCLE = $(SRC)/RECYCLE
IDRS_DRIVER = $(SRC)/IDRS_DRIVER

all:	test_real test_complex test_integer

clean:
	rm -f *.mod *.o test_real test_complex test_integer *.dat *.plt

test_real:test_real.o $(OBJ)
	$(F90) $(FFLAGS) test_real.o $(OBJ) -o test_real

test_complex:test_complex.o $(OBJ)
	$(F90) $(FFLAGS) test_complex.o $(OBJ) -o test_complex

test_integer:test_integer.o $(OBJ)
	$(F90) $(FFLAGS) test_integer.o $(OBJ) -o test_integer

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

recycle_module.o:  $(RECYCLE)/recycle_module.f90
	$(F90) $(FFLAGS) -c  $(RECYCLE)/recycle_module.f90

idrs_driver_module.o:	$(IDRS_DRIVER)/idrs_driver_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS_DRIVER)/idrs_driver_module.f90

test_real.o:test_real.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  test_real.f90

test_complex.o:test_complex.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  test_complex.f90

test_integer.o:test_integer.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  test_integer.f90

