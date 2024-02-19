F90 = gfortran 
#FFLAGS = -O3 -fcoarray=single -fcheck=all -std=f2018
FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = interface_module.o precision_module.o eigen_module.o user_module.o matrix_module.o idrs_module.o ritz_module.o
IDRS = ../IDRS_FORTRAN
MATRIX = ../MATRIX_MODULE
INPUT = ../INTERFACE

all:    room ms_room cdr ms_cdr cdr_time

clean:
	rm -f *.o *.mod room ms_room cdr ms_cdr cdr_time

room:room.o $(OBJ)
	$(F90) $(FFLAGS) room.o $(OBJ) -o room

ms_room:ms_room.o $(OBJ)
	$(F90) $(FFLAGS) ms_room.o $(OBJ) -o ms_room

cdr:cdr.o $(OBJ)
	$(F90) $(FFLAGS) cdr.o $(OBJ) -o cdr

ms_cdr:ms_cdr.o $(OBJ)
	$(F90) $(FFLAGS) ms_cdr.o $(OBJ) -o ms_cdr

cdr_time:cdr_time.o $(OBJ)
	$(F90) $(FFLAGS) cdr_time.o $(OBJ) -o cdr_time

precision_module.o:     $(MATRIX)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/precision_module.f90

eigen_module.o: $(IDRS)/eigen_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/eigen_module.f90

user_module.o:user_module.f90
	$(F90) $(FFLAGS) -c  user_module.f90

matrix_module.o:	$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/matrix_module.f90

interface_module.o:$(INPUT)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INPUT)/interface_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msqmridr_body.f90 
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

ritz_module.o:	$(IDRS)/ritz_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/ritz_module.f90

room.o:room.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  room.f90

ms_room.o:ms_room.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  ms_room.f90

cdr.o:cdr.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  cdr.f90

ms_cdr.o:ms_cdr.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  ms_cdr.f90

cdr_time.o:	cdr_time.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  cdr_time.f90
