#Standard gfortran, single processor. Should always work:
#F90 = gfortran
#FFLAGS = -O3 -std=f2018 -fcoarray=single

#gfortran, with co-arrays. Make sure that the communication library is installed and linked correctly.
#F90 = caf
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi

#ifort compiler with single node support for co-arrays:
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = interface_module.o precision_module.o dense_la_module.o mmio.o user_module.o mm_module.o matrix_module.o idrs_module.o recycle_module.o idrs_driver_module.o preconditioner_module.o

SRC = ../../MODULES
IDRS = $(SRC)/IDRS_ALGORITHMS
MATRIX = $(SRC)/MATRIX_OPERATIONS
INPUT = $(SRC)/INTERFACE
RECYCLE = $(SRC)/RECYCLE
DRIVER = $(SRC)/IDRS_DRIVER

all:	wedge_neumann wedge_sommerfeld room_neumann room_impedance stommel sag pagerank

clean:	
	rm -f *.mod *.o wedge_neumann wedge_sommerfeld room_neumann room_impedance stommel sag pagerank *.dat *.plt *.gif *.jpeg

wedge_neumann:	wedge_neumann.o $(OBJ)
	$(F90) $(FFLAGS) wedge_neumann.o $(OBJ) -o wedge_neumann

wedge_sommerfeld:	wedge_sommerfeld.o $(OBJ)
	$(F90) $(FFLAGS) wedge_sommerfeld.o $(OBJ) -o wedge_sommerfeld

room_neumann:	room_neumann.o $(OBJ)
	$(F90) $(FFLAGS) room_neumann.o $(OBJ) -o room_neumann

room_impedance:	room_impedance.o $(OBJ)
	$(F90) $(FFLAGS) room_impedance.o $(OBJ) -o room_impedance

stommel:	stommel.o $(OBJ)
	$(F90) $(FFLAGS) stommel.o $(OBJ) -o stommel

sag:	sag.o $(OBJ)
	$(F90) $(FFLAGS) sag.o $(OBJ) -o sag

pagerank:	pagerank.o $(OBJ)
	$(F90) $(FFLAGS) pagerank.o $(OBJ) -o pagerank

precision_module.o:	$(IDRS)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/precision_module.f90 

dense_la_module.o:	$(IDRS)/dense_la_module.f90 
	$(F90) $(FFLAGS) -c  $(IDRS)/dense_la_module.f90 

user_module.o:	$(MATRIX)/user_module.f90 
	$(F90) $(FFLAGS) -c  $(MATRIX)/user_module.f90 

mm_module.o:	mm_module.f90
	$(F90) $(FFLAGS) -c  mm_module.f90 

matrix_module.o:	$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c $(MATRIX)/matrix_module.f90 

preconditioner_module.o:	preconditioner_module.f90 
	$(F90) $(FFLAGS) -c preconditioner_module.f90 

interface_module.o:	$(INPUT)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INPUT)/interface_module.f90

idrs_module.o:	$(IDRS)/idrs_module.f90 $(IDRS)/idrs.inc $(IDRS)/qmridr.inc $(IDRS)/msidrs.inc $(IDRS)/msqmridr.inc
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

recycle_module.o:	$(RECYCLE)/recycle_module.f90
	$(F90) $(FFLAGS) -c  $(RECYCLE)/recycle_module.f90

idrs_driver_module.o:	$(DRIVER)/idrs_driver_module.f90
	$(F90) $(FFLAGS) -c  $(DRIVER)/idrs_driver_module.f90

wedge_sommerfeld.o:	wedge_sommerfeld.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  wedge_sommerfeld.f90

wedge_neumann.o:	wedge_neumann.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  wedge_neumann.f90

room_neumann.o:	room_neumann.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  room_neumann.f90

room_impedance.o:	room_impedance.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  room_impedance.f90

stommel.o:	stommel.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  stommel.f90

sag.o:	sag.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  sag.f90

pagerank.o:	pagerank.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  pagerank.f90

mmio.o:	mmio.f
	$(F90) -c  mmio.f
