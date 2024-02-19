F90 = gfortran
#FFLAGS = -O3 -fcoarray=single -fcheck=all -std=f2018
FFLAGS = -O3 -fcheck=all -std=f2018 -fimplicit-none -fcoarray=lib -I/opt/homebrew/include -L/opt/homebrew/lib -lcaf_mpi 
#FFLAGS = -O3 -fcheck=all -std=f2018 -fcoarray=lib -L${OPENCOARRAYS_ROOT}/lib64 -lcaf_mpi
#F90=ifort
#FFLAGS = -Ofast -coarray

OBJ = interface_module.o precision_module.o eigen_module.o mmio.o user_module.o mm_module.o matrix_module.o idrs_module.o ritz_module.o
IDRS = ../IDRS_FORTRAN
MATRIX = ../MATRIX_MODULE
INPUT = ../INTERFACE

all:	wedge_sommerfeld wedge_neumann ms_wedge_sommerfeld ms_wedge_neumann stommel sag 

clean:	
	rm -f *.mod *.o wedge_sommerfeld wedge_neumann ms_wedge_sommerfeld ms_wedge_neumann stommel sag

wedge_sommerfeld:	wedge_sommerfeld.o $(OBJ)
	$(F90) $(FFLAGS) wedge_sommerfeld.o $(OBJ) -o wedge_sommerfeld

wedge_neumann:	wedge_neumann.o $(OBJ)
	$(F90) $(FFLAGS) wedge_neumann.o $(OBJ) -o wedge_neumann

ms_wedge_sommerfeld:	ms_wedge_sommerfeld.o $(OBJ)
	$(F90) $(FFLAGS) ms_wedge_sommerfeld.o $(OBJ) -o ms_wedge_sommerfeld

ms_wedge_neumann:	ms_wedge_neumann.o $(OBJ)
	$(F90) $(FFLAGS) ms_wedge_neumann.o $(OBJ) -o ms_wedge_neumann

stommel:	stommel.o $(OBJ)
	$(F90) $(FFLAGS) stommel.o $(OBJ) -o stommel

sag:	sag.o $(OBJ)
	$(F90) $(FFLAGS) sag.o $(OBJ) -o sag

precision_module.o:	$(MATRIX)/precision_module.f90
	$(F90) $(FFLAGS) -c  $(MATRIX)/precision_module.f90 

eigen_module.o:	$(IDRS)/eigen_module.f90 
	$(F90) $(FFLAGS) -c  $(IDRS)/eigen_module.f90 

user_module.o:	$(MATRIX)/user_module.f90 
	$(F90) $(FFLAGS) -c  $(MATRIX)/user_module.f90 

mm_module.o:	mm_module.f90
	$(F90) $(FFLAGS) -c  mm_module.f90 

matrix_module.o:	$(MATRIX)/matrix_module.f90
	$(F90) $(FFLAGS) -c $(MATRIX)/matrix_module.f90 

interface_module.o:	$(INPUT)/interface_module.f90
	$(F90) $(FFLAGS) -c  $(INPUT)/interface_module.f90

idrs_module.o:	$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/msqmridr_body.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

ritz_module.o:	$(IDRS)/ritz_module.f90
	$(F90) $(FFLAGS) -c  $(IDRS)/ritz_module.f90

wedge_sommerfeld.o:	wedge_sommerfeld.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  wedge_sommerfeld.f90

wedge_neumann.o:	wedge_neumann.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  wedge_neumann.f90

ms_wedge_sommerfeld.o:	ms_wedge_sommerfeld.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  ms_wedge_sommerfeld.f90

ms_wedge_neumann.o:	ms_wedge_neumann.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  ms_wedge_neumann.f90

stommel.o:	stommel.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  stommel.f90

sag.o:	sag.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  sag.f90

mmio.o:	mmio.f
	$(F90) -c  mmio.f
