F90 = gfortran
IDRS = ../../IDRS_FORTRAN
OBJ = mmio.o matrix_module.o preconditioner_module.o general_module.o idrs_module.o
FFLAGS = -O -fcoarray=single -fcheck=all

all:	mm_idrs mm_msidrs

mm_idrs:	mm_idrs.o $(OBJ)
	$(F90) $(FFLAGS) mm_idrs.o $(OBJ) -o mm_idrs

mm_msidrs:	mm_msidrs.o $(OBJ)
	$(F90) $(FFLAGS) mm_msidrs.o $(OBJ) -o mm_msidrs

matrix_module.o:	matrix_module.f90
	$(F90) $(FFLAGS) -c  matrix_module.f90 

preconditioner_module.o:	preconditioner_module.f90
	$(F90) $(FFLAGS) -c  preconditioner_module.f90 

general_module.o:	general_module.f90
	$(F90) $(FFLAGS) -c  general_module.f90

idrs_module.o:	$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/msqmridr_body.f90 matrix_module.o preconditioner_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

mm_idrs.o:	mm_idrs.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  mm_idrs.f90

mm_msidrs.o:	mm_msidrs.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  mm_msidrs.f90

mmio.o:	mmio.f
	$(F90) -c  mmio.f
