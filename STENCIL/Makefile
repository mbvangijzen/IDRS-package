F90=gfortran
NFDIR = /opt/homebrew
OBJ = idrs_module.o matrix_module.o preconditioner_module.o general_module.o
IDRS = ../IDRS_FORTRAN
FFLAGS = -O -fcoarray=lib -I${NFDIR}/include -L${NFDIR}/lib -lcaf_mpi

stencil:stencil.o $(OBJ)
	$(F90) $(FFLAGS) stencil.o $(OBJ) -o stencil

matrix_module.o:matrix_module.f90
	$(F90) $(FFLAGS) -c  matrix_module.f90

preconditioner_module.o:preconditioner_module.f90
	$(F90) $(FFLAGS) -c  preconditioner_module.f90

general_module.o:general_module.f90
	$(F90) $(FFLAGS) -c  general_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msqmridr_body.f90 matrix_module.o preconditioner_module.o
	$(F90) $(FFLAGS) -c  $(IDRS)/idrs_module.f90

stencil.o:stencil.f90 $(OBJ)
	$(F90) $(FFLAGS) -c  stencil.f90
