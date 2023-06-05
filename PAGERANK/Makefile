F90=gfortran -O -fcoarray=single
OBJ = idrs_module.o matrix_module.o preconditioner_module.o
IDRS = ../IDRS_FORTRAN

pagerank:pagerank.o $(OBJ)
	$(F90) pagerank.o $(OBJ) -o pagerank

matrix_module.o:matrix_module.f90
	$(F90) -c  matrix_module.f90

preconditioner_module.o:preconditioner_module.f90
	$(F90) -c  preconditioner_module.f90

idrs_module.o:$(IDRS)/idrs_module.f90 $(IDRS)/idrs_body.f90 $(IDRS)/msidrs_body.f90 $(IDRS)/qmridr_body.f90 $(IDRS)/msqmridr_body.f90 matrix_module.o preconditioner_module.o
	$(F90) -c  $(IDRS)/idrs_module.f90

pagerank.o:pagerank.f90 $(OBJ)
	$(F90) -c  pagerank.f90
