# IDRS
## What is IDRS?
IDRS is a package for solving large and sparse nonsymmetric linear systems fo equations. The "IDR" in the names denotes Induced Dimension Reduction. The "S" stand for the iteration parameter s. This parameter determines the amount of vector operations and storage required by the method. A larger s normally increases the rate of convergence, but comes at the price of more memory usage and computations per iteration. For most problems a small value of s (for exampel 2 or 4) gives the best compromise. 
IDR-methods are Krylov subspace methods. They use a recursion depth of s+2. This in contrast to most popular Krylov methods for solvinf nonsymmetric linear systems: BiCGSTAB and GMRES. BiCGSTAB uses three-term recursions while GMRES uses long recursions (recursion depth is equal to the iteration number plus one). 
IDR-methods terminate using at most N+N/s iterations, with N the problem size. By comparison, BiCGSTAB terminates within 2N iterations and GMRES within N iterations.
## In what language is IDRS written?
IDRS comes in a Matlab version and a Fortran version. The Fortran version adheres to the 2018 standard. It is stand-alone, no external library is needed. Co-array Fortran is used for parallelisation. Co-array Fortran is part of the Fortran standard, so the package should run on any platform without modifications. However, some compilers (for example gfortran) require an external library for the programs to run in parallel.
## Which IDR algorithms are included in the package?
The package has implementations of the following algorithms:
- IDR(s): an efficient algorithm that uses bi-orthogonality conditions to compute new iteration vectors. Convergence may be erratic. The algorithms is optimised for parallel processing and has only one synchronisation point per iteration. Spectral information (ritz values and vectors) can be obtained, and an initial search space can be provided to speed-up the convergence. This is useful when solving sequences of linear systems.
  With s=1 a method that is mathematically equivalent to BiCGSTAB can be obtained. This IDR(1)-version of BiCGSTAB is exactly as efficient in terms of memory and computations as BiCGSTAB. However, BiCGSTAB might seem faster if one only checks the number of iterations. This is only cosmetic, since one BiCGSTAB iteration corresponds to two IDR iterations. 
- QMRIDR(s): a robust algorithm that uses orthogonalisation conditions to compute new iteration vectors. QMRIDR(s) uses "quasy" optimisation of the residual norm: a procedure that mimics the optimal minimisation of the residual norm that is used in GMRES. As long as the iteration number is smaller than s, QMRIDR(s) is equivalent to GMRES. However, QMRIDR(s) has extra overhead in memory and vector operations. Once the iterations number is larger than s, QMRIDR(s) can loosly be seen as restarted GMRES augmented with a restarting procedure that ensures finite termination.
  QMRIDR(s) is a flexible method, and can be used as an inner-outer iterative method. In the IDRS package, IDR(s) can be used as inner iterative method for performing preconditioning operations.
