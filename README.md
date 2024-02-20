# IDRS
## What is IDRS?
IDRS is a package for solving large and sparse nonsymmetric linear systems fo equations. The "IDR" in the names denotes Induced Dimension Reduction. The "S" stand for the iteration parameter s. This parameter determines the amount of vector operations and storage required by the method. A larger s normally increases the rate of convergence, but comes at the price of more memory usage and computations per iteration. For most problems a small value of s (for exampel 2 or 4) gives the best compromise. 
IDR-methods use a recursion depth of s+2. This in contrast to the well known three-term recursion method BiCGSTAB and the long recursion method GMRES. For special choices of the parameters IDR-methods can be obtained that are equivalent to BiCGSTAB and GMRES.
IDR-methods terminate using at most N+N/s iterations, with N the problem size. By comparison, BiCGSTAB terminates within 2N iterations and GMRES within N iterations.
