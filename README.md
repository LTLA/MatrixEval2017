# Evaluating the performance of different matrix types

This checks the performance of different matrix types for row/column access, using simulated and real data sets.
To run the tests on your machine, follow these instructions.

1. Install [_beachmat_](https://github.com/LTLA/beachmat).
This requires [_Rhdf5lib_](https://github.com/grimbough/Rhdf5lib) and _Rcpp_.
2. Enter `timings` and run `R CMD INSTALL --clean package`.
This requires _RcppArmadillo_ in addition to _beachmat_.
3. Enter `simulations` and execute all R scripts.
This may take some time and will generate timings (in milliseconds) for matrix access in a range of scenarios.
4. Enter `real` and download the count matrix for the Zeisel data set ([here](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt)).
Execute the R script to generate timings (in milliseconds) for matrix access to this data.
5. Enter `memory` and execute the R script to determine memory usage for each matrix representation.

