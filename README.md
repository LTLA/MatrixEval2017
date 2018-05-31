# Evaluating the performance of different matrix types

## Overview 

This repository contains code for the paper **beachmat: a Bioconductor C++ API for accessing high-throughput biological data from a variety of R matrix types** by Lun _et al._ 
(2018, _PLoS Computational Biology_; https://doi.org/10.1371/journal.pcbi.1006135).

The provided code will check the performance of different matrix types for row/column access, using simulated and real data sets.
To run the tests on your machine, please read the following instructions.

## Setup

1. Install [_beachmat_](https://bioconductor.org/packages/beachmat) from Bioconductor.
2. Enter `timings` and run `R CMD INSTALL --clean package`.
This requires installation of _RcppArmadillo_ and _RcppEigen_. 

## Running simulations
  
- `timings/` contains scripts for timings (in milliseconds) for accessing data from different matrix representations.
- `timings/chunking/` contains scripts for timing rechunking, as well as checking the chunk cache logic.
- `memory/` contains scripts for memory usage for different matrix representations.
- `miscellaneous` contains scripts to compare timings to R, and to verify the no-copy access method of _RcppArmadillo_ and _RcppEigen_.

## Running real analyses

### Zeisel dataset

Enter `real/zeisel` and download the [count matrix for the Zeisel data set](https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt).

- Execute the `zeisel_time.R` script to generate timings (in milliseconds) for matrix access to this data.
This will also determine memory usage for each matrix representation.
- Execute the `detection_stats.R` script to generate timings (in milliseconds) for computing various cell- or gene-based statistics from this data.

### 10X dataset

Enter `real/10X` and install [_TENxBrainData_](https://bioconductor.org/packages/TENxBrainData).
Read the `README.md` file for order of evaluation of the various Rmarkdown scripts.
