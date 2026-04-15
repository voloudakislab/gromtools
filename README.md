<br>

### Computation Tools for Genetically Regulated Omics with PGEN Streaming

The `gromtools` package enables high-throughput imputation of genetically
regulated omics values from chromosome-split PLINK2 genotype data using
additive linear SNP models.

Major functionality of the `gromtools` package:

+ [`read_db_dir()`](reference/read_db_dir.html)                Load model SQLite databases into a weights table
+ [`gromtools_impute()`](reference/gromtools_impute.html)  Build sparse weights and write `.grom`, `.gid`, and `.sid` outputs
+ [`gromtools_read()`](reference/gromtools_read.html)        Stream selected model, gene, and sample combinations from an existing `.grom` output

## Resources

+ Minimal [getting started vignette](articles/getting-started.html)

## Motivation

Large-scale genetically regulated omics workflows need a compact interface for
loading prediction weights, streaming genotype data, and writing output in a
format that can be read back efficiently without materializing the full matrix
in memory. `gromtools` is built around that workflow. It converts model weights
into a sparse binary representation, streams chromosome-split `.pgen` inputs,
computes predicted values using additive SNP effects, and writes the result to
an on-disk `.grom` matrix with matching `.gid` and `.sid` index files.

The package is designed around three user-facing steps: read model databases,
run imputation, and read back selected results. Example synthetic inputs for the
full workflow are bundled with the package.

## Install

```r
devtools::install_github("voloudakislab/gromtools")
```

## Technical notes

The package requires R, Rcpp, a working C++17 toolchain, GNU `make`, and
`data.table`. `configure` can optionally detect Intel MKL through `MKLROOT`,
but it can also build without MKL.
