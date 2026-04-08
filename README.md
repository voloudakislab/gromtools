# gromtools

`gromtools` is an R package for imputing genetically regulated omics values from PLINK2 genotype data.

In practice, the package:

- reads chromosome-split `.pgen` genotype files
- converts variant weights into an internal sparse binary format
- computes predicted omics values with additive linear SNP effects
- writes the result to an on-disk `.grom` matrix with matching sample `.sid` and feature `.gid` index files
- reads `.grom` matrix with optional subsetting of samples and features

## Installation

From the package root:

```bash
git clone voloudakislab/gromtools
cd gromtools
R CMD INSTALL .
```

## Notes

- The package requires R, Rcpp, a working C++17 toolchain, GNU `make` and R package `data.table`.
- `configure` can optionally detect Intel MKL through `MKLROOT`, but it can also build without MKL.
- Example input files are available under `inst/extdata/`.
