# gromtools

`gromtools` is an R package for imputing genetically regulated omics values from PLINK2 genotype data.

In practice, the package:

- reads chromosome-split `.pgen` genotype files
- converts variant weights into an internal sparse binary format
- computes predicted omics values with additive linear SNP effects
- writes the result to an on-disk `.grom` matrix with matching sample and feature index files

## Installation

From the package root:

```bash
R CMD INSTALL .
```

## Inspect Input

If you want to inspect the example genotype input files and weights table:

```r
pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
weights_path <- system.file("extdata", "synth_small_variant_weights.tsv", package = "gromtools")

list.files(pgen_dir)
weights <- read_table(weights_path)
head(weights)
```

## Minimal Imputation Example

For imputation on synthetic data, after you install run:

```r
library(gromtools)

pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
weights_path <- system.file("extdata", "synth_small_variant_weights.tsv", package = "gromtools")
out_dir <- file.path("tmp_grom_run")
grom_pfx <- file.path(out_dir, "synth_example")

gromtools_impute(
  weights_path = weights_path,
  grom_pfx = grom_pfx,
  pgen_dir = pgen_dir
)
```

This writes output files under `tmp_grom_run/`.

## Read grom output


## Notes

- The package requires R, Rcpp, a working C++17 toolchain, and GNU `make`.
- `configure` can optionally detect Intel MKL through `MKLROOT`, but it can also build without MKL.
- Example input files are available under `inst/extdata/`.

## Troubleshooting

If `R CMD build .` or `R CMD INSTALL .` fails with an error like `./configure: not found`, the `configure` script may have Windows CRLF line endings. Run:

```bash
sed -i 's/\r$//' configure
```
