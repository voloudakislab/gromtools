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

Use `gromtools_impute()` to impute genetically regulated -omics levels (e.g. gene expression).

The weights table must contain `model_ID`, `gene`, `chromosomes`, `rsid`, and `weight`. An `ancestry` column is optional and, when present, is preserved in the generated manifest metadata.

```r
library(gromtools)

pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
db_dir <- system.file("extdata", "synth_small_variant_weights_db", package = "gromtools")
weights_table <- read_db_dir(db_dir)
out_dir <- file.path("tmp_grom_run")
grom_pfx <- file.path(out_dir, "synth_example")

gromtools_impute(
  weights_table = weights_table,
  grom_pfx = grom_pfx,
  pgen_dir = pgen_dir
)
```

This writes output files under `tmp_grom_run/`.

## Read grom output

Use `gromtools_read()` to stream selected model / gene / individual combinations from an existing `.grom` output.

```r
library(gromtools)

# Change this path if you want to point to a different grom output.
grom_pfx <- "tmp_grom_run/synth_example"

# models, genes, and samples default to NULL; NULL loads all available entries.
grom_mat <- gromtools_read(
  grom_pfx = grom_pfx,
  models = c("AMR_subclass_IN_SST", "AMR_subclass_VLMC"),
  genes = c("ENSG00000003987", "ENSG00000053900"),
  samples = c(
    "sample0001", "sample0014", "sample0034", "sample0039", "sample0043",
    "sample0087", "sample0082", "sample0068", "sample0059", "sample0051"
  )
)

cat("\nSmall demo matrix from the .grom file:\n")
print(grom_mat)
```

Returned matrix columns are named as `model_ID_gene`, so the same gene can be loaded from multiple models without ambiguous duplicate column names.

## Notes

- The package requires R, Rcpp, a working C++17 toolchain, and GNU `make`.
- `configure` can optionally detect Intel MKL through `MKLROOT`, but it can also build without MKL.
- Example input files are available under `inst/extdata/`.

## Troubleshooting

If `R CMD build .` or `R CMD INSTALL .` fails with an error like `./configure: not found`, the `configure` script may have Windows CRLF line endings. Run:

```bash
sed -i 's/\r$//' configure
```
