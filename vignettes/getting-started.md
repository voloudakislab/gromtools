Getting started with gromtools
================

## Load the package

``` r
library(gromtools)
```

## Declare input paths

Declare PLINK2 chromosomes directory

``` r
pgen_dir <- system.file("extdata", "synthetic_chromosomes", package = "gromtools")
print(list.files(pgen_dir)[1:12])
#>  [1] "synth.chr1.pgen"  "synth.chr1.psam"  "synth.chr1.pvar"  "synth.chr10.pgen"
#>  [5] "synth.chr10.psam" "synth.chr10.pvar" "synth.chr11.pgen" "synth.chr11.psam"
#>  [9] "synth.chr11.pvar" "synth.chr12.pgen" "synth.chr12.psam" "synth.chr12.pvar"
```

## Read model weights SQL db files

``` r
db_directory <- system.file("extdata", "synth_small_variant_weights_db", package = "gromtools")
list.files(db_directory)
#> [1] "AMR_subclass_IN_SST.db" "AMR_subclass_VLMC.db"

# Loads db files in weights table format
model_weights_table <- read_db_dir(db_dir = db_directory)
head(model_weights_table)
#>               model_ID            gene        rsid chromosomes position
#> 1: AMR_subclass_IN_SST ENSG00000228594  rs28625089        chr1   904947
#> 2: AMR_subclass_IN_SST ENSG00000228594  rs12080505        chr1  1110226
#> 3: AMR_subclass_IN_SST ENSG00000228594   rs6604971        chr1  1135812
#> 4: AMR_subclass_IN_SST ENSG00000228594  rs13374146        chr1  1180851
#> 5: AMR_subclass_IN_SST ENSG00000228594   rs4970443        chr1  1381507
#> 6: AMR_subclass_IN_SST ENSG00000228594 rs368836064        chr1  1382885
#>    ref_allele eff_allele      weight
#> 1:          G          A -0.07537545
#> 2:          A          C  0.14694625
#> 3:          T          C  0.08669795
#> 4:          T          C -0.09174421
#> 5:          G          A -0.03473008
#> 6:          A          G -0.03415348
```

## Run imputation

Use `gromtools_impute()` to impute genetically regulated -omics levels
such as gene expression.

``` r

# Output files get written at `tmp_grom_run/synth_example.grom/gid/sid`
grom_pfx="tmp_grom_run/synth_example"

gromtools_impute(
  grom_pfx,
  weights_table = model_weights_table,
  pgen_dir = pgen_dir
)
```

## Read .grom file

Use `gromtools_read()` to stream selected model / gene / individual
combinations from an existing `.grom` output.

``` r
# models, genes, and samples default to NULL; NULL loads all available entries.
grom_mat <- gromtools_read(
  grom_pfx = grom_pfx,
  models = c("AMR_subclass_IN_SST", "AMR_subclass_VLMC"), # defaults to NULL; selects all models
  genes = c("ENSG00000003987", "ENSG00000053900"),        # defaults to NULL; selects all genes
  samples = c(                                            # defaults to NULL; selects all samples
    "sample0001", "sample0014", "sample0034", "sample0039", "sample0043"
  )
)


print(grom_mat)
#>            AMR_subclass_IN_SST_ENSG00000003987
#> sample0001                           0.3191238
#> sample0014                           0.4542425
#> sample0034                           0.4138578
#> sample0039                           0.2676161
#> sample0043                           0.3721308
#>            AMR_subclass_IN_SST_ENSG00000053900
#> sample0001                           0.2737647
#> sample0014                          -0.4481538
#> sample0034                          -0.2212396
#> sample0039                          -0.3217697
#> sample0043                           0.2609762
#>            AMR_subclass_VLMC_ENSG00000053900
#> sample0001                        -0.5058091
#> sample0014                        -0.9249124
#> sample0034                        -0.5259653
#> sample0039                        -0.9044481
#> sample0043                        -0.6188168
```
