\name{grom_pgen_test_helpers}
\alias{grom_test_readlist_buffer_contract}
\alias{grom_test_pgen_workspace_estimate}
\title{Internal PGEN Regression Test Helpers}
\usage{
grom_test_readlist_buffer_contract(
  pgen_file,
  raw_sample_ct,
  first_variant_1based,
  variant_count,
  sample_subset,
  meanimpute = TRUE,
  rounded_mean = FALSE
)

grom_test_pgen_workspace_estimate(
  pgen_file,
  raw_sample_ct,
  sample_subset,
  variant_1based,
  rounded_mean = FALSE
)
}
\arguments{
\item{pgen_file}{Path to a PGEN file.}
\item{raw_sample_ct}{Number of samples in the PGEN file.}
\item{first_variant_1based}{First variant to decode, using one-based indexing.}
\item{variant_count}{Number of contiguous variants to decode.}
\item{sample_subset}{Optional one-based sample subset.}
\item{meanimpute}{Whether missing genotypes are mean-imputed.}
\item{rounded_mean}{Whether imputed means are rounded.}
\item{variant_1based}{One-based variant index used for the workspace test.}
}
\description{
Internal helpers used by the package test suite to verify low-level PGEN
reader buffer validation and workspace-memory diagnostics.
}
\keyword{internal}
