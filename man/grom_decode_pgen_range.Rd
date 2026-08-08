\name{grom_decode_pgen_range}
\alias{grom_decode_pgen_range}
\alias{grom_decode_pgen_range_readlist}
\title{Decode a Contiguous PGEN Variant Range}
\usage{
grom_decode_pgen_range(
  pgen_file,
  raw_sample_ct,
  first_variant_1based,
  variant_count,
  sample_subset,
  meanimpute = TRUE,
  rounded_mean = FALSE,
  pgen_threads = 1L
)

grom_decode_pgen_range_readlist(
  pgen_file,
  raw_sample_ct,
  first_variant_1based,
  variant_count,
  sample_subset,
  meanimpute = TRUE,
  rounded_mean = FALSE
)
}
\arguments{
\item{pgen_file}{Path to a PLINK2 \code{.pgen} file.}
\item{raw_sample_ct}{Total number of samples in the PGEN file.}
\item{first_variant_1based}{First variant to decode, using 1-based PGEN indexing.}
\item{variant_count}{Number of contiguous variants to decode.}
\item{sample_subset}{Optional strictly increasing 1-based sample subset. Use \code{integer(0)} for all samples.}
\item{meanimpute}{Logical indicating whether missing dosages should be mean imputed.}
\item{rounded_mean}{Logical indicating whether imputed means should be rounded.}
\item{pgen_threads}{Number of independent PGEN loader threads requested by the direct range decoder.}
}
\description{
Decode a contiguous PGEN variant range into a column-major numeric matrix.
\code{grom_decode_pgen_range()} uses the direct range decoder used by the
parallel PGEN loader. \code{grom_decode_pgen_range_readlist()} uses the legacy
\code{ReadList()} path and is intended for regression testing.
}
\value{
A numeric matrix with samples in rows and variants in columns.
}
