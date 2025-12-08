library(Rcpp)
Rcpp::Rcpp.package.skeleton("gromtools")
Rcpp::compileAttributes("/sc/arion/projects/va-biobank/PROJECTS/gromtools/gromtools")
Rcpp::compileAttributes("/sc/arion/projects/va-biobank/PROJECTS/gromtools/gromtools")


.libPaths()
# See where pgenlibr is found on each lib path (if at all)
for (p in .libPaths()) {
  cat("\nLib path:", p, "\n")
  print(system.file(package = "pgenlibr", lib.loc = p))
}
library(GenomicDataStream)
install.packages("GenomicDataStream")
BiocManager::install("GabrielHoffman/GenomicDataStream")
