# install_gromtools.sh
#!/bin/bash
set -e
cd /sc/arion/projects/va-biobank/PROJECTS/gromtools/gromtools
R CMD build .
R CMD INSTALL gromtools_0.1.0.tar.gz -l ~/R/libs
