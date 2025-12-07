AC_PREREQ([2.69])
AC_INIT([gromtools],[0.1.0],[you@example.com])
AC_LANG([C++])

dnl ---- R + C++ toolchain ----
: ${R_HOME=`R RHOME`}
PATH_TO_CPP_COMPILER=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX([$PATH_TO_CPP_COMPILER])
AC_MSG_NOTICE([C++ compiler: $CXX])

dnl ---- Hardcoded PLINK2 root ----
PLINK2_ROOT="/sc/arion/projects/va-biobank/PROJECTS/libraries/plink-ng/2.0"
AC_MSG_NOTICE([Using PLINK2_ROOT = $PLINK2_ROOT])

dnl ---- Sanity check: pic_min must exist ----
if test ! -d "$PLINK2_ROOT/pic_min"; then
  AC_MSG_ERROR([
PLINK2_ROOT='$PLINK2_ROOT' does not contain a pic_min/ directory.

Please make sure this path is correct and that plink2 has been built
with pic_min objects there.
])
fi

AC_SUBST([PLINK2_ROOT])

dnl ---- Generate src/Makevars from src/Makevars.in ----
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
