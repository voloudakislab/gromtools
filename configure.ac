AC_PREREQ([2.69])
AC_INIT([gromtools],[0.1.2],[you@example.com])
AC_LANG([C++])

dnl ---- R + C++ toolchain ----
: ${R_HOME=`R RHOME`}
PATH_TO_CPP_COMPILER=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX([$PATH_TO_CPP_COMPILER])
AC_MSG_NOTICE([C++ compiler: $CXX])

dnl ---- Optional MKL detection ----
dnl User can set MKLROOT in the environment, e.g.
dnl   MKLROOT=/hpc/packages/.../mkl/2025.0 R CMD INSTALL ...
AC_ARG_VAR([MKLROOT],[Root directory of Intel MKL (optional)])

HAVE_MKL=no
MKL_CPPFLAGS=""
MKL_LIBS=""

if test "x$MKLROOT" != x; then
  AC_MSG_NOTICE([Checking for MKL under $MKLROOT])
  if test -f "$MKLROOT/include/mkl.h"; then
    HAVE_MKL=yes
    MKL_CPPFLAGS="-I$MKLROOT/include"

    dnl Simple libdir guess: prefer $MKLROOT/lib if it exists
    if test -d "$MKLROOT/lib"; then
      MKL_LIBDIR="$MKLROOT/lib"
    elif test -d "$MKLROOT/lib/intel64"; then
      MKL_LIBDIR="$MKLROOT/lib/intel64"
    else
      MKL_LIBDIR="$MKLROOT/lib"
    fi

    MKL_LIBS="-L$MKL_LIBDIR -lmkl_rt -lpthread -lm -ldl"
  else
    AC_MSG_WARN([mkl.h not found under $MKLROOT; continuing without MKL])
  fi
else
  AC_MSG_NOTICE([MKLROOT not set; building without MKL])
fi

if test "x$HAVE_MKL" = xyes; then
  AC_DEFINE([HAVE_MKL],[1],[Define to 1 if MKL is available])
  AC_MSG_NOTICE([MKL detected; will use MKL])
else
  AC_MSG_NOTICE([MKL not detected; falling back to generic BLAS/LAPACK])
fi

AC_SUBST([MKL_CPPFLAGS])
AC_SUBST([MKL_LIBS])

dnl ---- Generate src/Makevars from src/Makevars.in ----
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
