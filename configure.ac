AC_ARG_VAR( WARNINGS, keeps compiler warnings options )

AC_PREREQ([2.71])

AC_INIT([Spinc structures on real Bott manifolds],[0.8],[rafal.lutowski@ug.edu.pl])

# Checks for programs.
# AC_PROG_CXX
if test -z $CFLAGS; then
    CFLAGS='-Ofast -march=native'
fi
AC_PROG_CC
AC_PROG_CXX
AC_PROG_EGREP

AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
# AC_DEFINE([ENABLE_OPENMP],[1],[Compile operators with OpenMP support])
CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
# Oddly, linker needs same OpenMP flag(s) as compiler
# 20170909 Linux clang on travis does not resolve __kmpc* symbols, apparently needs -lgomp or -lomp
LDFLAGS="${LDFLAGS} ${OPENMP_CFLAGS}"

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T


# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([bzero clock_gettime gettimeofday])

AC_MSG_CHECKING([gcc warnings options])
AC_ARG_ENABLE(warnings,AS_HELP_STRING([--enable-warnings=wall|all|extra],[add various warnings options to CFLAGS [default=wall]]),enable_warnings=${enableval},enable_warnings=wall)
if ( test "x${enable_warnings}" = "xwall" ); then
    WARNINGS="-Wall"
    AC_MSG_RESULT([wall])
fi;
if ( test "x${enable_warnings}" = "xall" ); then
    WARNINGS="-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wuninitialized -Wstrict-prototypes"
    AC_MSG_RESULT([all])
fi;
if ( test "x${enable_warnings}" = "xextra" ); then
    WARNINGS="-Wall -Wextra -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wuninitialized -Wstrict-prototypes -Wconversion"
    AC_MSG_RESULT([extra])
fi;
AC_SUBST([WARNINGS])

AC_MSG_CHECKING([row width])
AC_ARG_ENABLE(short-row,AS_HELP_STRING([--enable-short-row],[use 32-bit row (default: 64-bit)]),ROW_WIDTH=32,ROW_WIDTH=64)
AC_MSG_RESULT([${ROW_WIDTH}])

AC_SUBST([ROW_WIDTH])

AC_MSG_CHECKING([progress enabled])
AC_ARG_ENABLE(progress, AS_HELP_STRING([--enable-progress],[give possibility to show progress of calculations]))
if ( test "x${enable_progress}" = "xyes" ); then
    ENABLE_PROGRESS=1
    AC_MSG_RESULT([yes])
else
    ENABLE_PROGRESS=0
    AC_MSG_RESULT([no])
fi;
AC_SUBST([ENABLE_PROGRESS])

ROOTDIR=$PWD
AC_SUBST([ROOTDIR])
AC_CONFIG_FILES([config.h Makefile gap.g])
AC_OUTPUT
