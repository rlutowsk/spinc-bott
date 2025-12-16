AC_ARG_VAR( WARNINGS, keeps compiler warnings options )

AC_PREREQ([2.71])

AC_INIT([Spinc structures on real Bott manifolds],[0.8],[rafal.lutowski@ug.edu.pl])

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AS_IF([test "x$ac_cv_env_CFLAGS_set" = x],
      [
        CFLAGS="-O3 -march=native"
      ])

AC_PROG_CXX
AC_PROG_EGREP

# Ensure m4 macros for pkg-config are included
m4_ifndef([PKG_PROG_PKG_CONFIG], [m4_fatal([you need pkg-config])])
PKG_PROG_PKG_CONFIG
if test "x$PKG_CONFIG" = "xfalse"; then
  AC_MSG_ERROR([The required tool pkg-config was not found.])
fi

# Opcja do wskazania ogólnego, niestandardowego prefixu dla bibliotek
AC_MSG_CHECKING([for a local library prefix])
AC_ARG_WITH(local-prefix,
  [AS_HELP_STRING([--with-local-prefix=PREFIX],
                  [use libraries from a local prefix (sets PKG_CONFIG_PATH, CPPFLAGS, LDFLAGS)])],
  [local_prefix="$withval"],
  [local_prefix=""])

if test -n "$local_prefix"; then
  AC_MSG_RESULT([$local_prefix])
  # 1. Configure PKG_CONFIG_PATH
  pkgconfig_dir="${local_prefix}/lib/pkgconfig"
  if test -d "$pkgconfig_dir"; then
    if test -n "$PKG_CONFIG_PATH"; then
      PKG_CONFIG_PATH="${pkgconfig_dir}:${PKG_CONFIG_PATH}"
    else
      PKG_CONFIG_PATH="${pkgconfig_dir}"
    fi
    # Wyeksportuj zmienną, aby była widoczna dla poleceń uruchamianych przez configure
    export PKG_CONFIG_PATH
  fi

  # 2. Configure CPPFLAGS
  local_incdir="${local_prefix}/include"
  CPPFLAGS="${CPPFLAGS}${CPPFLAGS:+ }-I${local_incdir}"

  # 3. Configure LDFLAGS
  local_libdir="${local_prefix}/lib"
  LDFLAGS="${LDFLAGS}${LDFLAGS:+ }-L${local_libdir} -Wl,-rpath,${local_libdir}"
  RPATH="-Wl,-rpath,${local_libdir}"
  AC_SUBST(RPATH)
else
  AC_MSG_RESULT([no])
fi

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

PKG_CHECK_MODULES([NAUTY], [libnauty >= 2.9.1],
                  [AC_SUBST(NAUTY_CFLAGS)]
                  [AC_SUBST(NAUTY_LIBS)],
                  [AC_MSG_ERROR([nauty 2.9.1 or newer is required.])])

PKG_CHECK_MODULES([XXHASH], [libxxhash >= 0.8.1],
                  [AC_SUBST(XXHASH_CFLAGS)]
                  [AC_SUBST(XXHASH_LIBS)],
                  [AC_MSG_ERROR([xxhash library required])])

AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
# AC_DEFINE([ENABLE_OPENMP],[1],[Compile operators with OpenMP support])
# CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
# Oddly, linker needs same OpenMP flag(s) as compiler
# 20170909 Linux clang on travis does not resolve __kmpc* symbols, apparently needs -lgomp or -lomp
# LDFLAGS="${LDFLAGS}"

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T


# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([bzero clock_gettime gettimeofday])

AC_CHECK_DECLS([CLOCK_MONOTONIC_RAW],
               [AC_DEFINE([HAVE_CLOCK_MONOTONIC_RAW], [1], [Define to 1 if CLOCK_MONOTONIC_RAW is available])],
               [],
               [#include <time.h>]) # header to check

# Set default options (final)
BUILD_MODE="final"
BUILD_MODE_CFLAGS="${CFLAGS}" # Variable to store mode-specific CFLAGS
BUILD_MODE_LDFLAGS="" # Variable to store mode-specific LDFLAGS

# -----------------------------------------------------------
# User option definitions
# -----------------------------------------------------------
AC_ARG_ENABLE([profile],
  [AS_HELP_STRING([--enable-profile], [Enable profiling build mode (implies -pg)])],
  [test "x$enable_profile" = "xyes" && BUILD_MODE="profile"])

AC_ARG_ENABLE([sanitize],
  [AS_HELP_STRING([--enable-sanitize], [Enable sanitizing build mode])],
  [test "x$enable_sanitize" = "xyes" && BUILD_MODE="sanitize"])

AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug], [Enable debug build mode])],
  [test "x$enable_debug" = "xyes" && BUILD_MODE="debug"])

# -----------------------------------------------------------
# Logic for mutual exclusivity and settings
# -----------------------------------------------------------

AC_MSG_NOTICE([build mode set to: $BUILD_MODE])

AS_IF(
  [test "x$BUILD_MODE" = "xfinal"],
    [], # No special flags
  [test "x$BUILD_MODE" = "xprofile"],
    [AC_DEFINE([PROFILE], [1], [Define to 1 for profiling mode.])
     BUILD_MODE_CFLAGS="-pg"
     BUILD_MODE_LDFLAGS="-pg"],
  [test "x$BUILD_MODE" = "xsanitize"],
    [AC_DEFINE([SANITIZE], [1], [Define to 1 for sanitizing mode.])
     BUILD_MODE_CFLAGS="-g -O2 -fsanitize=address -fno-omit-frame-pointer"
     BUILD_MODE_LDFLAGS="-fsanitize=address,undefined"],
  [test "x$BUILD_MODE" = "xdebug"],
    [AC_DEFINE([DEBUG], [1], [Define to 1 for debug mode.])
     BUILD_MODE_CFLAGS="-g -O0"],
  [AC_MSG_ERROR([Internal error: Invalid build mode: $BUILD_MODE])]
)

AC_SUBST([BUILD_MODE_CFLAGS])
AC_SUBST([BUILD_MODE_LDFLAGS])

AC_ARG_ENABLE([assert],
  [AS_HELP_STRING([--disable-assert], [Disable assertion build mode (define NDEBUG)])],
  [enable_assert="$enableval"],
  [enable_assert=yes])


AS_IF(
  [test "x$enable_assert" = "xno"], 
    [AC_DEFINE([NDEBUG], [1], [Define to disable asserts and enable optimizations.])]
)


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

ROOTDIR=$PWD
AC_SUBST([ROOTDIR])
AC_CONFIG_FILES([Makefile gap.g])

AC_OUTPUT
