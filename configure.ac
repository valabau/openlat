#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(openlat, 0.1, valabau@gmail.com)
AC_CONFIG_SRCDIR([src/include/openlat/compat.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([src/include/openlat/config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(openlat, 0.1)


LT_INIT
LTDL_INIT


YACC=bison
# Checks for programs.
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_LANG(C++)
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_YACC
AM_PROG_LEX

# boost
AX_BOOST_BASE(1.41)
AX_BOOST_UNIT_TEST_FRAMEWORK
AM_CONDITIONAL([MAKE_TEST], [test x"$BOOST_UNIT_TEST_FRAMEWORK_LIB" != x])


# doxygen
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(foobar, doxygen.cfg, doc)


# optional configure options
AC_ARG_WITH(debug, AC_HELP_STRING(--with-debug, \
[enable debug]), \
[CPPFLAGS="-DDEBUG=1" \
CXXFLAGS="-g -O0"])

AC_ARG_WITH(gprof, AC_HELP_STRING(--with-gprof, \
[enable profiler]), \
[CPPFLAGS="-DDEBUG=1" \
CXXFLAGS+=" -pg"])



# Checks for libraries.
AM_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([fst/fstlib.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.



# Check for compilation options
# if error on this function then AX_CHECK_COMPILE_FLAG is not defined.
# please, install autoconf-archive (sudo apt-get install autoconf-archive) 
 AX_CHECK_COMPILER_FLAGS([-std=c++11], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++11"], 
    [AX_CHECK_COMPILER_FLAGS([-std=c++0x], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++0x"])])
#AX_CHECK_COMPILER_FLAGS([-std=c++98], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++98"],)


#COMPILANCE_CXXFLAGS="${CPPCXXFLAGS} -D_POSIX_C_SOURCE=200809L"
AX_CHECK_COMPILER_FLAGS([-pedantic], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -pedantic"])
#AX_CHECK_COMPILER_FLAGS([-fpermissive], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -fpermissive"])
AX_CHECK_COMPILER_FLAGS([-Wall], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wall"])
AX_CHECK_COMPILER_FLAGS([-Wextra], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wextra"])
AX_CHECK_COMPILER_FLAGS([-Wshadow], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wshadow"])
AX_CHECK_COMPILER_FLAGS([-Winline], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Winline"])
AX_CHECK_COMPILER_FLAGS([-Wunused], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused"])
AX_CHECK_COMPILER_FLAGS([-Wunused-variable], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-variable"])
AX_CHECK_COMPILER_FLAGS([-Wunused-parameter], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-parameter"])
AX_CHECK_COMPILER_FLAGS([-Wunused-function], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-function"])

AM_CONDITIONAL(DEBUG, [test x"$debug" = x"true"])
AS_IF(test x"$DEBUG" = x"true", WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -DDEBUG -DLT_DEBUG_LOADERS"], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -DNDEBUG"])
#AS_IF(test x"$DEBUG" = x"true", OPTIMIZATION_CXXFLAGS="${OPTIMIZATION_CXXFLAGS} -DDEBUG -g3 -O0"], [OPTIMIZATION_CXXFLAGS="${OPTIMIZATION_CXXFLAGS} -O2 -ffast-math"])



AC_SUBST([WARNING_CXXFLAGS], ["${WARNING_CXXFLAGS}"])
AC_SUBST([COMPILANCE_CXXFLAGS])
AC_SUBST([OPTIMIZATION_CXXFLAGS])


# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/bin/Makefile])
AC_CONFIG_FILES([src/test/Makefile])


AC_OUTPUT
