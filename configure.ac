#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(openlat, 0.1, valabau@gmail.com)
AC_CONFIG_SRCDIR([src/include/openlat/compat.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([src/include/openlat/config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE


LT_INIT
LTDL_INIT


# Checks for programs.
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_CXX
AC_LANG(C++)
AC_PROG_INSTALL
AC_PROG_MAKE_SET
#AM_PROG_LEX
AC_CHECK_PROGS(LEX,flex lex, no)
if test "x$LEX" = "xno" ; then
    AC_MSG_ERROR([lex is required])
fi
#AC_PROG_YACC
AC_CHECK_PROGS(YACC,'bison -y' byacc yacc,no)
if test "x$YACC" = "xno" ; then
    AC_MSG_ERROR([bison is required])
fi
BISON_VERSION=`bison --version | awk 'NR==1{print $NF}'`
AC_SUBST([BISON_VERSION])

# boost
AX_BOOST_BASE(1.41)
AX_BOOST_UNIT_TEST_FRAMEWORK
AM_CONDITIONAL([MAKE_TEST], [test x"$BOOST_UNIT_TEST_FRAMEWORK_LIB" != x])
AX_BOOST_IOSTREAMS # BOOST_IOSTREAMS_LIB


# doxygen
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(openlat, doxygen.cfg, doc)


# optional configure options
AC_ARG_WITH(debug, AC_HELP_STRING(--with-debug, \
[enable debug]), \
[ DEBUG=true \
CXXFLAGS="-g3 -O0"])

AC_ARG_WITH(gprof, AC_HELP_STRING(--with-gprof, \
[enable profiler]), \
[ PROFILE=true \
CPPFLAGS="-DDEBUG=1" \
CXXFLAGS+=" -pg"])

AC_ARG_WITH(thot, AC_HELP_STRING(--with-thot, [enable thot]), [ HAVE_THOT=true ])
AM_CONDITIONAL([HAVE_THOT], [test "x$HAVE_THOT" != x])
if test "x$HAVE_THOT" != "x" ; then
	AC_DEFINE([HAVE_THOT], [1], [Define if Thot is found ])
fi

#CPPFLAGS+=" -DGIT_SHA1=\"\\\"\`git describe --abbrev=64 --always --dirty\`\\\"\""

# Checks for libraries.
AM_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([fst/fstlib.h])

AC_CHECK_HEADERS([execinfo.h])
AC_CHECK_TOOL([ADDR2LINE], [addr2line], [:])

if test "x$ADDR2LINE" != "x" ; then
	AC_DEFINE_UNQUOTED(ADDR2LINEBIN, ["$ADDR2LINE"], [Define if add2line is found ])
fi
# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.

# check for Python
AM_PATH_PYTHON
AX_PYTHON_DEVEL([>= 2.5])
AM_CONDITIONAL([HAVE_PYTHON_DEVEL], [test x"$PYTHON_EXTRA_LIBS" != x])

if [ test "x$PYTHON_EXTRA_LIBS" != x ]; then
  AC_DEFINE([HAVE_PYTHON], [1], [Whether Python devel is available])
  
  # Swig python bindings
  AX_PKG_SWIG(1.3.40, [], [ AC_MSG_WARN([SWIG is required to build python bindings]) ])
  AX_SWIG_ENABLE_CXX
  AX_SWIG_PYTHON

  AM_CONDITIONAL([HAVE_SWIG], [test "x$SWIG" != x])
  
  AM_CONDITIONAL([HAVE_BOOST_PYTHON], test "x$boost_python_cpplib_avail_cv_" = xtrue)
fi



# Check for compilation options
# if error on this function then AX_CHECK_COMPILE_FLAG is not defined.
# please, install autoconf-archive (sudo apt-get install autoconf-archive) 
# AX_CHECK_COMPILE_FLAG([-std=c++11], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++11"], 
#    [AX_CHECK_COMPILE_FLAG([-std=c++0x], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++0x"])])
#AX_CHECK_COMPILE_FLAG([-std=c++98], [COMPILANCE_CXXFLAGS="${COMPILANCE_CXXFLAGS} -std=c++98"],)


AX_CXX_COMPILE_STDCXX_11([noext])
#COMPILANCE_CXXFLAGS="${CPPCXXFLAGS} -D_POSIX_C_SOURCE=200809L"
AX_CHECK_COMPILE_FLAG([-pedantic], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -pedantic"])
#AX_CHECK_COMPILE_FLAG([-fpermissive], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -fpermissive"])
AX_CHECK_COMPILE_FLAG([-Wall], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wall"])
AX_CHECK_COMPILE_FLAG([-Wextra], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wextra"])
AX_CHECK_COMPILE_FLAG([-Wshadow], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wshadow"])
# enable this on profiling since it will show which inlines where not made
#AX_CHECK_COMPILE_FLAG([-Winline], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Winline"])
#AX_CHECK_COMPILE_FLAG([-Wunused], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused"])
#AX_CHECK_COMPILE_FLAG([-Wunused-variable], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-variable"])
#AX_CHECK_COMPILE_FLAG([-Wunused-parameter], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-parameter"])
#AX_CHECK_COMPILE_FLAG([-Wunused-function], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -Wunused-function"])

AS_IF(test x"$DEBUG" = x"true" -o "$PROFILE" = x"true" , [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -DDEBUG -DLT_DEBUG_LOADERS"], [WARNING_CXXFLAGS="${WARNING_CXXFLAGS} -DNDEBUG"])
#AS_IF(test x"$DEBUG" = x"true", [OPTIMIZATION_CXXFLAGS="${OPTIMIZATION_CXXFLAGS} -DDEBUG -g3 -O0"], [OPTIMIZATION_CXXFLAGS="${OPTIMIZATION_CXXFLAGS} -O2 -ffast-math"])



AC_SUBST([WARNING_CXXFLAGS], ["${WARNING_CXXFLAGS}"])
AC_SUBST([COMPILANCE_CXXFLAGS])
AC_SUBST([OPTIMIZATION_CXXFLAGS])


# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/bin/Makefile])
AC_CONFIG_FILES([src/lib/Makefile])
AC_CONFIG_FILES([src/test/Makefile])


AC_OUTPUT
