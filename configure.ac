### configure.ac                                        -*- Autoconf -*-
###
### Process this file with autoconf to produce a configure script.
###
### Evaporative Cooling (EC)
### Bill White (c) 2005-2012
### http://insilico.utulsa.edu/evaporative-cooling

AC_PREREQ(2.61)

AC_INIT(ec, 1.0, bill-whit@utulsa.edu, http://insilico.utulsa.edu/evaporative-cooling)

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR(src/library/EvaporativeCooling.cpp)
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE(ec, 1.0)

AC_SUBST(SRCDIR, src)
AC_SUBST(LIBDIR, library)

AC_CANONICAL_BUILD
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC
AC_F77_LIBRARY_LDFLAGS
AC_PROG_LN_S
AC_PROG_CXX
AC_PROG_LIBTOOL

# Use libtool for all builds
LT_INIT

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for header files.
AC_CHECK_HEADERS([math.h stdlib.h string.h sys/time.h])
# AC_CHECK_FUNCS([round pow sqrt floor memset], [], AC_MSG_ERROR([required function missing]))

# makeinfo
AC_CHECK_PROG([MAKEINFO], [makeinfo], [makeinfo], AC_MSG_ERROR([could not find required program: makeinfo]))

AC_MSG_NOTICE([--------------------------------------------------])

# add parameter(s)
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[Compile code for debugging.]), [debug_ec=yes])

# Set variables for basic compile debug/release
if test "$debug_ec" != 'yes' ; then
    CPPFLAGS="-O2 -g0 -Wall -DNDEBUG $CPPFLAGS"
    CXXFLAGS="-O2 -g0 -Wall -DNDEBUG $CXXFLAGS"
else
    CPPFLAGS="-O0 -gdwarf-2 -Wall $CPPFLAGS -D__DEBUG__"
    CXXFLAGS="-O0 -gdwarf-2 -Wall $CXXFLAGS -D__DEBUG__"
fi
CXXFLAGS="$CXXFLAGS -D__NOPLUGIN__ -I/usr/local/include"
CXXFLAGS="$CXXFLAGS -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen"
CPPFLAGS="$CPPFLAGS -D__NOPLUGIN__ -I/usr/local/include"
CPPFLAGS="$CPPFLAGS -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen"

AC_CHECK_LIB(pthread, pthread_create)
LDFLAGS="-lpthread $LDFLAGS"

# OpenMP
CPPFLAGS="$CPPFLAGS -fopenmp"
CXXFLAGS="$CXXFLAGS -fopenmp"
LDFLAGS="$LDFLAGS -fopenmp -lgomp"

# check for RandomJungle dependency
AC_CHECK_LIB([rjungle], [main])
AC_CHECK_LIB([lr], [mk_lr_options])
RJ_LIBS="-lrjungle -llr" 
AC_SUBST([RJ_LIBS])
LDFLAGS="$LDFLAGS -lrjungle -llr"

# boost.m4
BOOST_REQUIRE([1.34])
BOOST_PROGRAM_OPTIONS
BOOST_FIND_HEADER(boost/program_options.hpp, AC_MSG_ERROR([could not find required headers of Boost c++. Install boost c++ libraries.],))

# gsl.m4
AX_PATH_GSL(1.8,, AC_MSG_ERROR([could not find required version of GSL]))
GSL_LIBS='-lgsl -lgslcblas'
AC_SUBST([GSL_LIBS])
LDFLAGS="$LDFLAGS -lgsl -lgslcblas"

# zlib.m4
CHECK_ZLIB()

# libxml2.m4
AM_PATH_XML2(2.0.0,, AC_MSG_ERROR([could not find required version of XML2]))
#
# Check for XML delegate library.
#
AC_ARG_WITH([xml],
    [AC_HELP_STRING([--without-xml],
                    [disable XML support])],
    [with_xml=$withval],
    [with_xml=$have_x])

if test "$with_xml" != 'yes' ; then
    DISTCHECK_CONFIG_FLAGS="${DISTCHECK_CONFIG_FLAGS} --with-xml=$with_xml "
fi

have_xml='no'
XML_LIBS=''
if test "$with_xml" != 'no'; then
    AC_MSG_CHECKING([for XML])
    AC_MSG_RESULT([])
    PERSIST_LDFLAGS=$LDFLAGS
    PERSIST_CPPFLAGS=$CPPFLAGS
    xml2_config=''
    AC_CHECK_PROGS(xml2_config,xml2-config,)dnl
    if test -n "$xml2_config"; then
        # Debian installs libxml headers under /usr/include/libxml2/libxml with
        # the shared library installed under /usr/lib, whereas the package
        # installs itself under $prefix/libxml and $prefix/lib.
        xml2_prefix=`xml2-config --prefix`
        if test -d "${xml2_prefix}/include/libxml2"; then
            CPPFLAGS="$CPPFLAGS -I${xml2_prefix}/include/libxml2"
        fi
        if test "${xml2_prefix}" != '/usr'; then
            LDFLAGS="$LDFLAGS -L${xml2_prefix}/lib"
        fi
    fi
    failed=0
    passed=0
    AC_CHECK_HEADER(libxml/parser.h,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
    AC_CHECK_LIB(xml2,xmlParseExternalEntity,passed=`expr $passed + 1`,failed=`expr $failed + 1`,)
    AC_MSG_CHECKING([if XML package is complete])
    if test $passed -gt 0; then
        if test $failed -gt 0; then
            AC_MSG_RESULT([no -- some components failed test])
            have_xml='no (failed tests)'
            LDFLAGS="$PERSIST_LDFLAGS"
            CPPFLAGS="$PERSIST_CPPFLAGS"
        else
            XML_LIBS='-lxml2'
            LIBS="$XML_LIBS $LIBS"
            AC_DEFINE(XML_DELEGATE,1,Define if you have XML library)
            AC_MSG_RESULT([yes])
            have_xml='yes'
        fi
    else
        AC_MSG_RESULT([no])
    fi
fi

AM_CONDITIONAL(XML_DELEGATE,test "$have_xml" = 'yes')
AC_SUBST(XML_LIBS)

######################################################################
# DOXYGEN SUPPORT
######################################################################

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg])

######################################################################

#LIBS="$LIBS -lm -lgfortran"
#LIBDIR="$LIBDIR -L/usr/local/lib"

AC_MSG_NOTICE([--------------------------------------------------])

# Process makefiles
AC_CONFIG_FILES([Makefile src/library/Makefile src/Makefile examples/Makefile])
AC_OUTPUT
