AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

my_include_paths = -I$(top_srcdir)/src -I$(top_srcdir)/src/library

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS = $(my_include_paths) -x f77-cpp-input

my_ld_flags = $(LDFLAGS) -lec -lgsl -lgslcblas $(BOOST_PROGRAM_OPTIONS_LIBS) @FLIBS@
my_ld_paths = -L$(top_srcdir)/src/$(LIBDIR)  

noinst_PROGRAMS = example1 example2 example3 example4
example1_SOURCES=example1.cpp
example2_SOURCES=example2.cpp
example3_SOURCES=example3.cpp
example4_SOURCES=example4.cpp

# libtool libraries
example1_LDADD=$(my_ld_flags) $(my_ld_paths)
example2_LDADD=$(my_ld_flags) $(my_ld_paths)
example3_LDADD=$(my_ld_flags) $(my_ld_paths)
example4_LDADD=$(my_ld_flags) $(my_ld_paths)
