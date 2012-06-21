AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

my_include_paths = -I$(top_srcdir)/src -I$(top_srcdir)/src/library

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS=$(my_include_paths) -x f77-cpp-input

my_ld_paths = $(LDFLAGS) -L$(top_srcdir)/src/$(LIBDIR)
my_ld_paths += -L/usr/local/lib
my_ld_flags = -lrjungle -llr -lec @FLIBS@ -lgsl -lgslcblas -lboost_program_options

# make ec and ecm (iMproved!)
bin_PROGRAMS = ec ecm ec_static

ec_SOURCES=EvaporativeCoolingCLI.cpp

# libtool libraries
ec_LDADD=$(my_ld_paths) $(my_ld_flags) 

ec_static_SOURCES=EvaporativeCoolingCLI.cpp
# libtool libraries
ec_static_LDFLAGS=-all-static -static-libgcc -static-libstdc++
ec_static_LDADD=$(my_ld_paths) $(my_ld_flags) 

ecm_SOURCES=EvaporativeCoolingCLINewAPI.cpp
# libtool libraries
ecm_LDADD=$(my_ld_paths) $(my_ld_flags) 

# additional include paths to support the library headers
INCLUDES = $(myincludes)

# the subdirectories of src to go into
SUBDIRS=library
