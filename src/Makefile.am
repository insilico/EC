AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

my_include_paths = -I$(top_srcdir)/src -I$(top_srcdir)/src/library

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS = $(my_include_paths) -x f77-cpp-input

my_ld_flags = -lrjungle -llr -lec -lgsl -lgslcblas -lboost_program_options @FLIBS@
my_ld_paths = -L$(top_srcdir)/src/$(LIBDIR) $(LDFLAGS) 
# my_ld_paths += -L/usr/local/lib

# make ec and ecm (iMproved!)
bin_PROGRAMS = ec ecm ec_static

ec_SOURCES = EvaporativeCoolingCLI.cpp

# libtool libraries
ec_LDADD = $(my_ld_flags) $(my_ld_paths)

ec_static_SOURCES=EvaporativeCoolingCLI.cpp
# libtool libraries
ec_static_LDFLAGS = -all-static -static-libgcc -static-libstdc++
ec_static_LDADD = $(my_ld_flags) $(my_ld_paths)

ecm_SOURCES = EvaporativeCoolingCLINewAPI.cpp
# libtool libraries
ecm_LDADD = $(my_ld_flags) $(my_ld_paths)

# additional include paths to support the library headers
INCLUDES = $(myincludes)

# the subdirectories of src to go into
SUBDIRS = library
