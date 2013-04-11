AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

my_include_paths = -I$(top_srcdir)/src -I$(top_srcdir)/src/library

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS = $(my_include_paths)

my_ld_flags = $(LDFLAGS) -lec -lgsl -lgslcblas $(BOOST_PROGRAM_OPTIONS_LIBS)
my_ld_paths = -L$(top_srcdir)/src/$(LIBDIR)  

# make ec and ecm (iMproved!)
bin_PROGRAMS = ec ecm

ec_SOURCES = EvaporativeCoolingCLI.cpp

# libtool libraries
ec_LDADD = $(my_ld_flags) $(my_ld_paths)

ecm_SOURCES = EvaporativeCoolingCLINewAPI.cpp
# libtool libraries
ecm_LDADD = $(my_ld_flags) $(my_ld_paths)

# additional include paths to support the library headers
INCLUDES = $(myincludes)

# the subdirectories of src to go into
SUBDIRS = library
