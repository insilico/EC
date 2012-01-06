LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@
mycppflags=$(OPENMP_CFLAGS) $(BOOST_CPPFLAGS) -D__NOPLUGIN__

### SHARED
# the program to build
# bin_PROGRAMS = ec ec_static
bin_PROGRAMS = ec

# compiler flags
ec_CPPFLAGS = $(mycppflags)
ec_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
ec_LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS)

# prog's sources
ec_SOURCES = EvaporativeCoolingCLI.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefor relative path where it will be built!)
ec_LDADD = $(OPENMP_CFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) -L$(top_srcdir)/src/$(LIBDIR) -lec-1.0

# compiler flags
#ec_static_CPPFLAGS = $(mycppflags)

# prog's sources
#ec_static_SOURCES = EvaporativeCoolingCLI.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefor relative path where it will be built!)
#ec_static_LDADD = -L$(top_srcdir)/src/$(LIBDIR) -lec-1.0

#ec_static_LDFLAGS = -all-static $(OPENMP_CFLAGS)

# additional include pathes to support the library headers
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/$(LIBDIR)

# the subdirectories of src to go into
SUBDIRS = library

