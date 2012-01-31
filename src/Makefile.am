AUTOMAKE_OPTIONS = foreign
AM_CXXFLAGS = -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen @OS_CXX@
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

# Random Jungle flags
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
# project flags
AM_CXXFLAGS += -I$(top_srcdir)/src/library -I$(top_srcdir)/src -I$(top_srcdir)
# library flags
AM_CXXFLAGS += $(OPENMP_CXXFLAGS) $(BOOST_CPPFLAGS)

# make ec and ecm (iMproved!)
bin_PROGRAMS = ec ecm

ec_SOURCES=EvaporativeCoolingCLI.cpp
# libtool libraries
ec_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec -lrjungle -llr

ecm_SOURCES=EvaporativeCoolingCLINewAPI.cpp
# libtool libraries
ecm_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec -lrjungle -llr

# additional include paths to support the library headers
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/$(LIBDIR)

# the subdirectories of src to go into
SUBDIRS=library
