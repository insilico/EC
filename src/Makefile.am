AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@

myincludes = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_srcdir)/src/$(LIBDIR)
mycppflags = $(OPENMP_CXXFLAGS) $(BOOST_CPPFLAGS) $(myincludes)

# make ec and ecm (iMproved!)
bin_PROGRAMS = ec ecm ec_static

ec_SOURCES=EvaporativeCoolingCLI.cpp
# libtool libraries
ec_LDADD=$(OPENMP_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) $(RJ_LIBS) -lec

ec_static_SOURCES=EvaporativeCoolingCLI.cpp
# libtool libraries
ec_static_LDFLAGS=-all-static -static-libgcc -static-libstdc++
ec_static_LDADD=$(OPENMP_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) $(RJ_LIBS) -lec

ecm_SOURCES=EvaporativeCoolingCLINewAPI.cpp
# libtool libraries
ecm_LDADD=$(OPENMP_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) $(RJ_LIBS) -lec

# additional include paths to support the library headers
INCLUDES = $(myincludes)

# the subdirectories of src to go into
SUBDIRS=library
