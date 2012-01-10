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

bin_PROGRAMS = example1 example2 example3 example4
example1_SOURCES=example1.cpp
example2_SOURCES=example2.cpp
example3_SOURCES=example3.cpp
example4_SOURCES=example4.cpp

# libtool libraries
example1_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec-1.0 \
-lrjungle-1.2.362 -llr-1.2.362
example2_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec-1.0 \
-lrjungle-1.2.362 -llr-1.2.362
example3_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec-1.0 \
-lrjungle-1.2.362 -llr-1.2.362
example4_LDADD=$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
-L$(top_srcdir)/src/$(LIBDIR) -lec-1.0 \
-lrjungle-1.2.362 -llr-1.2.362

# additional include paths to support the library headers
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/$(LIBDIR)
