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
example1_LDADD=-lec -lrjungle -llr \
$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
example2_LDADD=-lec -lrjungle -llr \
$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
example3_LDADD=-lec -lrjungle -llr \
$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
example4_LDADD=-lec -lrjungle -llr \
$(OPENMP_CXXFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
