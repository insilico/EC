AUTOMAKE_OPTIONS = foreign
AM_CXXFLAGS = -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen @OS_CXX@
LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@

# Random Jungle flags
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
# project flags
AM_CXXFLAGS += -I$(top_srcdir)/src -I$(top_srcdir)
# library flags
AM_CXXFLAGS += $(OPENMP_CXXFLAGS) $(BOOST_CPPFLAGS)

# shared library build
lib_LTLIBRARIES=libec.la
libec_la_SOURCES=\
ArffDataset.cpp ChiSquared.cpp CleanSnpDataset.cpp Dataset.cpp \
DatasetInstance.cpp DistanceMetrics.cpp EvaporativeCooling.cpp \
FilesystemUtils.cpp Insilico.cpp PlinkBinaryDataset.cpp \
PlinkDataset.cpp PlinkRawDataset.cpp RReliefF.cpp RandomJungle.cpp \
ReliefF.cpp Statistics.cpp

# the list of header files that belong to the library
libec_la_HEADERS= \
ArffDataset.h ChiSquared.h CleanSnpDataset.h Dataset.h DatasetInstance.h \
Debugging.h DistanceMetrics.h EvaporativeCooling.h FilesystemUtils.h \
GSLRandomBase.h GSLRandomFlat.h Insilico.h PlinkBinaryDataset.h PlinkDataset.h \
PlinkRawDataset.h RReliefF.h RandomJungle.h ReliefF.h

libec_la_LDFLAGS = -release 1.0.0

# where to install the headers on the system
libec_ladir=$(includedir)/ec

#library_includedir=$(includedir)
#library_include_HEADERS = ec/ArffDataset.h ec/ChiSquared.h ec/CleanSnpDataset.h \
ec/Dataset.h ec/DatasetInstance.h ec/Debugging.h ec/DistanceMetrics.h \
ec/EvaporativeCooling.h ec/FilesystemUtils.h ec/GSLRandomBase.h \
ec/GSLRandomFlat.h ec/Insilico.h ec/PlinkBinaryDataset.h ec/PlinkDataset.h \
ec/PlinkRawDataset.h ec/RReliefF.h ec/RandomJungle.h ec/ReliefF.h \
ec/Statistics.h ec/StringUtils.h ec/best_n.h
