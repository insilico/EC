AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@

myincludes = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_srcdir)/src/library
mycppflags = $(OPENMP_CXXFLAGS) $(BOOST_CPPFLAGS) $(myincludes)

# shared library build
lib_LTLIBRARIES=libec.la
libec_la_SOURCES=\
Dataset.cpp ArffDataset.cpp PlinkDataset.cpp  PlinkBinaryDataset.cpp \
PlinkRawDataset.cpp DgeData.cpp BirdseedData.cpp DatasetInstance.cpp DistanceMetrics.cpp \
FilesystemUtils.cpp Statistics.cpp ChiSquared.cpp Insilico.cpp \
ReliefF.cpp RReliefF.cpp RandomJungle.cpp EvaporativeCooling.cpp

# the list of header files that belong to the library
libec_la_HEADERS= \
ArffDataset.h ChiSquared.h Dataset.h DatasetInstance.h \
Debugging.h DistanceMetrics.h EvaporativeCooling.h FilesystemUtils.h \
GSLRandomBase.h GSLRandomFlat.h Insilico.h PlinkBinaryDataset.h PlinkDataset.h \
PlinkRawDataset.h RReliefF.h RandomJungle.h ReliefF.h DgeData.h Statistics.h \
StringUtils.h best_n.h BirdseedData.h

libec_la_CPPFLAGS = $(mycppflags)
libec_la_LDFLAGS = -release 1.0.0

# where to install the headers on the system
libec_ladir=$(includedir)/ec

#library_includedir=$(includedir)
#library_include_HEADERS = ec/ArffDataset.h ec/ChiSquared.h ec/CleanSnpDataset.h \
ec/Dataset.h ec/DatasetInstance.h ec/Debugging.h ec/DistanceMetrics.h \
ec/EvaporativeCooling.h ec/FilesystemUtils.h ec/GSLRandomBase.h \
ec/GSLRandomFlat.h ec/Insilico.h ec/PlinkBinaryDataset.h ec/PlinkDataset.h \
ec/PlinkRawDataset.h ec/RReliefF.h ec/RandomJungle.h ec/ReliefF.h \
ec/Statistics.h ec/StringUtils.h ec/best_n.h ec/BirdseedData.h
