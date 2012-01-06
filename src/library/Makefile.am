LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@
myincludes = -I$(top_srcdir)/src -I$(top_srcdir)
mycppflags=$(OPENMP_CFLAGS) $(myincludes)
AM_FFLAGS=$(myincludes) -x f77-cpp-input

myheaders= \
ArffDataset.h ChiSquared.h CleanSnpDataset.h Dataset.h DatasetInstance.h \
Debugging.h DistanceMetrics.h EvaporativeCooling.h FilesystemUtils.h \
GSLRandomBase.h GSLRandomFlat.h Insilico.h PlinkBinaryDataset.h PlinkDataset.h \
PlinkRawDataset.h RReliefF.h RandomJungle.h ReliefF.h

mysources= \
ArffDataset.cpp ChiSquared.cpp CleanSnpDataset.cpp Dataset.cpp \
DatasetInstance.cpp DistanceMetrics.cpp EvaporativeCooling.cpp \
FilesystemUtils.cpp Insilico.cpp PlinkBinaryDataset.cpp \
PlinkDataset.cpp PlinkRawDataset.cpp RReliefF.cpp RandomJungle.cpp \
ReliefF.cpp Statistics.cpp

# shared library build
lib_LTLIBRARIES=libec-1.0.la
libec_1_0_la_SOURCES = $(mysources)

# compiler flags
libec_1_0_la_CPPFLAGS=$(mycppflags)

# the list of header files that belong to the library
libec_1_0_la_HEADERS=$(myheaders)

# where to install the headers on the system
libec_1_0_ladir=$(includedir)/ec

#library_includedir=$(includedir)
#library_include_HEADERS = ec/ArffDataset.h ec/ChiSquared.h ec/CleanSnpDataset.h \
ec/Dataset.h ec/DatasetInstance.h ec/Debugging.h ec/DistanceMetrics.h \
ec/EvaporativeCooling.h ec/FilesystemUtils.h ec/GSLRandomBase.h \
ec/GSLRandomFlat.h ec/Insilico.h ec/PlinkBinaryDataset.h ec/PlinkDataset.h \
ec/PlinkRawDataset.h ec/RReliefF.h ec/RandomJungle.h ec/ReliefF.h \
ec/Statistics.h ec/StringUtils.h ec/best_n.h
