AUTOMAKE_OPTIONS = foreign
LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@

myincludes = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_srcdir)/src/library
mycppflags = $(CXXFLAGS) $(BOOST_CPPFLAGS) $(myincludes)
AM_CXXFLAGS = $(my_cpp_flags) -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
 
# shared library build
lib_LTLIBRARIES=libec.la
libec_la_SOURCES=\
Insilico.cpp DistanceMetrics.cpp Statistics.cpp \
Dataset.cpp ArffDataset.cpp PlinkDataset.cpp  PlinkBinaryDataset.cpp \
PlinkRawDataset.cpp DgeData.cpp BirdseedData.cpp DatasetInstance.cpp \
AttributeRanker.cpp ChiSquared.cpp RandomJungle.cpp Deseq.cpp \
ReliefF.cpp RReliefF.cpp SNReliefF.cpp ReliefFSeq.cpp \
EvaporativeCooling.cpp

# the list of header files that belong to the library
libec_la_HEADERS= \
ArffDataset.h ChiSquared.h Dataset.h DatasetInstance.h DistanceMetrics.h \
EvaporativeCooling.h GSLRandomBase.h GSLRandomFlat.h Insilico.h \
PlinkBinaryDataset.h PlinkDataset.h PlinkRawDataset.h AttributeRanker.h \
RReliefF.h SNReliefF.h ReliefFSeq.h RandomJungle.h ReliefF.h Deseq.h \
DgeData.h Statistics.h StringUtils.h BestN.h BirdseedData.h

libec_la_CPPFLAGS = $(mycppflags)
libec_la_LDFLAGS = $(LDFLAGS) -L/usr/local/lib -lxml2 -lz -lpthread -lm \
-release 1.0.0 -lgsl -lgslcblas -lgfortran -lgomp

# where to install the headers on the system
libec_ladir=$(includedir)/ec

#library_includedir=$(includedir)
#library_include_HEADERS = ec/ArffDataset.h ec/ChiSquared.h ec/CleanSnpDataset.h \
ec/Dataset.h ec/DatasetInstance.h ec/Debugging.h ec/DistanceMetrics.h \
ec/EvaporativeCooling.h ec/FilesystemUtils.h ec/GSLRandomBase.h \
ec/GSLRandomFlat.h ec/Insilico.h ec/PlinkBinaryDataset.h ec/PlinkDataset.h \
ec/PlinkRawDataset.h ec/RReliefF.h ec/RandomJungle.h ec/ReliefF.h \
ec/Statistics.h ec/StringUtils.h ec/best_n.h ec/BirdseedData.h
