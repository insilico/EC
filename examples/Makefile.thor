# -----------------------------------------------------------------------------
# Makefile.thor - Bill White - 1/10/11
#
# ReliefF Projects for the McKinney insilico Lab
# Builds standard binaries and shared libraries on Linux (thor).
# -----------------------------------------------------------------------------

PROJECT_NAME = insilico

# compile and link settings ----------------------------------------------------
C++C = g++

# Compile and link flags
C++C_COMPILE_FLAGS =-Wall -ansi -g -fopenmp -Wno-long-long -I. \
	-I/usr/include/libxml2 -I/usr/local/include -DHAVE__BOOL -DNDEBUG \
	-D__NOPLUGIN__ -fPIC -D_FILE_OFFSET_BITS=64
C++C_LINK_FLAGS = -lm -lgsl -lgslcblas -lboost_program_options \
	-lboost_filesystem  -lboost_system \
	-lgomp -lxml2 -llapack -lblas -L/usr/local/lib -lrjungle -lz -llr
C++C_LIB_LINK_FLAGS = -shared

RELIEFF_LIB_NAME = librelieff.so
EC_LIB_NAME = libec.so

include Makefile.common

