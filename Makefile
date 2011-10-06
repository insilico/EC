# -----------------------------------------------------------------------------
# Makefile - Bill White - 7/12/11
#
# Evaporative Cooling Projects for the McKinney insilico Lab
#
# -----------------------------------------------------------------------------

PROJECT_NAME = ec

# --------------------------------------------------- compile and link settings
C++C = g++

# Compile and link flags
C++C_COMPILE_FLAGS =-Wall -ansi -O2 -g -DNDEBUG -fopenmp \
	-I../afni_src -Wno-long-long -I../cpprelieff -I. \
	-I/usr/local/include/rjungle -I/usr/local/include -DHAVE__BOOL \
	-D__NOPLUGIN__ -I/usr/include/libxml2
C++C_LINK_FLAGS = -lm -lgsl -lgslcblas -lboost_program_options -lgomp -lmri \
  -L$(HOME)/abin -lf2c -lxml2 -L/usr/local/lib -lrjungle -lz -llr

# -------------------------------------------------------------------------- EC
EC_PROG_NAME = ec
EC_OBJS = EvaporativeCoolingCLI.o EvaporativeCooling.o ../cpprelieff/Dataset.o \
	../cpprelieff/DatasetInstance.o ../cpprelieff/Statistics.o \
	../cpprelieff/ArffDataset.o ../cpprelieff/PlinkDataset.o \
	../cpprelieff/PlinkRawDataset.o ../cpprelieff/PlinkBinaryDataset.o \
	../cpprelieff/DistanceMetrics.o ../cpprelieff/ChiSquared.o \
	../cpprelieff/FilesystemUtils.o ../cpprelieff/ReliefF.o \
	../cpprelieff/RReliefF.o ../cpprelieff/CleanSnpDataset.o

# ----------------------------------------------------------------- BUILD RULES

all: $(EC_PROG_NAME)
	@echo "$(EC_PROG_NAME) is now up to date."

$(EC_PROG_NAME): $(EC_OBJS)
	$(C++C) -o $(EC_PROG_NAME) $(EC_OBJS) $(C++C_LINK_FLAGS)

%.o: %.C
	$(C++C) -c $(C++C_COMPILE_FLAGS) $<

%.o: %.cpp
	$(C++C) -c $(C++C_COMPILE_FLAGS) $<

%.o: %.c
	$(C++C) -c $(C++C_COMPILE_FLAGS) $<

clean:
	rm -f *.o core *~

reallyclean:
	rm -f *.o core $(EC_PROG_NAME) \
	
tarball:
	rm -f *.o core $(EC_PROG_NAME)
	tar cvfpj ~/archive/$(PROJECT_NAME)_`date '+%Y%m%d'`.tar.bz2 *
