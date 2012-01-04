#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Simple_test.txt
CND_DISTDIR=dist

# Include project Makefile
include Makefile-cppec.mk

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES=


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	cd . && $(MAKE) 

# Subprojects
.build-subprojects:
	cd ../cpprelieff_fucku && ${MAKE}  -f cpprelieff_fucku-Makefile.mk CONF=Default

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	cd . && $(MAKE) clean

# Subprojects
.clean-subprojects:
	cd ../cpprelieff_fucku && ${MAKE}  -f cpprelieff_fucku-Makefile.mk CONF=Default clean
