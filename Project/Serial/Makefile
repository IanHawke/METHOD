# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
CC = nvcc

# Module directory
MODULE_DIR = ./Src

# Rootfinder directory
RTFIND_DIR = ./CminpackLibrary

# Include directory
INC_DIR = ./Include

# Cminpack include directory
RTFIND_INC_DIR = ./CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./CminpackLibrary/Src

# Test directory location
TEST_DIR = ../Tests

# Doxygen documentation directory
DOC_DIR = ../Doxumentation

# C++ compiler flags
CXXFLAGS = -Xcompiler -fopenmp -Xcompiler -Wall

# NVIDIA compiler flags
NVFLAGS = -std=c++11 -Wno-deprecated-gpu-targets

# Sources
SRCS = main.cu \
			 simulation.cu \
			 srmhd.cu \
			 twoFluidEMHD.cu \
			 initFunc.cu \
			 simData.cu \
			 weno.cu \
			 RK2.cu \
			 rkSplit.cu \
			 backwardsRK.cu \
			 backRKArgs.cu \
			 SSP2.cu \
			 IMEX2Args.cu \
			 SSP3.cu \
			 IMEX3Args.cu \
			 boundaryConds.cu \
			 saveData.cu \
			 fluxVectorSplitting.cu \
			 srrmhd.cu


# Headers
HDRS = ${SRCS:.cu=.h} cudaErrorCheck.h

# Objects
OBJS = ${SRCS:.cu=.o}

# Rootfinder objects
RTFIND_OBJS = $(RTFIND_SRC_DIR)/dogleg.o \
							$(RTFIND_SRC_DIR)/dpmpar.o \
							$(RTFIND_SRC_DIR)/enorm.o \
							$(RTFIND_SRC_DIR)/fdjac1.o \
							$(RTFIND_SRC_DIR)/hybrd.o \
							$(RTFIND_SRC_DIR)/hybrd1.o \
							$(RTFIND_SRC_DIR)/qform.o \
							$(RTFIND_SRC_DIR)/qrfac.o \
							$(RTFIND_SRC_DIR)/r1mpyq.o \
							$(RTFIND_SRC_DIR)/r1updt.o

# Programmes
PROGS = ${SRCS:.cu=}

# Main programme executable
EXEC = main

# Rootfinder executables
RTFIND = buildRootfinder

# Doxygen documentation
DOXY = doxumentation

# Build everything
build : $(DOXY) $(RTFIND) $(EXEC)

# Make and run all tests on simulation scripts and run main
all : $(RTFIND) $(EXEC) $(DOXY)
	cd $(TEST_DIR) && $(MAKE) test
	@$(MAKE) run

# Build and run main (does not run any tests if called explicitly)
run : $(RTFIND) $(EXEC)
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ''
	@echo ''
	@echo '######################################################'
	@echo '#              Executing main programme              #'
	@echo '######################################################'
	@./$(EXEC)

clean :
	rm	-f	$(EXEC)	$(OBJS)	*.gch


#################
# Build objects #
#################

simData.o : $(MODULE_DIR)/simData.cu $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cu $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

twoFluidEMHD.o : $(MODULE_DIR)/twoFluidEMHD.cu $(INC_DIR)/twoFluidEMHD.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

srmhd.o : $(MODULE_DIR)/srmhd.cu $(INC_DIR)/srmhd.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

srrmhd.o : $(MODULE_DIR)/srrmhd.cu $(INC_DIR)/srrmhd.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cu	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

main.o : $(MODULE_DIR)/main.cu $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

weno.o : $(MODULE_DIR)/weno.cu $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

RK2.o : $(MODULE_DIR)/RK2.cu $(INC_DIR)/RK2.h
		$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

rkSplit.o : $(MODULE_DIR)/rkSplit.cu $(INC_DIR)/rkSplit.h $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

backwardsRK.o : $(MODULE_DIR)/backwardsRK.cu $(INC_DIR)/backwardsRK.h $(INC_DIR)/backRKArgs.h $(INC_DIR)/rkSplit.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

backRKArgs.o : $(MODULE_DIR)/backRKArgs.cu $(INC_DIR)/backRKArgs.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP2.o : $(MODULE_DIR)/SSP2.cu $(INC_DIR)/SSP2.h $(INC_DIR)/IMEX2Args.h $(INC_DIR)/timeInt.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

IMEX2Args.o : $(MODULE_DIR)/IMEX2Args.cu $(INC_DIR)/IMEX2Args.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP3.o : $(MODULE_DIR)/SSP3.cu $(INC_DIR)/SSP3.h $(INC_DIR)/IMEX3Args.h $(INC_DIR)/timeInt.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

IMEX3Args.o : $(MODULE_DIR)/IMEX3Args.cu $(INC_DIR)/IMEX3Args.h $(INC_DIR)/IMEX2Args.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cu $(INC_DIR)/boundaryConds.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

saveData.o : $(MODULE_DIR)/saveData.cu $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cu $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)


# Executable
main : $(OBJS) $(RTFIND_OBJS)
	@$(CC)	$^	-o $@ $(CXXFLAGS) $(NVFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects

doxumentation:
	@cd $(DOC_DIR) && doxygen Doxyfile

# Tidy up and re-do
rebuild : clean all