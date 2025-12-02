
COMP := $(if $(COMP),$(COMP),1)
$(info COMP = $(COMP))

PERIODIC := $(if $(PERIODIC),$(PERIODIC),0)
$(info PERIODIC = $(PERIODIC))

HDF5 := $(if $(HDF5),$(HDF5),0)
$(info HDF5 = $(HDF5))

#### This is an example for a local machine ####
ifeq ($(COMP),1)          
 FC=gfortran
 FFLAGS=-O2 -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -funroll-loops
 LIBS=
 INC=
endif

#### This is an example for a local machine DEBUGGING ####
ifeq ($(COMP),2)          
 FC=gfortran
 FFLAGS=-O1 -g -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
 LIBS=
 INC=
endif
##########################################################################

# DIRECTORIES
SRCDIR     := src
BINDIR     := bin

# -J is the directory where the .mod files are stored (BIN)
FFLAGS  += -J $(BINDIR) 

# Including HDF5 library if needed
ifeq ($(HDF5),1)
 INC += -I/usr/include/hdf5/serial/
 LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_serial_fortran -lhdf5hl_fortran
endif

# cosmokdtree precompiled variables
LONGINT = 1
DOUBLEPRECISION = 0
DIMEN = 3

FLAGS_ALL = -cpp
FLAGS_ALL += $(FFLAGS)
FLAGS_ALL += -Dperiodic=$(PERIODIC) -Dlongint=$(LONGINT) -Ddoubleprecision=$(DOUBLEPRECISION) -Ddimen=$(DIMEN)
FLAGS_ALL += -Duse_hdf5=$(HDF5)
 
# EXECUTABLE
EXEC=avism.x

# OBJECTS
OBJ=commondata.o kdtree.o particles.o voidfinding.o avism.o

# Just to show the final executable when typing 'make'
all: $(EXEC)

# make sure that ./bin/ folder exists before any .o goes there
$(addprefix $(BINDIR)/, $(OBJ)): | $(BINDIR)

$(BINDIR):
	@mkdir -p $@

# COMPILATION
$(EXEC): $(addprefix $(BINDIR)/, $(OBJ))
	$(FC) $(FLAGS_ALL) $(addprefix $(BINDIR)/, $(OBJ)) -o $(EXEC) $(LIBS)

# Rule for commondata.o
$(BINDIR)/commondata.o: $(SRCDIR)/commondata.f90
	$(FC) $(FLAGS_ALL) $(INC) -c -o $(BINDIR)/commondata.o $(SRCDIR)/commondata.f90

# Rule for kdtree.o
$(BINDIR)/kdtree.o: $(SRCDIR)/kdtree.f90
	$(FC) $(FLAGS_ALL) $(INC) -c -o $(BINDIR)/kdtree.o $(SRCDIR)/kdtree.f90

# Rule for particles.o
$(BINDIR)/particles.o: $(SRCDIR)/particles.f90
	$(FC) $(FLAGS_ALL) $(INC) -c -o $(BINDIR)/particles.o $(SRCDIR)/particles.f90

# Rule for voidfinding.o
$(BINDIR)/voidfinding.o: $(SRCDIR)/voidfinding.f90
	$(FC) $(FLAGS_ALL) $(INC) -c -o $(BINDIR)/voidfinding.o $(SRCDIR)/voidfinding.f90


# Rule for avism.o 
$(BINDIR)/avism.o: $(SRCDIR)/avism.f90 $(BINDIR)/commondata.o $(BINDIR)/kdtree.o $(BINDIR)/particles.o $(BINDIR)/voidfinding.o
	$(FC) $(FLAGS_ALL) $(INC) -c -o $(BINDIR)/avism.o $(SRCDIR)/avism.f90

# CLEAN
clean:
	rm -f $(SRCDIR)/*.mod $(SRCDIR)/*.o
	rm -f $(BINDIR)/*.o $(BINDIR)/*.mod $(EXEC) *.mod

info:
	@echo ""
	@echo "***********************************************************"
	@echo "***                    AVISM                            ***"
	@echo "***********************************************************"
	@echo "* Óscar Monllor-Berbegal et al., Univ. de València, 2025  *"
	@echo "***********************************************************"
	@echo "*** FLAGS for compiling ***"
	@echo "- COMP (optional, default: 1): 1 is a normal run, 2 is for debugging"
	@echo "- PERIODIC (optional, default: 0): set to 1 if you want periodic boundary conditions"
	@echo "- HDF5 (optional, default: 0): set to 1 if you want to use HDF5 for reading data (e.g. AREPO)"
