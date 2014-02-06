###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = mpsqmc
CPPSRC	= main.cpp\
           Random.cpp\
           Operator.cpp\
           OpSx.cpp\
           OpSy.cpp\
           OpSz.cpp\
           OpI.cpp\
           Op0.cpp\
           MPO.cpp\
           AFMPO.cpp\
           MPStensor.cpp\
           TwoSiteObject.cpp\
           MPSstate.cpp\
           Walker.cpp\
           TrotterJ1J2.cpp\
           J1J2MPO.cpp\
           HamMPO.cpp\
           AFQMC.cpp\
           WorkSpace.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

THIS_FOLDER = .

INCLUDE = ./include

LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread

CC  = icc
CXX = icpc

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS = -I$(INCLUDE) -O3 -ipo -openmp
LDFLAGS	= -O3 -ipo -openmp


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(THIS_FOLDER)/$(BINNAME)
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(THIS_FOLDER)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(THIS_FOLDER)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(THIS_FOLDER)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

#-----------------------------------------------------------------------------
# Make the documentation
#----------------------------------------------------------------------------
doc:
	@doxygen .doc-config

# ====================== End of file 'makefile.in' ========================== #
