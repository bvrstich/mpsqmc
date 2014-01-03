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
           HeisenbergMPO.cpp\
           AFMPO.cpp\
           MPStensor.cpp\
           TwoSiteObject.cpp\
           MPSstate.cpp\
           Walker.cpp\
           TrotterHeisenberg.cpp\
           AFQMC.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

THIS_FOLDER = .

INCLUDE = ./include

LIBS = -lblas -llapack

CC  = clang
CXX = clang++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS = -g -I$(INCLUDE) -I$(INCLUDE2) -Wall 
LDFLAGS	= -g -Wall


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
