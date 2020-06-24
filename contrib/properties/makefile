# Compiler
#
#FC=gfortran
#
# Instructions
#
#FFLAGS=-fcheck=all -fbacktrace -Wall
#
# libraries
#
LDFLAGS=
#=====================================================================================
OBJ_SRC=input_module.o variaveis_module.o rdf_module.o
#
OBJ_MAIN=main.o
#=====================================================================================
#
$(OBJ_SRC):
	$(FC) -c $(FFLAGS) $*.f90
#
$(OBJ_MAIN):
	$(FC) -c $(FFLAGS) $*.f90
#
all: $(OBJ_SRC) $(OBJ_MAIN)
	$(FC) -o hproperties $(LDFLAGS) $(OBJ_SRC) $(OBJ_MAIN)
#
clean:
	rm -f *.o *.mod hproperties
#==========================================================================
