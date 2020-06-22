EX = hproperties
TYPE = comp
#===============================================================
OBJ_MOD = input_module.o variaveis_module.o rdf_module.o
OBJ_SRC = main.o
#===============================================================
gfortran:
	$(MAKE) LD="gfortran -o" FC="gfortran" \
	FFLAGS="-c -fcheck=all -Wall" \
	LDFLAGS="" \
	EX=$(EX) $(TYPE)
#
comp:	$(OBJ_MOD) $(OBJ_SRC)
	$(LD) $(EX) $(OBJ_MOD) $(OBJ_SRC) $(LDFLAGS)
#
clean:
	rm -f $(EX) *.mod *.o
#=====================================================================
# Declare dependencies
%.o:
	$(FC) $(FFLAGS) $*.f90
