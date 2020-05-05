EX = gaffields
TYPE = comp
#===============================================================
OBJ_MOD = input.o populacao.o objetivo.o selecao.o utils.o \
	operadores.o
OBJ_OUT = outfile.o
OBJ_SRC = main.o
#===============================================================
gfortran:
	$(MAKE) LD="gfortran -o" FC="gfortran" \
	FFLAGS="-c -fcheck=all -Wall" \
	LDFLAGS="" \
	EX=$(EX) $(TYPE)
#
comp:	$(OBJ_MOD) $(OBJ_OUT) $(OBJ_SRC)
	$(LD) $(EX) $(OBJ_MOD) $(OBJ_OUT) $(OBJ_SRC) $(LDFLAGS)
#
clean:
	rm -f $(EX) *.mod *.o *.dat
#=====================================================================
# Declare dependencies
%.o:
	$(FC) $(FFLAGS) $*.f90
