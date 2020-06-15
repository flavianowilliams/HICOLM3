EX = hftir
TYPE = seq
SHELL=/bin/bash
#===============================================================
MOD_MAIN = sizes.o error.o system.o utils.o cmmolec.o
MOD_MOD = utilsrange.o dipmol.o corr.o tcf.o ftir.o printfile.o
MOD_SRC = ftirclass.o
#===============================================================
intel:
	$(MAKE) LD="gfortran -o" FC="gfortran" \
	FFLAGS="-fcheck=all -Wall" \
	EX=$(EX) $(TYPE)
#
seq:	$(MOD_MAIN) $(MOD_MOD) $(MOD_SRC)
	$(LD) $(EX) $(MOD_MAIN) $(MOD_MOD) $(MOD_SRC)
#
clean:
	rm -f $(MOD_MAIN) $(MOD_MOD) $(MOD_SRC) $(EX) *.dat *.out *.mod
#=====================================================================
# Declare dependencies
.o:
	$(SHELL) $(FC) $(FFLAGS) $*.f
