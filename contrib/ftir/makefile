#=====================================================================================
OBJ_SRC=sizes.o error.o system.o utils.o cmmolec.o utilsrange.o dipmol.o corr.o tcf.o \
	ftir.o printfile.o
#
OBJ_MAIN=ftirclass.o
#=====================================================================================
#
$(OBJ_SRC):
	$(FC) -c $(FFLAGS) $*.f
#
$(OBJ_MAIN):
	$(FC) -c $(FFLAGS) $*.f
#
all: $(OBJ_SRC) $(OBJ_MAIN)
	$(FC) -o hftir $(LDFLAGS) $(OBJ_SRC) $(OBJ_MAIN)
#
clean:
	rm -f *.o *.mod hftir
#==========================================================================
