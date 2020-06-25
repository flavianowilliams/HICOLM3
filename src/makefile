#=====================================================================================
OBJ_SYS=sistema.o elements.o utils.o amber.o input.o alloc_arrays.o estrutura.o \
	brillouin.o
#
OBJ_MD=bonds_module.o bends_module.o dihedral_module.o \
	vdw_module.o coulomb_module.o neighbour_list.o \
	force_field.o ensemble_nve.o ensemble_nvt.o ensemble_npt.o molecular_dynamics.o
#
OBJ_OPT=optimize.o
#
OBJ_MAIN=outfile.o main.o
#=====================================================================================
ODIR1=main
ODIR2=md
ODIR3=opt
ODIR4=main
#
$(OBJ_SYS):
	$(FC) -c $(FFLAGS) $(ODIR1)/$*.f90
#
$(OBJ_MD):
	$(FC) -c $(FFLAGS) $(ODIR2)/$*.f90
#
$(OBJ_OPT):
	$(FC) -c $(FFLAGS) $(ODIR3)/$*.f90
#
$(OBJ_MAIN):
	$(FC) -c $(FFLAGS) $(ODIR4)/$*.f90
#
all: $(OBJ_SYS) $(OBJ_MD) $(OBJ_OPT) $(OBJ_MAIN)
	$(FC) -o HICOLM $(LDFLAGS) $(OBJ_SYS) $(OBJ_MD) $(OBJ_OPT) $(OBJ_MAIN)
#
clean:
	rm -f $(ODIR1)/*.o $(ODIR2)/*.o $(ODIR3)/*.o $(ODIR4)/*.o *.o *.mod *.dat *.out *.md
#==========================================================================
