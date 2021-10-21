#=====================================================================================
OBJ_UTILS = constants_module.o
#
OBJ_SYSTEM = structure_module.o molecule_module.o zmatrix_module.o system_module.o
#
OBJ_TOPOLOGY = charmm_module.o opls_module.o forcefield_module.o
#
OBJ_PREPARE = prepare_module.o
#
OBJ_MD = input_module.o atoms_module.o neighbourlist_module.o bonds_module.o \
	angles_module.o idihedrals_module.o dihedrals_module.o coulomb_module.o \
	vanderwaals_module.o interaction_module.o thermodynamics_module.o \
	ensemble_module.o moleculardynamics_module.o
#
OBJ_OPTIMIZE = interopt_module.o gradientdescent_module.o optimize_module.o
#
OBJ_MAIN = main.o
#=====================================================================================
ODIR1=mm/utils
ODIR2=mm/system
ODIR3=mm/topology
ODIR4=mm/prepare
ODIR5=mm/md
ODIR6=mm/optimize
ODIR7=mm/main
#
$(OBJ_UTILS):
	$(FC) -c $(FFLAGS) $(ODIR1)/$*.f08
#
$(OBJ_SYSTEM):
	$(FC) -c $(FFLAGS) $(ODIR2)/$*.f08
#
$(OBJ_TOPOLOGY):
	$(FC) -c $(FFLAGS) $(ODIR3)/$*.f08
#
$(OBJ_PREPARE):
	$(FC) -c $(FFLAGS) $(ODIR4)/$*.f08
#
$(OBJ_MD):
	$(FC) -c $(FFLAGS) $(ODIR5)/$*.f08
#
$(OBJ_OPTIMIZE):
	$(FC) -c $(FFLAGS) $(ODIR6)/$*.f08
#
$(OBJ_MAIN):
	$(FC) -c $(FFLAGS) $(ODIR7)/$*.f08
#
all: $(OBJ_UTILS) $(OBJ_SYSTEM) $(OBJ_TOPOLOGY) $(OBJ_PREPARE) $(OBJ_MD) $(OBJ_OPTIMIZE) $(OBJ_MAIN)
	$(FC) -o HICOLM3 $(LDFLAGS) $(OBJ_UTILS) $(OBJ_SYSTEM) $(OBJ_TOPOLOGY) $(OBJ_PREPARE) $(OBJ_MD) $(OBJ_OPTIMIZE) $(OBJ_MAIN)
#
clean:
#	rm -f $(ODIR1)/*.o $(ODIR2)/*.o $(ODIR3)/*.o $(ODIR4)/*.o  $(ODIR5)/*.o  $(ODIR6)/*.o $(ODIR7)/*.o  *.o *.mod *.md
	rm -f $(OBJ_MD) $(OBJ_OPTIMIZE) $(OBJ_MAIN) *.md
#==========================================================================
