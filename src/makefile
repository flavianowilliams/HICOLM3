#=====================================================================================
OBJ_GEN = constants_module.o structure_module.o system_module.o molecule_module.o
#
OBJ_PRP = zmatrix_module.o amber_module.o forcefield_module.o prepare_module.o
#
OBJ_MD = input_module.o atoms_module.o neighbourlist_module.o interaction_module.o \
				moleculardynamics_module.o
#
OBJ_OPT = optimize.o
#
OBJ_MAIN = main.o
#=====================================================================================
ODIR1=main
ODIR2=prepare
ODIR3=md
ODIR4=opt
#
$(OBJ_GEN):
	$(FC) -c $(FFLAGS) $(ODIR1)/$*.f08
#
$(OBJ_PRP):
	$(FC) -c $(FFLAGS) $(ODIR2)/$*.f08
#
$(OBJ_MD):
	$(FC) -c $(FFLAGS) $(ODIR3)/$*.f08
#
$(OBJ_OPT):
	$(FC) -c $(FFLAGS) $(ODIR4)/$*.f08
#
$(OBJ_MAIN):
	$(FC) -c $(FFLAGS) $(ODIR1)/$*.f08
#
#all: $(OBJ_SYS) $(OBJ_MD) $(OBJ_OPT) $(OBJ_MAIN)
#	$(FC) -o HICOLM $(LDFLAGS) $(OBJ_SYS) $(OBJ_MD) $(OBJ_OPT) $(OBJ_MAIN)
all: $(OBJ_GEN) $(OBJ_PRP) $(OBJ_MD) $(OBJ_MAIN)
	$(FC) -o HICOLM $(LDFLAGS) $(OBJ_GEN) $(OBJ_PRP) $(OBJ_MD) $(OBJ_MAIN)
#
clean:
	rm -f $(ODIR1)/*.o $(ODIR2)/*.o $(ODIR3)/*.o $(ODIR4)/*.o *.o *.mod *.dat *.out *.md
#==========================================================================
