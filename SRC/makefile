# Diretorio de instalacao do executavel
ROOT=.
#
# Compilador
FC=gfortran
#
# Instrucoes do compilador
FFLAGS=-fcheck=all -fbacktrace -Wall
#
# bibliotecas
LDFLAGS=
#=====================================================================================
OBJ_SYS = sistema.o elements.o utils.o amber.o input.o alloc_arrays.o estrutura.o \
		brillouin.o
#
OBJ_MD = vdw_module.o coulomb_module.o neighbour_list.o \
		bonds_module.o bends_module.o dihedral_module.o \
		force_field.o ensemble_nve.o \
		ensemble_nvt.o ensemble_npt.o molecular_dynamics.o
#
OBJ_MAIN = outfile.o main.o
#=====================================================================================
ODIR1=$(ROOT)/MAIN
ODIR2=$(ROOT)/MD
ODIR3=$(ROOT)/MAIN
#
OBJ1 = $(patsubst %,$(ODIR1)/%,$(OBJ_SYS))
OBJ2 = $(patsubst %,$(ODIR2)/%,$(OBJ_MD))
OBJ3 = $(patsubst %,$(ODIR3)/%,$(OBJ_MAIN))
#
$(ODIR1)/%.o:
		$(FC) -c $(FFLAGS) -o $@ $(ODIR1)/$*.f90
#
$(ODIR2)/%.o:
		$(FC) -c $(FFLAGS) -o $@ $(ODIR2)/$*.f90
#
$(ODIR3)/%.o:
		$(FC) -c $(FFLAGS) -o $@ $(ODIR3)/$*.f90
#
HICOLM:	$(OBJ1) $(OBJ2) $(OBJ3)
		$(FC) -o $@ $^ $(LDFLAGS)
#
clean:
		rm -f $(ODIR1)/*.o $(ODIR2)/*.o $(ODIR3)/*.o *.o *.mod *.dat *.out *.md HICOLM
#==========================================================================
