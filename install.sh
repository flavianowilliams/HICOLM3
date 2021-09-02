#!/bin/bash
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package and utilities
#  version: 3.0.0
#  license: MIT license
#
path=`pwd`
#
# -- definind installation and auxiliary directory --
#
#echo
#echo "Please, type the compiler or press ENTER (default: gfortran)"
#read compiler
#
if [ -z $compiler ]
then
    compiler="gfortran"
fi
#
export FC="$compiler"
#
#echo
#echo "Please, type the instructions of compilation or press ENTER"
#read instructions
#
if [ -z $instructions ]
then
    instructions="-fcheck=all -fbacktrace -Wall"
#    instructions=""
fi
#
export FFLAGS="$instructions"
#
#echo
#echo "Please, type the installation directory or press ENTER (default: /usr/local/bin)"
#read exe_dir
#
if [ -z $exe_dir ]
then
    exe_dir="/usr/local/bin"
fi
#
export exe_dir=$exe_dir
#
#echo
#echo "Please, type the auxiliary directory or press ENTER (default: /usr/local/share)"
#read aux_dir
#
if [ -z $aux_dir ]
then
    aux_dir="/usr/local/share"
fi
#
export aux_dir=$aux_dir
#
# -- removing old directories and files
#
if [ -d "$aux_dir/HICOLM" ]
then
    rm -rf $aux_dir/HICOLM
fi
#
# -- creating new directories
#
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/amber
#
# --copying auxiliary files--
#
echo
echo -e "\e[33m   Preparing to install\e[0m"
echo
echo -e "\e[33m-> Moving files\e[0m"
#
cp -r $path/hicolm/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
#
#sh ./hicolm/install.sh # -- HICOLM install script
#
<<<<<<< HEAD
#echo -e "\e[33m-> Preparing R environment\e[0m"
#Rscript $path/contrib/R/prepare.R
#echo
#
# --preparing script to call HICOLM executable
#
if [ -f "$exe_dir/hicolm" ]
then
    rm $exe_dir/hicolm
fi
#
touch $exe_dir/hicolm
#
echo "#!/bin/sh
#
if [ ! -d '/tmp/amber' ]
then
    cp -r $aux_dir/HICOLM/amber /tmp/amber
else
    if [ ! -f '/tmp/amber/amber_bonds.prm' ]
    then
        cp -r $aux_dir/HICOLM/amber/amber_bonds.prm /tmp/amber/amber_bonds.prm
    fi
    if [ ! -f '/tmp/amber/amber_angles.prm' ]
    then
        cp -r $aux_dir/HICOLM/amber/amber_angles.prm /tmp/amber/amber_angles.prm
    fi
    if [ ! -f '/tmp/amber/amber_dihedrals_general.prm' ]
    then
        cp -r $aux_dir/HICOLM/amber/amber_dihedrals_general.prm /tmp/amber/amber_dihedrals_general.prm
    fi
    if [ ! -f '/tmp/amber/amber_dihedrals_proper.prm' ]
    then
        cp -r $aux_dir/HICOLM/amber/amber_dihedrals_proper.prm /tmp/amber/amber_dihedrals_proper.prm
    fi
    if [ ! -f '/tmp/amber/amber_vdw.prm' ]
    then
        cp -r $aux_dir/HICOLM/amber/amber_vdw.prm /tmp/amber/amber_vdw.prm
    fi
fi
#
$exe_dir/HICOLM.bin" >> $exe_dir/hicolm
#
#case "$supp" in
#    yes|YES|Yes)
#        echo "
#if [ -d \"/home/\$USER/.hicolm\" ]
=======
#if [ ! -f "$path/hicolm/src/HICOLM" ]
>>>>>>> development
#then
#    echo -e "\e[31mError in compiling HICOLM. Installation aborted!"
#    exit
#fi
sh ./hicolm3/install.sh # -- HICOLM3 install script
#
if [ ! -f "$path/hicolm3/src/HICOLM3" ]
then
    echo -e "\e[31mError in compiling HICOLM3. Installation aborted!"
    exit
fi
#
echo -e "\e[32m-> SUCCESS!\e[0m"
echo
echo -e "\e[33m   Installation directory:\e[0m" $exe_dir
echo -e "\e[33m      Auxiliary directory:\e[0m" $aux_dir
echo -e "\e[33m          Compiling rules:\e[0m" $compiler $instructions
echo
echo -e "\e[33m   \"Thank you for choosing HICOLM3.\e[0m"
echo -e "\e[33m    To start the older version, just type \e[31mhicolm\e[33m in the terminal, and to start the newest one type \e[31mhicolm3\e[33m.\e[0m"
echo -e "\e[33m    Warning! The older version is deprecated and it will be maintained until march, 2022.\\e[0m"
echo -e "\e[33m    HAVE A GREAT JOB!\"\e[0m"
echo
#
