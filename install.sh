#!/bin/sh
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package and utilitaries
#  version: 2.1.0
#  license: MIT license
#
#
path=`pwd`
#
# -- definind installation directory --
#
exe_dir="/usr/local/bin"
aux_dir="/usr/local/share"
#
# -- uninstalling program
#
if [ -f "$exe_dir/hicolm" ]
then
    rm $exe_dir/hicolm
    rm $exe_dir/x2x
fi
#
if [ -d "$aux_dir/HICOLM" ]
then
    rm -r $aux_dir/HICOLM
fi
#
# --creating directory--
#
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/amber
#
# --compilling HICOLM--
#
cd $path/src
#
make clean
make
#
mv $path/src/HICOLM $exe_dir/hicolm
#
make clean
#
# --compilling x2x--
#
cd $path/contrib/x2x
#
gfortran x2x.f90 -o x2x
#
mv $path/contrib/x2x/x2x $exe_dir/x2x
#
# --compilling properties--
#
cd $path/contrib/properties
#
make clean
make
#
mv $path/contrib/properties/properties $exe_dir/properties
#
make clean
#
# --copying auxiliary files--
#
cp $path/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
#
