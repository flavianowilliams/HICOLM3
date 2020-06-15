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
cd $path/SRC
#
make clean
make
#
# --compilling x2x--
#
cd $path/UTILS/x2x
#
gfortran x2x.f90 -o x2x
#
# --copying files--
#
cp $path/CONTRIB/AMBER/*.prm $aux_dir/HICOLM/amber/.
cp $path/SRC/HICOLM $exe_dir/hicolm
cp $path/CONTRIB/x2x/x2x $exe_dir/x2x
#
# --cleaning directories
#
