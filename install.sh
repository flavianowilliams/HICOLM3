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
fi
#
if [ -f "$exe_dir/hsystem" ]
then
    rm $exe_dir/hsystem
fi
#
if [ -f "$exe_dir/hproperties" ]
then
    rm $exe_dir/hproperties
fi
#
if [ -f "$exe_dir/hftir" ]
then
    rm $exe_dir/hftir
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
# --compilling system--
#
cd $path/contrib/system
#
gfortran system.f90 -o hsystem
#
mv $path/contrib/system/hsystem $exe_dir/hsystem
#
# --compilling properties--
#
cd $path/contrib/properties
#
make clean
make
#
mv $path/contrib/properties/hproperties $exe_dir/hproperties
#
make clean
#
# --compilling ftir--
#
cd $path/contrib/ftir
#
make clean
make
#
mv $path/contrib/ftir/hftir $exe_dir/hftir
#
make clean
#
# --copying auxiliary files--
#
cp $path/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
#
