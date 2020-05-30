#!/bin/sh
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package
#  version: 3.0.0
#  license: MIT license
#
# --compiling codes--
#
make clean
make
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
# --copying files--
#
cp AMBER/*.prm $aux_dir/HICOLM/amber/.
cp HICOLM $exe_dir/hicolm
#
make clean
