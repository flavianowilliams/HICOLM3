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
# -- definind installation and auxiliary directory --
#
#echo
#echo "Please, type the installation directory or presse ENTER to accept the default (/usr/local/bin)"
#read exe_dir
#
if [ -z $exe_dir ]
then
    exe_dir="/usr/local/bin"
fi
#
#echo "Please, type the installation directory or presse ENTER to accept the default option (/usr/local/share)"
#read aux_dir
#
if [ -z $aux_dir ]
then
    aux_dir="/usr/local/share"
fi
#
# --creation directories
#
if [ -d "$aux_dir/HICOLM" ]
then
    rm -r $aux_dir/HICOLM
fi
#
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/amber
#
# -- installing HICOLM --
#
echo
echo "\e[33mInstalling HICOLM in $exe_dir\e[0m"
echo
#
if [ -f "$exe_dir/hicolm" ]
then
    rm $exe_dir/hicolm
fi
#
cd $path/src
make clean
make
if [ ! -f "HICOLM" ]
then
    echo "\e[31mError in compiling HICOLM. The installation is going to finish!"
    exit
fi
mv $path/src/HICOLM $exe_dir/hicolm
make clean
#
# -- installing hsystem --
#
echo
echo "\e[33mInstalling utilitaries in $exe_dir\e[0m"
echo
#
if [ -f "$exe_dir/hsystem" ]
then
    rm $exe_dir/hsystem
fi
#
cd $path/contrib/system
gfortran system.f90 -o hsystem
if [ ! -f "hsystem" ]
then
    echo "\e[31mError in compiling hsystem. The installation is going to finish!"
    exit
fi
mv $path/contrib/system/hsystem $exe_dir/hsystem
#
# -- installing hproperties --
#
if [ -f "$exe_dir/hproperties" ]
then
    rm $exe_dir/hproperties
fi
#
cd $path/contrib/properties
make clean
make
if [ ! -f "hproperties" ]
then
    echo "\e[31mError in compiling hproperties. The installation is going to finish!"
    exit
fi
mv $path/contrib/properties/hproperties $exe_dir/hproperties
make clean
#
# -- installing hftir --
#
if [ -f "$exe_dir/hftir" ]
then
    rm $exe_dir/hftir
fi
#
cd $path/contrib/ftir
make clean
make
if [ ! -f "hftir" ]
then
    echo "\e[31mError in compiling hftir. The installation is going to finish!"
    exit
fi
mv $path/contrib/ftir/hftir $exe_dir/hftir
make clean
#
# --copying auxiliary files--
#
echo
echo "\e[33mCopying auxiliary files to $aux_dir\e[0m"
echo
#
cp $path/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
#
echo "\e[32mSUCCESS!\e[0m"
echo
echo "\e[32mTo start, just type \e[31mhicolm\e[32m in the terminal.\e[0m"
echo
