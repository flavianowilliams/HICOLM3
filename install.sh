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
#echo
#echo "Please, type the auxiliary directory or press ENTER (default: /usr/local/share)"
#read aux_dir
#
if [ -z $aux_dir ]
then
    aux_dir="/usr/local/share"
fi
#
# -- removing old directories and files
#
if [ -d "$aux_dir/HICOLM3" ]
then
    rm -rf $aux_dir/HICOLM
fi
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/hicolm
mkdir $aux_dir/HICOLM/hicolm/R
mkdir $aux_dir/HICOLM/hicolm/R/report
mkdir $aux_dir/HICOLM/hicolm/amber
#
# -- installing HICOLM --
#
echo
echo -e "\e[33m-> Compiling HICOLM\e[0m"
echo
#
# -- installing older version of HICOLM
#
sh ./hicolm/install.sh
#
