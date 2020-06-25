#!/bin/bash
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package and utilitaries
#  version: 2.1.0
#  license: MIT license
#
path=`pwd`
#
# -- definind installation and auxiliary directory --
#
echo
echo "Please, type the compiler or press ENTER (default: gfortran)"
read compiler
#
if [ -z $compiler ]
then
    compiler="gfortran"
fi
#
export FC="$compiler"
#
echo
echo "Please, type the instructions of compilation or press ENTER"
read instructions
#
if [ -z $instructions ]
then
#    instructions="-fcheck=all -fbacktrace -Wall"
    instructions=""
fi
#
export FFLAGS="$instructions"
#
echo
echo "Please, type the installation directory or press ENTER (default: /usr/local/bin)"
read exe_dir
#
if [ -z $exe_dir ]
then
    exe_dir="/usr/local/bin"
fi
#
echo
echo "Please, type the auxiliary directory or press ENTER (default: /usr/local/share)"
read aux_dir
#
if [ -z $aux_dir ]
then
    aux_dir="/usr/local/share"
fi
#
# --creating directories
#
if [ -d "$aux_dir/HICOLM" ]
then
    rm -rf $aux_dir/HICOLM
fi
#
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/amber
#
# -- installing HICOLM --
#
echo
echo -e "\e[33m-> Compiling HICOLM\e[0m"
echo
#
if [ -f "$exe_dir/HICOLM.bin" ]
then
    rm -f $exe_dir/HICOLM.bin
fi
#
cd $path/src
#
if [ -f "HICOLM" ]
then
    rm HICOLM
fi
make -s clean
make -s all
if [ ! -f "HICOLM" ]
then
    echo -e "\e[31mError in compiling HICOLM. The installation will be finish!"
    exit
fi
make -s clean
#
# -- installing hsystem --
#
echo -e "\e[33m-> Compiling utilities\e[0m"
echo
#
if [ -f "$exe_dir/hsystem" ]
then
    rm $exe_dir/hsystem
fi
#
cd $path/contrib/system
if [ -f "hsystem" ]
then
    rm hsystem
fi
$compiler system.f90 $instructions -o hsystem
#
if [ ! -f "hsystem" ]
then
    echo
    echo -e "\e[31mError in compiling hsystem. The installation will be finish!"
    echo
    exit
fi
#
# -- installing hproperties --
#
if [ -f "$exe_dir/hproperties" ]
then
    rm $exe_dir/hproperties
fi
#
cd $path/contrib/properties
#
if [ -f "hproperties" ]
then
    rm hproperties
fi
make -s clean
make -s all
if [ ! -f "hproperties" ]
then
    echo
    echo -e "\e[31mError in compiling hproperties. The installation will be finish!"
    echo
    exit
fi
make -s clean
#
# -- installing hftir --
#
if [ -f "$exe_dir/hftir" ]
then
    rm $exe_dir/hftir
fi
#
cd $path/contrib/ftir
#
if [ -f "hftir" ]
then
    rm hftir
fi
make -s clean
make -s all
if [ ! -f "hftir" ]
then
    echo
    echo -e "\e[31mError in compiling hftir. The installation will be finish!"
    echo
    exit
fi
make -s clean
#
# --copying auxiliary files--
#
echo -e "\e[33m-> Copying auxiliary files\e[0m"
echo
#
cp $path/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
#
echo -e "\e[33m-> Preparing to install\e[0m"
echo
echo -e "\e[33m   Installation directory:\e[0m" $exe_dir
echo -e "\e[33m   Auxiliary directory:\e[0m" $aux_dir
echo
mv $path/src/HICOLM $exe_dir/HICOLM.bin
mv $path/contrib/ftir/hftir $exe_dir/hftir
mv $path/contrib/properties/hproperties $exe_dir/hproperties
mv $path/contrib/system/hsystem $exe_dir/hsystem
#
# --creating executing script--
#
if [ -f "$exe_dir/hicolm" ]
then
    rm $exe_dir/hicolm
fi
#
touch $exe_dir/hicolm
#
echo "#!/bin/sh
if [ ! -d '/tmp/amber' ]
then
    cp -r $aux_dir/HICOLM/amber /tmp/amber
fi
$exe_dir/HICOLM.bin" >> $exe_dir/hicolm
#
chmod +x $exe_dir/hicolm
#
echo -e "\e[32m-> SUCCESS!\e[0m"
echo
echo -e "\e[33m   To start, just type \e[31mhicolm\e[33m in the terminal.\e[0m"
echo
