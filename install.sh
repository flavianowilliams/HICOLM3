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
echo
echo "Please, type the compiler or press ENTER (default: gfortran)"
read compiler
#
if [ -z $compiler ]
then
    compiler="gfortran"
fi
#
echo
echo "Please, type the instructions of compilation or press ENTER"
read instructions
#
if [ -z $instructions ]
then
    instructions="-fcheck=all -fbacktrace -Wall"
#    instructions=""
fi
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
    rm $exe_dir/HICOLM.bin
fi
#
cd $path/src
#
sed -i -e "s/.*FC=.*/FC=$compiler/" makefile
sed -i -e "s/.*FFLAGS=.*/FFLAGS=$instructions/" makefile
#
make clean
make all
if [ ! -f "HICOLM" ]
then
    echo "\e[31mError in compiling HICOLM. The installation will be finish!"
    exit
fi
mv $path/src/HICOLM $exe_dir/HICOLM.bin
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
$compiler system.f90 $instructions -o hsystem
#
if [ ! -f "hsystem" ]
then
    echo "\e[31mError in compiling hsystem. The installation will be finish!"
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
#
sed -i -e "s/.*FC=.*/FC=$compiler/" makefile
sed -i -e "s/.*FFLAGS=.*/FFLAGS=$instructions/" makefile
#
make clean
make all
if [ ! -f "hproperties" ]
then
    echo "\e[31mError in compiling hproperties. The installation will be finish!"
    exit
fi
mv $path/contrib/properties/hproperties $exe_dir/hproperties
make clean
echo
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
sed -i -e "s/.*FC=.*/FC=$compiler/" makefile
sed -i -e "s/.*FFLAGS=.*/FFLAGS=$instructions/" makefile
#
make clean
make all
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
echo "\e[32mSUCCESS!\e[0m"
echo
echo "\e[32mTo start, just type \e[31mhicolm\e[32m in the terminal.\e[0m"
echo
