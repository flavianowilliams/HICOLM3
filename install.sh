#!/bin/bash
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package and utilities
#  version: 2.3.1
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
# -- removing old directories and files
#
if [ -d "$aux_dir/HICOLM" ]
then
    rm -rf $aux_dir/HICOLM
fi
mkdir $aux_dir/HICOLM
mkdir $aux_dir/HICOLM/R
mkdir $aux_dir/HICOLM/R/report
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
echo -e "\e[33m   Preparing to install\e[0m"
echo
echo -e "\e[33m-> Moving files\e[0m"
echo
#
cp -r $path/contrib/amber/*.prm $aux_dir/HICOLM/amber/.
cp -r $path/contrib/R/report/*.R $aux_dir/HICOLM/R/report/.
cp -r $path/contrib/R/report/*.Rmd $aux_dir/HICOLM/R/report/.
#
mv $path/src/HICOLM $exe_dir/HICOLM.bin
mv $path/contrib/ftir/hftir $exe_dir/hftir
mv $path/contrib/properties/hproperties $exe_dir/hproperties
mv $path/contrib/system/hsystem $exe_dir/hsystem
#
# --creating executing script--
#
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
#then
#    $exe_dir/HICOLM.bin | Rscript $aux_dir/HICOLM/R/time_series.R
#else
#    mkdir /home/\$USER/.hicolm
#    cp -r $aux_dir/HICOLM/R /home/\$USER/.hicolm/R
#    cp -r $aux_dir/HICOLM/amber /home/\$USER/.hicolm/amber
#    $exe_dir/HICOLM.bin | Rscript $aux_dir/HICOLM/R/time_series.R
#fi" >> $exe_dir/hicolm
#        ;;
#    no|NO|No|"")
#        echo "
#    $exe_dir/HICOLM.bin
#    ">> $exe_dir/hicolm
#esac
#
chmod +x $exe_dir/hicolm
#
# preparing scripts to get results
#
if [ -f "$exe_dir/hprepare" ]
then
    rm $exe_dir/hprepare
fi
#
echo "#!/bin/sh
#
echo
echo \"Updating R environment...\"
echo
#
if [ -d \"/home/\$USER/.hicolm\" ]
then
    rm -r /home/\$USER/.hicolm
    mkdir /home/\$USER/.hicolm
    cp -r $aux_dir/HICOLM/R /home/\$USER/.hicolm/R
    cp -r $aux_dir/HICOLM/amber /home/\$USER/.hicolm/amber
    echo
    echo \"Finish!\"
    echo
else
    mkdir /home/\$USER/.hicolm
    cp -r $aux_dir/HICOLM/R /home/\$USER/.hicolm/R
    cp -r $aux_dir/HICOLM/amber /home/\$USER/.hicolm/amber
    echo
    echo \"Finish!\"
    echo
fi">> $exe_dir/hprepare
#
chmod +x $exe_dir/hprepare
#
touch $exe_dir/hresults
#
if [ -f "$exe_dir/hresults" ]
then
    rm $exe_dir/hresults
fi
#
touch $exe_dir/hresults
#
echo "#!/bin/sh
#
# - check for auxiliary files and directories
#
if [ ! -d \"/home/\$USER/.hicolm\" ]
then
    echo \"Error to find the auxiliary directory! Running hprepare...\"
    $exe_dir/hprepare
else
    if [ ! -d \"/home/\$USER/.hicolm/R\" ]
    then
        echo \"Error to find the auxiliary directory! Running hprepare...\"
        $exe_dir/hprepare
    else
        if [ ! -d \"/home/\$USER/.hicolm/R/report\" ]
        then
            echo \"Error to find the auxiliary directory! Running hprepare...\"
            $exe_dir/hprepare
        fi
    fi
    if [ ! -d \"/home/\$USER/.hicolm/amber\" ]
    then
        echo \"Error to find the auxiliary directory! Running hprepare...\"
    fi
fi
#
# copying files to auxiliary directories
#
echo
echo 'Please, choose one of the following options:'
echo
echo '1 -> Thermodynamic variables'
echo '2 -> RDF and coordination number (incomplete)'
echo '3 -> Vibrational analysis (incomplete)'
echo
read option
if [ ! -d '1' ]
then
    cp HICOLM.out /home/\$USER/.hicolm/R/report/.
    cp thermodynamics.csv /home/\$USER/.hicolm/R/report/.
    cp atoms.csv /home/\$USER/.hicolm/R/report/.
    Rscript -e \"rmarkdown::render('/home/\$USER/.hicolm/R/report/report.Rmd')\"
    mv /home/\$USER/.hicolm/R/report/report.pdf .
    rm /home/\$USER/.hicolm/R/report/HICOLM.out
    rm /home/\$USER/.hicolm/R/report/thermodynamics.csv
    rm /home/\$USER/.hicolm/R/report/atoms.csv
    rm /home/\$USER/.hicolm/R/report/report.tex
fi" >> $exe_dir/hresults
#
chmod +x $exe_dir/hresults
#
echo -e "\e[32m-> SUCCESS!\e[0m"
echo
echo -e "\e[33m   Installation directory:\e[0m" $exe_dir
echo -e "\e[33m      Auxiliary directory:\e[0m" $aux_dir
echo -e "\e[33m          Compiling rules:\e[0m" $compiler $instructions
echo
echo -e "\e[33m   \"Thank you for choosing HICOLM.\e[0m"
echo -e "\e[33m    To start, just type \e[31mhicolm\e[33m in the terminal.\e[0m"
echo -e "\e[33m    Do not forget the input files HICOLM.in and HICOLM.sys.\\e[0m"
echo -e "\e[33m    HAVE A GREAT JOB!\"\e[0m"
echo
