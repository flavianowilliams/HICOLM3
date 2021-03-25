#!/bin/bash
#
#   author: Flaviano Williams Fernandes <flaviano.fernandes@ifpr.edu.br>
# describe: Installation of HICOLM package and utilities
#  version: 3.0.0
#  license: MIT license
#
path=`pwd`
path=$path/hicolm
#
# -- creating new directories
#
mkdir $aux_dir/HICOLM/hicolm
mkdir $aux_dir/HICOLM/hicolm/R
mkdir $aux_dir/HICOLM/hicolm/R/report
#
# -- installing hsystem --
#
#
echo
echo "\e[33m-> Compiling HICOLM routines\e[0m"
echo
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
    exit
fi
make -s clean
#
# -- installing hsystem --
#
echo "\e[33m-> Compiling utilities\e[0m"
echo
#
if [ -f "$exe_dir/hsystem" ]
then
    rm $exe_dir/hsystem
fi
#
cd $path/contrib/system
#
if [ -f "hsystem" ]
then
    rm hsystem
fi
$FC system.f90 $FFLAGS -o hsystem
#
if [ ! -f "hsystem" ]
then
    echo
    echo "\e[31mError in compiling hsystem. The installation will be finish!"
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
    echo "\e[31mError in compiling hproperties. The installation will be finish!"
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
#cd $path/contrib/ftir
#
#if [ -f "hftir" ]
#then
#    rm hftir
#fi
#make -s clean
#make -s all
#if [ ! -f "hftir" ]
#then
#    echo
#    echo -e "\e[31mError in compiling hftir. The installation will be finish!"
#    echo
#    exit
#fi
#make -s clean
#
cp $path/src/HICOLM $exe_dir/HICOLM.bin
cp -r $path/contrib/R/report/*.R $aux_dir/HICOLM/hicolm/R/report/.
cp -r $path/contrib/R/report/*.Rmd $aux_dir/HICOLM/hicolm/R/report/.
#mv $path/contrib/ftir/hftir $exe_dir/hftir
mv $path/contrib/properties/hproperties $exe_dir/hproperties
mv $path/contrib/system/hsystem $exe_dir/hsystem
#
# --creating executing script--
#
#echo -e "\e[33m-> Preparing R environment\e[0m"
#Rscript $path/contrib/R/prepare.R
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
