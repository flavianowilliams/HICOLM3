#!/bin/bash
#
path=`pwd`
path=$path/hicolm3
#
# -- installing HICOLM3 --
#
echo
echo "\e[33m-> Compiling HICOLM3 routines\e[0m"
echo
if [ -f "$exe_dir/HICOLM3.bin" ]
then
    rm -f $exe_dir/HICOLM3.bin
fi
#
cd $path/src
#
if [ -f "HICOLM3" ]
then
    rm HICOLM3
fi
make -s clean
make -s all
if [ ! -f "HICOLM3" ]
then
    exit
fi
make -s clean
#
# --copying auxiliary files--
#
echo "\e[33m   Preparing to install\e[0m"
echo
echo "\e[33m-> Moving files\e[0m"
echo
#
cp $path/src/HICOLM3 $exe_dir/HICOLM3.bin
#
# --creating executing script--
#
#echo -e "\e[33m-> Preparing R environment\e[0m"
#Rscript $path/contrib/R/prepare.R
#
# --preparing script to call HICOLM executable
#
if [ -f "$exe_dir/hicolm3" ]
then
    rm $exe_dir/hicolm3
fi
#
touch $exe_dir/hicolm3
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
$exe_dir/HICOLM3.bin" >> $exe_dir/hicolm3
#
chmod +x $exe_dir/hicolm3
#
