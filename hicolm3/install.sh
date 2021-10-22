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
if [ -d '/tmp/hicolm3' ]
then
    rm -r /tmp/hicolm3
fi
#
mkdir /tmp/hicolm3
#
if [ ! -d '/tmp/hicolm3/charmm' ]
then
    cp -r $aux_dir/HICOLM3/charmm /tmp/hicolm3/charmm
else
    if [ ! -f '/tmp/hicolm3/charmm/charmm_bonds.prm' ]
    then
        cp -r $aux_dir/HICOLM3/charmm/charmm_bonds.prm /tmp/hicolm3/charmm/charmm_bonds.prm
    fi
    if [ ! -f '/tmp/hicolm3/charmm/charmm_angles.prm' ]
    then
        cp -r $aux_dir/HICOLM3/charmm/charmm_angles.prm /tmp/hicolm3/charmm/charmm_angles.prm
    fi
    if [ ! -f '/tmp/hicolm3/charmm/charmm_dihedrals.prm' ]
    then
        cp -r $aux_dir/HICOLM3/charmm/charmm_dihedrals.prm /tmp/hicolm3/charmm/charmm_dihedrals.prm
    fi
    if [ ! -f '/tmp/hicolm3/charmm/charmm_idihedrals.prm' ]
    then
        cp -r $aux_dir/HICOLM3/charmm/charmm_idihedrals.prm /tmp/hicolm3/charmm/charmm_idihedrals.prm
    fi
    if [ ! -f '/tmp/hicolm3/charmm/charmm_vdw.prm' ]
    then
        cp -r $aux_dir/HICOLM3/charmm/charmm_vdw.prm /tmp/hicolm3/charmm/charmm_vdw.prm
    fi
fi
if [ ! -d '/tmp/hicolm3/opls' ]
then
    cp -r $aux_dir/HICOLM3/opls /tmp/hicolm3/opls
else
    if [ ! -f '/tmp/hicolm3/opls/opls_vdw.prm' ]
    then
        cp -r $aux_dir/HICOLM3/opls/opls_vdw.prm /tmp/hicolm3/opls/opls_vdw.prm
    fi
    if [ ! -f '/tmp/hicolm3/opls/opls_bonds.prm' ]
    then
        cp -r $aux_dir/HICOLM3/opls/opls_bonds.prm /tmp/hicolm3/opls/opls_bonds.prm
    fi
    if [ ! -f '/tmp/hicolm3/opls/opls_angles.prm' ]
    then
        cp -r $aux_dir/HICOLM3/opls/opls_angles.prm /tmp/hicolm3/opls/opls_angles.prm
    fi
fi
#
$exe_dir/HICOLM3.bin" >> $exe_dir/hicolm3
#
chmod +x $exe_dir/hicolm3
#
