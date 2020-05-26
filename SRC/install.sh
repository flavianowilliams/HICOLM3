#!/bin/sh
#
# --compiling codes--
#
make clean
make
#
# --creating directory--
#
mkdir /usr/local/share/HICOLM
mkdir /usr/local/share/HICOLM/amber
#
# --copying files--
#
cp HICOLM /usr/local/bin/hicolm
#
cp AMBER/*.prm /usr/local/share/HICOLM/amber/.
#
make clean
