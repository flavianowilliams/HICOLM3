reset
set termoption enhanced
#set encoding utf8
set term post eps solid size 6,8 font "Times-Roman,20"
set output "pdos-1layer.eps"
#set term png size 600,400 font "Times-Roman,14" enhanced
#set term wxt size 600,800 "Times-Roman,12"
set style arrow 1 lt 0 lw 4 nohead
set style line 2 lc rgb "black" lt 1 lw 4
set style line 3 lc rgb "red"   lt 1 lw 4
set style line 4 lc rgb "green" lt 1 lw 4
set xrange [-4:4]
set ylabel "PDOS"
set mytics 2
###################################################################
set multiplot layout 2,1
#
#configurando grafico PDOS do zinco
#
set format x ""
set yrange [0:1000]
set bmargin 0
set size 1,0.25
set origin 0,0.75
set arrow from 0,0 to 0,1000 arrowstyle 1
plot "PDOS-Zn1-1layer.dat" u ($1*27+3.07827):2 w l ls 2 title "Zn_{s}", \
     "PDOS-Zn1-1layer.dat" u ($1*27+3.07827):3 w l ls 3 title "Zn_{p}", \
     "PDOS-Zn1-1layer.dat" u ($1*27+3.07827):4 w l ls 4 title "Zn_{d}"
#
unset format
unset bmargin
unset yrange
unset arrow
#
set xlabel "E-E_{F} (eV)"
set yrange [1000:0]
set tmargin 0
set size 1,0.25
set origin 0,0.50
set arrow from 0,0 to 0,1000 arrowstyle 1
plot "zno_xo_DS2_DOS_AT0001" u ($1*27-2.9832):2 w l ls 2 notitle, \
     "zno_xo_DS2_DOS_AT0001" u ($1*27-2.9832):3 w l ls 3 notitle, \
     "zno_xo_DS2_DOS_AT0001" u ($1*27-2.9832):4 w l ls 4 notitle
#
unset tmargin
unset arrow
unset yrange
unset xlabel
#####################################################################
#configurando grafico PDOS do oxigenio
#
set format x ""
set yrange [0:2500]
set bmargin 0
set size 1,0.25
set origin 0,0.25
set arrow from 0,0 to 0,2500 arrowstyle 1
plot "PDOS-O1-1layer.dat" u ($1*27+3.07827):2 w l ls 2 title "Zn_{s}", \
     "PDOS-O1-1layer.dat" u ($1*27+3.07827):3 w l ls 3 title "Zn_{p}", \
     "PDOS-O1-1layer.dat" u ($1*27+3.07827):4 w l ls 4 title "Zn_{d}"
#
unset format
unset bmargin
unset yrange
unset xlabel
unset arrow
#
set xlabel "E-E_{F} (eV)"
set yrange [2500:0]
set tmargin 0
set size 1,0.25
set origin 0,0.0
set arrow from 0,0 to 0,2500 arrowstyle 1
plot "zno_xo_DS2_DOS_AT0003" u ($1*27-2.9832):2 w l ls 2 notitle, \
     "zno_xo_DS2_DOS_AT0003" u ($1*27-2.9832):3 w l ls 3 notitle, \
     "zno_xo_DS2_DOS_AT0003" u ($1*27-2.9832):4 w l ls 4 notitle
#
unset multiplot
#####################################################################