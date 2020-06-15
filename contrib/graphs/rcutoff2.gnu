reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
set term pngcairo enhanced size 600,400 font "Times-Roman,10"
set output "rcutoff.png"
set multiplot layout 2,1
#set xrange [0.5:3.5]
set xtics 0.5
set size 1,0.5
set style line 1 lc "black" lt 1 lw 1 pt 7 ps 1.2
#####################################################
# GAP x layer
set origin 0,0.5
unset xlabel
set format x ""
#set xlabel 'Graphene layer'
set bmargin 0
set format y "%5.2f"
#set yrange [1.16:1.65]
#set ytics 0.1
set ylabel 'Tempo (s)'
set key title 'Tempo de CPU'
plot  "rcutoff.dat" u 1:2 w lp ls 1 notitle
#####################################################
# Ebind x layer
unset format
unset bmargin
set tmargin 0
set origin 0,0
#set yrange [-8.25:-8.05]
#set ytics 0.04
set xlabel 'r_c (A)'
set format y "%5.2f"
set ylabel 'Energia (eV)'
set key title 'Correção de VdW'
plot  "rcutoff.dat" u 1:3 w lp ls 1 notitle
#####################################################
unset multiplot
