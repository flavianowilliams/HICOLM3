reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
set term png size 600,400 font "Times-Roman,12"
#set term wxt size 600,400 "Times-Roman,12"
set output "gap.png"
set multiplot layout 2,1
set xrange [0.5:3.5]
set xtics 1
set size 1,0.5
#####################################################
# GAP x layer
set origin 0,0.5
unset xlabel
set format x ""
#set xlabel 'Graphene layer'
set bmargin 0
set format y "%5.2f"
set yrange [1.16:1.65]
set ytics 0.1
set ylabel 'Energy (eV)'
set style line 1 lt 1 lw 1 pt 7 ps 1.2
set key title 'GAP'
plot  '-' u 1:2 w lp ls 1 notitle
1 1.61
2 1.34
3 1.22
e
#####################################################
# Ebind x layer
unset format
unset bmargin
set tmargin 0
set origin 0,0
set yrange [-8.25:-8.05]
set ytics 0.04
set xlabel 'Graphene layer'
set format y "%5.2f"
set ylabel 'Energy (eV)'
set style line 1 lt 1 lw 1 pt 7 ps 1.2
set key title 'Binding energy'
plot  '-' u 1:2 w lp ls 1 notitle
1 -8.06
2 -8.20
3 -8.23
e
#####################################################
unset multiplot
