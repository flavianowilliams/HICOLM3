reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
set term pngcairo enhanced size 600,400 font "Times-Roman,10"
set output "rcutoff.png"
set multiplot layout 2,1
set xrange [4:10]
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
set format y "%5.0f"
set yrange [70:250]
set ytics 50
set ylabel 'Tempo (s)'
set label 'Tempo de CPU' center at graph 0.8,0.2
plot  "rcutoff.dat" u 1:2 w lp ls 1 notitle
#####################################################
# Ebind x layer
unset format
unset bmargin
unset label
set tmargin 0
set origin 0,0
set yrange [-2:0]
set ytics 0.5
set xlabel 'r_c (A)'
set format y "%5.1f"
set ylabel 'Energia (eV)'
set label 'Correção de VdW' center at graph 0.8,0.2
plot  "rcutoff.dat" u 1:3 w lp ls 1 notitle
#####################################################
unset multiplot
