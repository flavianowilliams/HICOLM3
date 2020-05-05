reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
#set term pngcairo enhanced size 600,400 font "Times-Roman,10"
#set output "tstat.png"
set style line 1 lc "black" lt 1 lw 1 pt 7 ps 1.2
set style line 2 lc "black" lt 1 lw 1 pt 7 ps 1.2
set multiplot
#####################################################
# desvio padrao da energia
set xrange [0.0005:1]
set origin 0,0
set size 1,1
set xlabel 'pstat (ps)'
set logscale x
#set xtics 0.01
#set format y "%5.0f"
#set yrange [370:374]
#set ytics 0.1
set ylabel 'Flutuação na energia total (eV)'
plot  'pstat.dat' u 1:2 w lp ls 1 notitle,\
      'pstat.dat' u 1:2 smooth cspline notitle
#####################################################
# temperatura media
unset xlabel
set origin 0.2,0.5
set size 0.4,0.4
set format y "%5.0f"
set yrange [4:10]
set ytics 1
set ylabel '{/Symbol D}T (K)' offset 3,0
plot  'pstat.dat' u 1:3 w lp ls 2 notitle
unset multiplot
