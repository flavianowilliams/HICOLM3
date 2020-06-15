reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
set term pngcairo enhanced size 600,400 font "Times-Roman,10"
set output "tstat.png"
set style line 1 lc "black" lt 1 lw 1 pt 7 ps 1.0
set style line 2 lc "black" lt 1 lw 1 pt 7 ps 1.0
set multiplot
#####################################################
# desvio padrao da energia
set xrange [0.0005:1]
set origin 0,0
set size 1,1
set xlabel '{/Symbol t}_T (ps)'
set logscale x 10
set format y "%5.2f"
set yrange [0.05:0.35]
#set ytics 0.1
set samples 1000
set ylabel 'Flutuação na energia total (eV)'
plot  'tstat.dat' u 1:2 w p ls 1 notitle,\
      'tstat.dat' u 1:2 smooth sbezier ls 1 notitle
#####################################################
# temperatura media
unset xlabel
set origin 0.5,0.55
set size 0.4,0.4
set format y "%5.0f"
set yrange [0:15]
set ytics 3
set ylabel '{/Symbol D}T (K)' offset 3,0
plot  'tstat.dat' u 1:3 w lp ls 2 notitle
unset multiplot
