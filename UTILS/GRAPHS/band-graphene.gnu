reset
###################################################################
#Gerando output
###################################################################
set termoption enhanced
set encoding utf8
set term pngcairo size 600,400 font "Times-Roman,10" enhanced
#set term wxt size 600,400 "Times-Roman,12"
#set term post eps solid size 6,4 font "Times-Roman,20" enhanced
set output "TB.png"
###################################################################
#Parametros gerais!!!
###################################################################
set style arrow 1 lc rgb "black"   lt 0 lw 5.0 nohead
set style line 1 lc rgb   "red"    lt 1 lw 1.0
set style line 2 lc rgb "black"    pt 4 ps 0.5
set style line 3 lc rgb "red"      lt 2 lw 1.5 ps 1.5  dashtype 1
####################################################################
#grafico estrutura de bandas
####################################################################
#
set xrange [0.0:2.2]
set yrange [-20.0:20.0]
#set xlabel "k"
set ylabel "Energy (eV)"
set ytics 4.0
set grid xtics ls 1
set xtics ("K" 0.0, '{/Symbol \107}' 0.90, "M" 1.68, "K" 2.13)
#set arrow from 0.0 to 2.2 arrowstyle 1
#set title 'Grafeno'
#
plot  for [col=2:10] "TB.dat" u 1:col w l ls 1 notitle, \
      for [col=2:8] "bandas.dat" every 5 u 1:col w p ls 2 notitle
####################################################################
