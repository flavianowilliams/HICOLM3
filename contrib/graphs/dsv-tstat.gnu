reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
#set term pngcairo enhanced size 600,400 font "Times-Roman,10"
#set output "dsv-tstat.png"
#set xrange [0:0.6]
#set xtics 1
#####################################################
# GAP x layer
set xlabel 'tstat (ps)'
#set format y "%5.0f"
#set yrange [370:374]
#set ytics 0.1
set ylabel 'Desvio padrão (eV)'
set style line 1 lt 1 lw 1 pt 7 ps 1.2
plot  'dsv-tstat.dat' u 1:2 w lp ls 1 notitle
