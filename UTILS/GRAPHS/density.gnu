reset
set termoption enhanced
set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
#set term pngcairo size 600,400 font "Times-Roman,12"
#set term wxt size 600,400 "Times-Roman,12"
#set output "densidade.png"
#set xrange [0.5:3.5]
#set xtics 1
set style line 1 lc rgb "black" lt 5 lw 1
set style line 2 lc rgb "red" lt 5 lw 1
set style line 3 lc rgb "green" lt 4 lw 1
#####################################################
#set format x ""
set xlabel 'Temperatura (K)'
#set format y "%5.2f"
#set yrange [1.16:1.65]
#set ytics 0.1
set ylabel 'Densidade (kg/m^3)'
set style line 1 lt 1 lw 1 pt 7 ps 1.2
plot "densidade-exp.dat" u 1:2 w lp ls 1 title "Experimental", \
     "densidade-mod1.dat" u 1:2 w lp ls 2 title "Simulado", \
     "densidade-mod2.dat" u 1:2 w lp ls 3 title "Simulado"
