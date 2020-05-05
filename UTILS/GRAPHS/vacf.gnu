reset
set termoption enhanced
set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo enhanced size 600,400 font "Times-Roman,14"
#set term pngcairo size 600,400 font "Times-Roman,10" enhanced
#set output "vacf.png"
set style line 1 lw 1 lc rgb "black" dt 1
set style line 2 lw 1 lc rgb "red" dt 1
set style line 3 lw 1 lc rgb "green" dt 1
#set xrange [0:0.1]
set xlabel "Time (ps)"
#set yrange [0:1000]
set ylabel "VACF"
plot "vacf.dat" u 1:2 w l ls 1 title "stretch", \
     "vacf.dat" u 1:3 w l ls 2 title "stretch", \
     "vacf.dat" u 1:4 w l ls 3 title "libration", \
     "vacf.dat" u 1:5 w l ls 3 title "1", \
     "vacf.dat" u 1:6 w l ls 3 title "2"
