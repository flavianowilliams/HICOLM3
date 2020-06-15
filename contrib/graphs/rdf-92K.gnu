reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo enhanced size 600,400 font "Times-Roman,14"
set term pngcairo size 600,400 font "Times-Roman,10" enhanced
set output "rdf-92K.png"
#set style arrow 1 lc rgb "red" lt 1 lw 1 nohead
set style line 1 lc rgb "black" lt 7 lw 0.25
set style line 2 lc rgb "black" lt 5 lw 1.5
#set xrange [300:400]
set xlabel "r (A)"
#set yrange [0:1000]
set ylabel "RDF"
#set arrow from 0,371.43 to 100,371.43 arrowstyle 1 front
plot "rdf-Ar-exp.dat" u 1:2 w p ls 1 title "Experimental", \
     "rdf-simulado-92K.dat" u 1:2 w l ls 2 title "Calculado"
