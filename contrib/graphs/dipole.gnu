reset
set termoption enhanced
#set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo enhanced size 600,400 font "Times-Roman,14"
#set term pngcairo size 600,400 font "Times-Roman,10" enhanced
#set output "dipole.png"
set style arrow 1 lc rgb "red" lt 1 lw 1 nohead
set style line 1 lc rgb "black" lt 7 lw 0.25
set style line 2 lc rgb "black" lt 5 lw 1
#set xrange [300:400]
set xlabel "Dipole (D)"
#set yrange [0:1000]
set ylabel "Probability"
set key off

stats 'probability.dat' using 2:3 name "A"
stats 'probability.dat' using 3 name "B"
#set label 1 gprintf("Mean = %g", mean_y) at 2,0.35
set arrow 2 from A_pos_max_y,0 to A_pos_max_y,B_max arrowstyle 1 front

plot 'probability.dat' u 8:9 w l ls 1 title "Dipolo", \
     'probability.dat' u 8:9 ls 2 smooth bezier notitle
