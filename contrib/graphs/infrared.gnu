reset
set termoption enhanced
set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo enhanced size 600,400 font "Times-Roman,14"
set term pngcairo size 600,400 font "Times-Roman,10" enhanced
set output "infrared.png"
set style arrow 1 lc rgb "gray" lw 1 dt 3 nohead
set style arrow 2 lc rgb "black" lw 1 dt 1 filled size screen 0.02,15,45
set style line 1 lw 1 lc rgb "black" dt 1
set style line 2 lw 1 lc rgb "red" dt 1
set style line 3 lw 1 lc rgb "green" dt 1
########################################################################
set multiplot
########################################################################
#plot 1
########################################################################
set origin 0,0
set xrange [0:4000]
set xlabel "Wavenumber (cm^{-1})"
set yrange [0:0.08]
set ylabel "Intensity"
unset ytics
set key on at 3100,0.03
set arrow from 3610,0 to 3610,0.08 arrowstyle 1 front
set arrow from 3733,0 to 3733,0.08 arrowstyle 1 front
set arrow from 1455,0 to 1455,0.08 arrowstyle 1 front
set arrow from 1525,0 to 1525,0.08 arrowstyle 1 front

set arrow from 1525,0.048 to 1525,0.04 arrowstyle 2 front
set label "1525" left at 1525,0.048 font "Times-Roman,8" offset 0.1,0.5
set arrow from 1455,0.058 to 1455,0.05 arrowstyle 2 front
set label "1455" center at 1455,0.058 font "Times-Roman,8" offset 0.1,0.5
set arrow from 3733,0.064 to 3733,0.056 arrowstyle 2 front
set label "3733" center at 3733,0.064 font "Times-Roman,8" offset 0.1,0.5
set arrow from 3610,0.069 to 3610,0.061 arrowstyle 2 front
set label "3610" center at 3610,0.069 font "Times-Roman,8" offset 0.1,0.5

plot "infrared.dat" u 1:3 w l ls 1 title "stretch", \
     "infrared.dat" u 1:4 w l ls 2 title "libration", \
     "infrared.dat" u 1:5 w l ls 3 title "oxygen", \
     "infrared.dat" u 1:6 w l ls 4 title "hydrogen"
########################################################################
#plot 1
########################################################################
set origin 0.5,0.55
set size 0.35,0.35
set key off
unset label
unset arrow
unset xlabel
unset ylabel
unset xrange
set yrange [-1:1]
set xtics 0.2
set ytics 1
plot "vacf.dat" u 1:3 w l ls 1 notitle, \
     "vacf.dat" u 1:4 w l ls 2 notitle, \
     "vacf.dat" u 1:5 w l ls 3 notitle, \
     "vacf.dat" u 1:6 w l ls 4 notitle
########################################################################
unset multiplot
