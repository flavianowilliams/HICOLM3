reset
set termoption enhanced
set encoding utf8
#set term post eps solid size 6,4 font "Times-Roman,20"
#set term pdfcairo solid size 6,4 font "Roman,12"
set term pngcairo size 600,900 font "Times-Roman,10"
#set term wxt size 600,400 "Times-Roman,12"
set output "densidade-press.png"
set style line 1 lc rgb "black" lt 1 lw 1 pt 7 ps 1.2
set style line 2 lc rgb "red" lt 7 lw 1
#####################################################
set xrange [50:1200]
set xtics 100
set yrange [10:1200]
set ytics 150
set ylabel 'Densidade (kg/m^3)'
unset xlabel
set multiplot layout 4,1
set origin 0,0.75
set bmargin 0
set tmargin 1
set label 'r_c = 10 A' at graph 0.05,0.75
set label '373 K' at graph 0.85,0.7
set label '298 K' at graph 0.85,0.82
set label '273 K' at graph 0.85,0.9
unset xlabel
set format x ""
plot "densidade-273K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-298K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-373K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-273K-10A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-298K-10A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-373K-10A.dat" u 1:2:3 w yerrorbars ls 2 notitle
set origin 0,0.5
set bmargin 0
set tmargin 0
unset label
set label 'r_c = 8 A' at graph 0.05,0.75
set label '373 K' at graph  0.85,0.7
set label '298 K' at graph 0.85,0.82
set label '273 K' at graph 0.85,0.9
plot "densidade-273K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-298K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-373K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-273K-8A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-298K-8A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-373K-8A.dat" u 1:2:3 w yerrorbars ls 2 notitle
set origin 0,0.25
#set bmargin 1
#set tmargin 0
unset label
set label 'r_c = 6 A' at graph 0.05,0.75
set label '373 K' at graph  0.85,0.7
set label '298 K' at graph 0.85,0.82
set label '273 K' at graph 0.85,0.9
plot "densidade-273K-exp.dat" u 1:2 w lp ls 1 notitle, \
     "densidade-298K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-373K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-273K-6A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-298K-6A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-373K-6A.dat" u 1:2:3 w yerrorbars ls 2 notitle
set origin 0,0.0
#set tmargin 0
set bmargin 3
unset label
set label 'r_c = 4 A' at graph 0.05,0.75
set label '373 K' at graph  0.85,0.7
set label '298 K' at graph 0.85,0.82
set label '273 K' at graph 0.85,0.9
set xlabel 'Pressão (atm)'
set format x "%5.0f"
plot "densidade-273K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-298K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-373K-exp.dat" u 1:2 w lp ls 1 notitle,\
     "densidade-273K-4A.dat" u 1:2:3 w yerrorbars ls 2 notitle,\
     "densidade-298K-4A.dat" u 1:2:3 w yerrorbars ls 2 notitle
unset multiplot
