set autoscale
reset
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ticslevel 0.0
set termoption enhanced
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
#set zrange [-2:2]
set ztics 1.0
set xlabel "k_{x}"
set ylabel "k_{y}"
set zlabel "E(k)"
set term pngcairo size 600,400 font "Times-Roman,10" enhanced
#set terminal jpeg font "Times-Roman,14"
#set terminal postscript eps color "Times-Roman" 14
#set term wxt size 600,400 "Times-Roman,12"
set dgrid3d 250,250
set hidden3d
#set view 90,0
set nokey
set contour base
set cntrparam levels incremental -22, 0.05, 6
#set pm3d map
#show palette fit2rgbformulae
set palette
#set palette defined (0. "white", 0.1 "black", 1. "black")
do for [col=4:4] {
   set title sprintf('Band %d',col)
   set output sprintf('bandas-3d%d.png',col)
   splot "bandas-3d.dat" u 1:2:col+3 w pm3d}
