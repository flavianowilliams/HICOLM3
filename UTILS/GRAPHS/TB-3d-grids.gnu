set autoscale
reset
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set ticslevel 0.0
set termoption enhanced
set xrange [-2.0:2.0]
set yrange [-2.0:2.0]
set term pngcairo size 600,400 font "Times-Roman,6" enhanced
#set terminal jpeg font "Times-Roman,14"
#set terminal postscript eps color "Times-Roman" 14
#set term wxt size 600,400 "Times-Roman,12"
############################################################
set dgrid3d 250,250
unset hidden3d
set nokey
set contour base
set cntrparam levels incremental -22, 0.05, 6
set pm3d map
#show palette fit2rgbformulae
set palette
#set palette defined (0. "white", 0.1 "black", 1. "black")
############################################################
set output 'TB-3d-valencia.png'
set multiplot layout 2,2
############################################################
#banda 1
col=1
set size 0.45,0.5
set origin 0,0.5
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 2
col=2
set size 0.45,0.5
set origin 0.5,0.5
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 3
col=3
set size 0.45,0.5
set origin 0.,0.
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 4
col=4
set size 0.45,0.5
set origin 0.5,0.
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
############################################################
unset multiplot
set output 'TB-3d-conducao.png'
set multiplot layout 2,2
############################################################
#banda 5
col=5
set size 0.45,0.5
set origin 0,0.5
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 6
col=6
set size 0.45,0.5
set origin 0.5,0.5
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 7
col=7
set size 0.45,0.5
set origin 0.,0.
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
#
#banda 8
col=8
set size 0.45,0.5
set origin 0.5,0.
set title sprintf('Band %d',col)
splot "bandas-3d.dat" u 1:2:col+3 w pm3d
############################################################
unset multiplot
