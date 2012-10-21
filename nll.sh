set term postscript landscape enhanced 22 color
set output 'nll.ps'

set pm3d map ftriangles interpolate 0,0
set title "-2{/Symbol D}ProfileLL"
set xlabel '{/Symbol d}_{qq}'
set ylabel 'R_{ag}'
set xtics 0.05
set ytics 0.1

set xrange[-1:-0.85]
set yrange[0.3:0.7]

splot 'nll.txt' using 1:2:(2*$3) lc 0 notitle
