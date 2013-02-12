set term postscript landscape enhanced 22 color
set output 'nll.ps'

#set pm3d map ftriangles interpolate 0,0
set title "ProfileNLL contours"
set xlabel '{/Symbol d}_{qq}'
set ylabel 'R_{ag}'
#set xtics 0.05
set ytics 0.1

set xrange [-1:0]

plot 'nll_1.150.txt' using 2:1 with lines lt 1 lc 1 lw 3 title '1.15', 'nll_3.000.txt' using 2:1 with lines lt 1 lc 3 lw 3 title '3.00'
