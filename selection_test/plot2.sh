set term postscript enhanced color 22 size 8in, 8in
set output 'plot2.ps'
set key bottom

set arrow 1 from 40,1 to 40,0 head size 2,40
set arrow 2 from 100,1 to 100,0 head filled size 2,30
set grid
set key top
set xlabel '{/Symbol e}_{inner} (%)'
set ylabel '{/Symbol e}_{outer} (%)'
plot [0:100] [0:10] 'slopes.txt' \
       u 9:10 w lines ls 2 lc 0 notitle, \
    '' u ($9==100?$9:1/0):10 pt 7 lc 0 ps 2 title '{/Symbol e}_{inner} = 1', \
    '' u ($11>0&&$11<0.001?$9:1/0):10 pt 2 lc 0 ps 2 title '{/Symbol a}_{inner} = 0'

set arrow 1 from 0,0.25 to 0,0 head size 0.06,35
set arrow 2 from 0.05,0.25 to 0.05,0 head filled size 0.06,35
set xlabel '{/Symbol a}_{inner} (%)'
set ylabel '{/Symbol a}_{outer} (%)'
plot [-2:0.5] [0:2.5]\
    '' u 11:12 w lines ls 2 lc 0 notitle, \
    '' u ($9==100?$11:1/0):12 lc 0 pt 7 ps 2 title '{/Symbol e}_{inner} = 1', \
    '' u ($11>0&&$11<0.001?$11:1/0):12 lc 0 pt 2 ps 2 title '{/Symbol a}_{inner} = 0'

set arrow 1 from 40,0.25 to 40,0 head size 2,40
set arrow 2 from 100,0.25 to 100,0 head filled size 2,30
set key right 
set xlabel '{/Symbol e}_{inner} (%)'
set ylabel '({/Symbol a}^{extrapolated} - {/Symbol a}_0) / ({/Symbol a} - {/Symbol a}_0)'
plot [0:100] [0:3] \
    '' u 1:2 w lines lt 2 lc 1 lw 3 notitle, \
    '' u 1:3 w lines lt 2 lc 3 lw 3 notitle, \
    -5 lt 2 lc 0 lw 3 title 'scale {/Symbol a}_{inner}', \
    '' u 1:4 w lines lt 4 lc 1 lw 3 notitle, \
    '' u 1:5 w lines lt 4 lc 3 lw 3 notitle, \
    -5 lt 4 lc 0 lw 3 title 'scale {/Symbol a}_{outer}', \
    '' u 1:6 w lines lt 1 lc 1 lw 3 notitle, \
    '' u 1:7 w lines lt 1 lc 3 lw 3 notitle, \
    -5 lt 1 lc 0 lw 3 title 'scale  both ', \
    -5 lt 1 lc 1 lw 20 title 'unfolding ', \
    -5 lt 1 lc 3 lw 20 title 'templates '

set arrow 1 from 40,-.75 to 40,-1 head size 2,40
set arrow 2 from 100,-.75 to 100,-1 head filled size 2,30
set xlabel '{/Symbol e}_{inner} (%)'
set ylabel '({/Symbol a}^{extrapolated} - {/Symbol a}) / {/Symbol a_0}'
plot [0:100] [-1:2] \
    '' u 1:13 w lines lt 2 lc 1 lw 3 notitle, \
    '' u 1:14 w lines lt 2 lc 3 lw 3 notitle, \
    -500 lt 2 lc 0 lw 3 title '-1 x {/Symbol a}_{inner} ', \
    '' u 1:15 w lines lt 4 lc 1 lw 3 notitle, \
    '' u 1:16 w lines lt 4 lc 3 lw 3 notitle, \
    -500 lt 4 lc 0 lw 3 title '-1 x {/Symbol a}_{outer} ', \
    '' u 1:17 w lines lt 1 lc 1 lw 3 notitle, \
    '' u 1:18 w lines lt 1 lc 3 lw 3 notitle, \
    -500 lt 1 lc 0 lw 3 title '-1 x both   ', \
    -500 lt 1 lc 1 lw 20 title 'unfolding ', \
    -500 lt 1 lc 3 lw 20 title 'templates '
    
