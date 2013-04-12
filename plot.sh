set term postscript landscape enhanced color 22
set output 'plot.ps'

set grid
set key bottom left
set xlabel 'f_{qq} A_{c}^{y(qq)}'
set ylabel 'f_{qg} A_{c}^{y(qg)}'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_full.txt' u 2:3 pt 6 lc 7 notitle, \
    'data/asymmetry_full_points.txt' w lines lt 4 lc 7 notitle, \
    '' u 3:4 w lines lt 2 lc 7 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 7 title 'full'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_hiM.txt' u 2:3 pt 6 lc 5 notitle, \
    'data/asymmetry_hiM_points.txt' w lines lt 4 lc 5 notitle, \
    '' u 3:4 w lines lt 2 lc 5 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 5 title 'mass>450GeV'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_loM.txt' u 2:3 pt 6 lc 3 notitle, \
    'data/asymmetry_loM_points.txt' w lines lt 4 lc 3 notitle, \
    '' u 3:4 w lines lt 2 lc 3 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 3 title 'mass<450GeV'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_hiY.txt' u 2:3 pt 6 lc 1 notitle, \
    'data/asymmetry_hiY_points.txt' w lines lt 4 lc 1 notitle, \
    '' u 3:4 w lines lt 2 lc 1 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 1 title 'tanhY>0.5'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_loY.txt' u 2:3 pt 6 lc 4 notitle, \
    'data/asymmetry_loY_points.txt' w lines lt 4 lc 4 notitle, \
    '' u 3:4 w lines lt 2 lc 4 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 4 title 'tanhY<0.5'


plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_full_points.txt' u 5:6 w lines lt 1 lw 3 lc 7 title 'full', \
    'data/asymmetry_hiM_points.txt' u 5:6 w lines lt 1 lw 3 lc 5 title 'mass>450GeV', \
    'data/asymmetry_loM_points.txt' u 5:6 w lines lt 1 lw 3 lc 3 title 'mass<450GeV', \
    'data/asymmetry_hiY_points.txt' u 5:6 w lines lt 1 lw 3 lc 1 title 'tanhY>0.5', \
    'data/asymmetry_loY_points.txt' u 5:6 w lines lt 1 lw 3 lc 4 title 'tanhY<0.5'

plot [-0.03:0.03] [-0.007:0.007] \
    'data/asymmetry_full_points.txt' w lines lt 4 lc 7 title 'full', \
    'data/asymmetry_hiM_points.txt' w lines lt 4 lc 5 title 'mass>450GeV', \
    'data/asymmetry_loM_points.txt' w lines lt 4 lc 3 title 'mass<450GeV', \
    'data/asymmetry_hiY_points.txt' w lines lt 4 lc 1 title 'tanhY>0.5', \
    'data/asymmetry_loY_points.txt' w lines lt 4 lc 4 title 'tanhY<0.5'

