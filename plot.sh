set term postscript landscape enhanced color 22
set output 'plot.ps'

set grid
set key bottom left
set xlabel 'f_{qq} A_{c}^{y(qq)}'
set ylabel 'f_{qg} A_{c}^{y(qg)}'


qq = 0.003953110899140473
qg = 0.0021230524809500123
ep = 0.0001

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_full.txt' u 2:3 pt 6 lc 7 notitle, \
    'data/asymmetry_full_points.txt' w lines lt 4 lc 7 notitle, \
    '' u 3:4 w lines lt 2 lc 7 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 7 title 'full', \
    '' u 7:8 w lines lt 1 lc 7 notitle, \
    '' u 9:10 w lines lt 1 lw 4 lc 7 notitle, \
    '' u 11:12 w lines lt 1 lc 2 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 2 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_hiM.txt' u 2:3 pt 6 lc 5 notitle, \
    'data/asymmetry_hiM_points.txt' w lines lt 4 lc 5 notitle, \
    '' u 3:4 w lines lt 2 lc 5 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 5 title 'mass>450GeV', \
    '' u 7:8 w lines lt 1 lc 7 notitle, \
    '' u 9:10 w lines lt 1 lw 4 lc 5 notitle, \
    '' u 11:12 w lines lt 1 lc 2 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 2 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_loM.txt' u 2:3 pt 6 lc 3 notitle, \
    'data/asymmetry_loM_points.txt' w lines lt 4 lc 3 notitle, \
    '' u 3:4 w lines lt 2 lc 3 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 3 title 'mass<450GeV', \
    '' u 7:8 w lines lt 1 lc 7 notitle, \
    '' u 9:10 w lines lt 1 lw 4 lc 3 notitle, \
    '' u 11:12 w lines lt 1 lc 2 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 2 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_hiY.txt' u 2:3 pt 6 lc 1 notitle, \
    'data/asymmetry_hiY_points.txt' w lines lt 4 lc 1 notitle, \
    '' u 3:4 w lines lt 2 lc 1 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 1 title '|tanh y|>0.5', \
    '' u 7:8 w lines lt 1 lc 7 notitle, \
    '' u 9:10 w lines lt 1 lw 4 lc 1 notitle, \
    '' u 11:12 w lines lt 1 lc 2 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 2 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_loY.txt' u 2:3 pt 6 lc 4 notitle, \
    'data/asymmetry_loY_points.txt' w lines lt 4 lc 4 notitle, \
    '' u 3:4 w lines lt 2 lc 4 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 4 title '|tanh y|<0.5', \
    '' u 7:8 w lines lt 1 lc 7 notitle, \
    '' u 9:10 w lines lt 1 lw 4 lc 4 notitle, \
    '' u 11:12 w lines lt 1 lc 2 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 2 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'

eps=0.0003
off=0.001

set arrow from qq+qg,0 to qq+qg,-(off+6*eps) lc 2 nohead

plot [-0.035:0.035] [-0.01:0.01] \
    'data/asymmetry_full_points.txt' u 5:6 w lines lt 1 lw 3 lc 7 title 'full', \
    '' u 9:($10-1*eps-off) w lines lt 1 lw 7 lc 7 notitle, \
    'data/asymmetry_hiM_points.txt' u 5:6 w lines lt 1 lw 3 lc 5 title 'mass>450GeV', \
    '' u 9:($10-2*eps-off) w lines lt 1 lw 7 lc 5 notitle, \
    'data/asymmetry_loM_points.txt' u 5:6 w lines lt 1 lw 3 lc 3 title 'mass<450GeV', \
    '' u 9:($10-3*eps-off) w lines lt 1 lw 7 lc 3 notitle, \
    'data/asymmetry_hiY_points.txt' u 5:6 w lines lt 1 lw 3 lc 1 title '|tanh y|>0.5', \
    '' u 9:($10-4*eps-off) w lines lt 1 lw 7 lc 1 notitle, \
    'data/asymmetry_loY_points.txt' u 5:6 w lines lt 1 lw 3 lc 4 title '|tanh y|<0.5', \
    '' u 9:($10-4*eps-off) w lines lt 1 lw 7 lc 4 notitle, \
    (x+ep >= qq && x <= qg + qq ? qq + qg - x : 1/0) lc 2 lt 1 notitle, \
    "<echo '0.003953110899140473 0.0021230524809500123'" w points ls 1 lc 2 ps 3 lw 4 title 'POWHEG CT10'
