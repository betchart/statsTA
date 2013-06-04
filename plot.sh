set term postscript enhanced color 22 size 7in, 7in
set output 'plot.ps'

set bmargin 3
set tmargin 3
set grid x2tics ytics
set key bottom right
set key samplen 1.5 font "Helvetica, 14"
set x2label 'A@_{c}^{y(qq)}  (%)' offset 0,-0.7
set ylabel 'A@_{c}^{y(qg)}  (%)' offset 1,0
set xlabel 'A@_{c}^{y}  (%)' offset 0,0.8


qq = 0.003953110899140473
qg = 0.0021230524809500123
ep = 0.0001
S = 100

set ytics 0.01
set mytics 2
set format y "%.0t"
set format x2 "%.0t"
set x2tics -0.03,0.01,0.03 offset 0,-0.5
set xtics ("-3" -0.02, "-2" -0.01, "-1" 0, "0" 0.01, "1" 0.02, "2" 0.03, "3" 0.04) rotate by 15
set size square

hor=0.04
ver=0.01
off = hor+0.005
theta = atan2(ver**2,hor**2)

altox(x,o) = (tan(theta)*o + x) / (1+tan(theta))
altoy(x,o) = (-tan(theta)*(o-x)) / (1+tan(theta))

altx(x) = altox(x,hor)
alty(x) = altoy(x,hor)
aplot(x) = (x-hor)*tan(theta)
tic(x,o) = (o-x<aplot(x)?o-x:1/0)

set style fill solid 1.0 noborder
set sample 1001

col=7
plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_full.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'Full Selection', \
    'data/asymmetry_full_points.txt' w lines lt 4 lc col title 'PLL contour: 1.14', \
    '' u 3:4 w lines lt 2 lc col title 'Systematic Unc.', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CL', \
    '' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 5:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($5+$6,off)):(altoy($5+$6,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($7+$8,off)):(altoy($7+$8,off)) w lines lt 1 lw 3 lc 2 notitle

col=5
plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_hiM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'm_{tt} > 450 GeV', \
    'data/asymmetry_hiM_points.txt' w lines lt 4 lc col title 'PLL contour: 1.14', \
    '' u 3:4 w lines lt 2 lc col title 'Systematic Unc.', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CL', \
    '' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 5:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($5+$6,off)):(altoy($5+$6,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($7+$8,off)):(altoy($7+$8,off)) w lines lt 1 lw 3 lc 2 notitle

col=3
plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_loM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'm_{tt} < 450 GeV', \
    'data/asymmetry_loM_points.txt' w lines lt 4 lc col title 'PLL contour: 1.14', \
    '' u 3:4 w lines lt 2 lc col title 'Systematic Unc.', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CL', \
    '' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 5:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($5+$6,off)):(altoy($5+$6,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($7+$8,off)):(altoy($7+$8,off)) w lines lt 1 lw 3 lc 2 notitle

col=1
plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_hiY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'tanh|y_{tt}| > 0.5', \
    'data/asymmetry_hiY_points.txt' w lines lt 4 lc col title 'PLL contour: 1.14', \
    '' u 3:4 w lines lt 2 lc col title 'Systematic Unc.', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CL', \
    '' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 5:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($5+$6,off)):(altoy($5+$6,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($7+$8,off)):(altoy($7+$8,off)) w lines lt 1 lw 3 lc 2 notitle

col=4
plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_loY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'tanh|y_{tt}| < 0.5', \
    'data/asymmetry_loY_points.txt' w lines lt 4 lc col title 'PLL contour: 1.14', \
    '' u 3:4 w lines lt 2 lc col title 'Systematic Unc.', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CL', \
    '' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 5:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($5+$6,off)):(altoy($5+$6,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($7+$8,off)):(altoy($7+$8,off)) w lines lt 1 lw 3 lc 2 notitle


lwidth=5
ep = 1.5*ep

set key title '68% CL' spacing 0.85

plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    'data/asymmetry_full_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 7:(ver-5*ep) w lines lt 1 lc 2 lw 50 notitle, \
    '' u (-hor+5*ep*hor/ver):8 w lines lt 1 lc 2 lw 50 notitle, \
    '' u (altox($7+$8,off-10*ep)):(altoy($7+$8,off-10*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($7+$8,off+20*ep)):(altoy($7+$8,off+20*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($7+$8,off+50*ep)):(altoy($7+$8,off+50*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($7+$8,off+80*ep)):(altoy($7+$8,off+80*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 7 title 'Full Selection', \
    '' u 5:(ver-9*ep) w lines lt 1 lc 7 lw lwidth notitle, \
    '' u (-hor+9*ep*hor/ver):6 w lines lt 1 lc 7 lw lwidth notitle, \
    '' u (altox($5+$6,off+80*ep)):(altoy($5+$6,off+80*ep)) w lines lt 1 lc 7 lw lwidth notitle, \
    'data/asymmetry_hiM_points.txt' u 5:6 w lines lt 1 lw 3 lc 5 title 'm_{tt} > 450 GeV', \
    '' u 5:(ver-7*ep) w lines lt 1 lc 5 lw lwidth notitle, \
    '' u (-hor+7*ep*hor/ver):6 w lines lt 1 lc 5 lw lwidth notitle, \
    '' u (altox($5+$6,off+50*ep)):(altoy($5+$6,off+50*ep)) w lines lt 1 lc 5 lw lwidth notitle, \
    'data/asymmetry_loM_points.txt' u 5:6 w lines lt 1 lw 3 lc 3 title 'm_{tt} < 450 GeV', \
    '' u 5:(ver-5*ep) w lines lt 1 lc 3 lw lwidth notitle, \
    '' u (-hor+5*ep*hor/ver):6 w lines lt 1 lc 3 lw lwidth notitle, \
    '' u (altox($5+$6,off+20*ep)):(altoy($5+$6,off+20*ep)) w lines lt 1 lc 3 lw lwidth notitle, \
    'data/asymmetry_hiY_points.txt' u 5:6 w lines lt 1 lw 3 lc 1 title 'tanh|y_{tt}| > 0.5', \
    '' u 5:(ver-3*ep) w lines lt 1 lc 1 lw lwidth notitle, \
    '' u (-hor+3*ep*hor/ver):6 w lines lt 1 lc 1 lw lwidth notitle, \
    '' u (altox($5+$6,off-10*ep)):(altoy($5+$6,off-10*ep)) w lines lt 1 lc 1 lw lwidth notitle, \
    'data/asymmetry_loY_points.txt' u 5:6 w lines lt 1 lw 3 lc 4 title 'tanh|y_{tt}| < 0.5', \
    '' u 5:(ver-1*ep) w lines lt 1 lc 4 lw lwidth notitle, \
    '' u (-hor+1*ep*hor/ver):6 w lines lt 1 lc 4 lw lwidth notitle, \
    '' u (altox($5+$6,off-10*ep)):(altoy($5+$6,off-10*ep)) w lines lt 1 lc 4 lw lwidth notitle
