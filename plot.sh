set term postscript enhanced color 22 size 8in, 8in
set output 'graphics/plot.ps'

set bmargin 1
set tmargin 1
set rmargin 4
set lmargin 6
set grid x2tics ytics
set key invert
set key top right
set key width -20
set key samplen 1.5 font "Helvetica, 14"
set x2label 'A@_{c}^{y(qq)}  (%)' offset 0,-0.7
set ylabel 'A@_{c}^{y(qg)}  (%)' offset 1,0
set xlabel 'A@_{c}^{y}  (%)' offset 0,0.8


kr2012 = 0.0102
ep = 0.0001
S = 100

set ytics 0.01
set mytics 2
set format y "%.0t"
set format x2 "%.0t"
set x2tics -0.03,0.01,0.03 offset 0,-0.5
set xtics ("-3" -0.02, "-2" -0.01, "-1" 0, "0" 0.01, "1" 0.02, "2" 0.03, "3" 0.04) rotate by 25
set size square

hor=0.02
ver=0.01
off = hor+0.002
theta = atan2(ver**2,hor**2)

altox(x,o) = (tan(theta)*o + x) / (1+tan(theta))
altoy(x,o) = (-tan(theta)*(o-x)) / (1+tan(theta))

altx(x) = altox(x,hor)
alty(x) = altoy(x,hor)
aplot(x) = (x-hor)*tan(theta)
tic(x,o) = (o-x<aplot(x)?o-x:1/0)
above(x,y)= y>aplot(x)?y:1/0

set linetype 5 lc rgb "#008B8B"
set style fill solid 1.0 noborder
set sample 1001

itemize = 1

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
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_full_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u (itemize?$15:1/0):(itemize?$16:1/0) w lines lt 4 lc col title 'Sys (MC stat)', \
    '' u (itemize?$3:$17):(itemize?$4:$18) w lines lt 2 lc col title 'Systematic', \
    '' w lines lt 1 lc col title 'PLL = 1.14', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CI', \
    '' u 9:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($11,off)):(altoy($11,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($14,off)):(altoy($14,off)) w lines lt 1 lw 3 lc 2 notitle, \
    'output/asymmetry_full.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'Full Selection'

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
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_hiM_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u (itemize?$15:1/0):(itemize?$16:1/0) w lines lt 4 lc col title 'Sys (MC stat)', \
    '' u (itemize?$3:$17):(itemize?$4:$18) w lines lt 2 lc col title 'Systematic', \
    '' w lines lt 1 lc col title 'PLL = 1.14', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CI', \
    '' u 9:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($11,off)):(altoy($11,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($14,off)):(altoy($14,off)) w lines lt 1 lw 3 lc 2 notitle, \
    'output/asymmetry_hiM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'm_{tt} > 450 GeV'

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
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_loM_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u (itemize?$15:1/0):(itemize?$16:1/0) w lines lt 4 lc col title 'Sys (MC stat)', \
    '' u (itemize?$3:$17):(itemize?$4:$18) w lines lt 2 lc col title 'Systematic', \
    '' w lines lt 1 lc col title 'PLL = 1.14', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CI', \
    '' u 9:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($11,off)):(altoy($11,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($14,off)):(altoy($14,off)) w lines lt 1 lw 3 lc 2 notitle, \
    'output/asymmetry_loM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'm_{tt} < 450 GeV'

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
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_hiY_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u (itemize?$15:1/0):(itemize?$16:1/0) w lines lt 4 lc col title 'Sys (MC stat)', \
    '' u (itemize?$3:$17):(itemize?$4:$18) w lines lt 2 lc col title 'Systematic', \
    '' w lines lt 1 lc col title 'PLL = 1.14', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CI', \
    '' u 9:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($11,off)):(altoy($11,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($14,off)):(altoy($14,off)) w lines lt 1 lw 3 lc 2 notitle, \
    'output/asymmetry_hiY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'tanh|y_{tt}| > 0.5'

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
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_loY_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u (itemize?$15:1/0):(itemize?$16:1/0) w lines lt 4 lc col title 'Sys (MC stat)', \
    '' u (itemize?$3:$17):(itemize?$4:$18) w lines lt 2 lc col title 'Systematic', \
    '' w lines lt 1 lc col title 'PLL = 1.14', \
    '' u 5:6 w lines lt 1 lw 3 lc col title '68% CI', \
    '' u 9:(ver-5*ep) w lines lt 1 lc col lw 10 notitle, \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 3 notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc col lw 10 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 3 notitle, \
    '' u (altox($11,off)):(altoy($11,off)) w lines lt 1 lw 10 lc col notitle, \
    '' u (altox($14,off)):(altoy($14,off)) w lines lt 1 lw 3 lc 2 notitle, \
    'output/asymmetry_loY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc col title 'tanh|y_{tt}| < 0.5'


lwidth=5
ep = 0.5*ep

set key spacing 0.8
unset key


plot [-hor:hor] [-ver:ver] \
    aplot(x) w filledcurves x1 lc rgb "white" notitle, \
    aplot(x) lt 1 lc rgb "white" lw 0 title '68% CI', \
    aplot(x) lt 1 lc 7 lw 3 notitle, \
    tic(x,-0.03) notitle lt 4 lc 7, \
    tic(x,-0.02) notitle lt 4 lc 7, \
    tic(x,-0.01) notitle lt 4 lc 7, \
    tic(x,0) notitle lt 4 lc 7, \
    tic(x,0.01) notitle lt 4 lc 7, \
    tic(x,0.02) notitle lt 4 lc 7, \
    tic(x,0.03) notitle lt 4 lc 7, \
    above(x, kr2012 - x) lt 1 lc rgb "#E0E0E0" lw 11 title '{}_{KR2012}', \
    'output/asymmetry_full_points.txt' u 7:8 w filledcu closed lc 2 title '{}_{POWHEG-CT10}', \
    '' u 12:(ver-5*ep) w lines lt 1 lc 2 lw 50 notitle, \
    '' u (-hor+5*ep*hor/ver):13 w lines lt 1 lc 2 lw 50 notitle, \
    '' u (altox($14,off-10*ep)):(altoy($14,off-10*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+20*ep)):(altoy($14,off+20*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+50*ep)):(altoy($14,off+50*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+80*ep)):(altoy($14,off+80*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u 5:6 w lines lt 1 lw 3 lc 7 title 'Full Selection', \
    '' u 9:(ver-9*ep) w lines lt 1 lc 7 lw lwidth notitle, \
    '' u (-hor+9*ep*hor/ver):10 w lines lt 1 lc 7 lw lwidth notitle, \
    '' u (altox($11,off+80*ep)):(altoy($11,off+80*ep)) w lines lt 1 lc 7 lw lwidth notitle, \
    'output/asymmetry_hiM_points.txt' u 5:6 w lines lt 1 lw 3 lc 5 title 'm_{tt} > 450 GeV', \
    '' u 9:(ver-7*ep) w lines lt 1 lc 5 lw lwidth notitle, \
    '' u (-hor+7*ep*hor/ver):10 w lines lt 1 lc 5 lw lwidth notitle, \
    '' u (altox($11,off+50*ep)):(altoy($11,off+50*ep)) w lines lt 1 lc 5 lw lwidth notitle, \
    'output/asymmetry_loM_points.txt' u 5:6 w lines lt 1 lw 3 lc 3 title 'm_{tt} < 450 GeV', \
    '' u 9:(ver-5*ep) w lines lt 1 lc 3 lw lwidth notitle, \
    '' u (-hor+5*ep*hor/ver):10 w lines lt 1 lc 3 lw lwidth notitle, \
    '' u (altox($11,off+20*ep)):(altoy($11,off+20*ep)) w lines lt 1 lc 3 lw lwidth notitle, \
    'output/asymmetry_hiY_points.txt' u 5:6 w lines lt 1 lw 3 lc 1 title 'tanh|y_{tt}| > 0.5', \
    '' u 9:(ver-3*ep) w lines lt 1 lc 1 lw lwidth notitle, \
    '' u (-hor+3*ep*hor/ver):10 w lines lt 1 lc 1 lw lwidth notitle, \
    '' u (altox($11,off-10*ep)):(altoy($11,off-10*ep)) w lines lt 1 lc 1 lw lwidth notitle, \
    'output/asymmetry_loY_points.txt' u 5:6 w lines lt 1 lw 3 lc 4 title 'tanh|y_{tt}| < 0.5', \
    '' u 9:(ver-1*ep) w lines lt 1 lc 4 lw lwidth notitle, \
    '' u (-hor+1*ep*hor/ver):10 w lines lt 1 lc 4 lw lwidth notitle, \
    '' u (altox($11,off-10*ep)):(altoy($11,off-10*ep)) w lines lt 1 lc 4 lw lwidth notitle


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
    above(x, kr2012 - x) lt 1 lc rgb "gray" lw 1 title '{}_{KR2012}', \
    'output/asymmetry_full.txt' u ($0==2?$2:1/0):3 pt 7 ps 0.6 lc 2 title '{}_{POWHEG-CT10}', \
    '' u ($0==3?$2:1/0):3 pt 6 ps 0.6 lc 2 notitle, \
    above(x,-x+0.181988641018*0.0116658564583+0.13387790945*0.0295313129862) lc 2 lt 1 notitle, \
    above(x,-x+0.158720638333*0.0107352062128+0.144493371919*0.0301769084344) lc 2 lt 3 notitle, \
    '' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc 0 title 'Full Selection', \
    '' u ($0==1?$2:1/0):3 pt 6 ps 0.6 lc 0 notitle , \
    above(x,-x+0.00392834766936+0.00371431842306) lc 0 lt 1 notitle, \
    above(x,-x+0.00470930384515+0.0030793670305) lc 0 lt 3 notitle, \
    'output/asymmetry_hiM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc 5 title 'm_{tt} > 450 GeV', \
    '' u ($0==1?$2:1/0):3 pt 6 ps 0.6 lc 5 notitle , \
    above(x,-x+0.00441966295901+0.00461882402026) lc 5 lt 1 notitle, \
    above(x,-x+0.0045882695171+0.00394332174142) lc 5 lt 3 notitle, \
    'output/asymmetry_loM.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc 3 title 'm_{tt} < 450 GeV', \
    '' u ($0==1?$2:1/0):3 pt 6 ps 0.6 lc 3 notitle , \
    above(x,-x-0.00258507695697+0.00329948404313) lc 3 lt 1 notitle, \
    above(x,-x-0.00302884613032+0.00294388451244) lc 3 lt 3 notitle, \
    'output/asymmetry_hiY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc 1 title 'tanh|y_{tt}| > 0.5', \
    '' u ($0==1?$2:1/0):3 pt 6 ps 0.6 lc 1 notitle , \
    above(x,-x+-0.000555240911894+0.00225455451716) lc 1 lt 1 notitle, \
    above(x,-x+-0.000108439852522+0.00198727490563) lc 1 lt 3 notitle, \
    'output/asymmetry_loY.txt' u ($0==0?$2:1/0):3 pt 7 ps 0.6 lc 4 title 'tanh|y_{tt}| < 0.5', \
    '' u ($0==1?$2:1/0):3 pt 6 ps 0.6 lc 4 notitle , \
    above(x,-x+0.0188165322137+0.0044505868299) lc 4 lt 1 notitle, \
    above(x,-x+0.0228061761161+0.00332675255366) lc 4 lt 3 notitle, \
    'output/asymmetry_full_points.txt' u (altox($14,off-10*ep)):(altoy($14,off-10*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+20*ep)):(altoy($14,off+20*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+50*ep)):(altoy($14,off+50*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($14,off+80*ep)):(altoy($14,off+80*ep)) w lines lt 1 lw 10 lc 2 notitle, \
    '' u (altox($11,off+80*ep)):(altoy($11,off+80*ep)) w lines lt 1 lc 7 lw lwidth notitle, \
    'output/asymmetry_hiM_points.txt' u (altox($11,off+50*ep)):(altoy($11,off+50*ep)) w lines lt 1 lc 5 lw lwidth notitle, \
    'output/asymmetry_loM_points.txt' u (altox($11,off+20*ep)):(altoy($11,off+20*ep)) w lines lt 1 lc 3 lw lwidth notitle, \
    'output/asymmetry_hiY_points.txt' u (altox($11,off-10*ep)):(altoy($11,off-10*ep)) w lines lt 1 lc 1 lw lwidth notitle, \
    'output/asymmetry_loY_points.txt' u (altox($11,off-10*ep)):(altoy($11,off-10*ep)) w lines lt 1 lc 4 lw lwidth notitle

