set term postscript eps enhanced 22
set grid

set xlabel '{/Symbol d}_{q@^{/=18-}q}'
dqq_best=0.075
dqq_pll05=0.21
set arrow 4 from 0,5 to 0,-2 head filled lc 2 lw 2 lt 1

set ylabel 'Profile -logL({/Symbol d}_{q@^{/=18-}q})'
set arrow 1 from dqq_best,0. to dqq_best,-2 head filled 
set arrow 2 from dqq_best+dqq_pll05,0.5 to dqq_best+dqq_pll05,-2 head lt 2
set arrow 3 from dqq_best-dqq_pll05,0.5 to dqq_best-dqq_pll05,-2 head lt 2

set output 'dqq_pll.eps'
plot [] [-2:6]\
    1/0 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'dqq_scan.txt' using 1:($2+6792851.8620) w lines ls 1 lw 2 lc 3 title 'PLL({/Symbol d}_{q@^{/=18-}q})', \
    0.5 lt 2 lc 0 notitle
    

set ylabel '{/Symbol a}_{L}'

set arrow 1 from dqq_best,-1 to dqq_best,-2 head filled 
set arrow 2 from dqq_best+dqq_pll05,3 to dqq_best+dqq_pll05,-2 head lt 2
set arrow 3 from dqq_best-dqq_pll05,3 to dqq_best-dqq_pll05,-2 head lt 2

set arrow 6 from 0,1 to -0.595,1 head filled lt 1 lw 2 lc 2
set arrow 5 from -0.5,0 to -0.595,0 head filled lt 1 lw 2 lc 0
set label '{}_{no asymmetry}' at -0.585,-0.1

set style fill transparent pattern 4 bo

set output 'dqq_scan.eps'
plot [] [-2:6]\
    (x>dqq_best-dqq_pll05&&x<dqq_best+dqq_pll05)?3:1/0 w filledcurves y1=-2 lw 1 lt 2 lc rgb "#eeeeee" notitle, \
    1 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'dqq_scan.txt' using 1:3 w lines lc 1 lt 1 lw 3 title 'PLL_{/Symbol d}({/Symbol a}_L)=0.0', \
    'dqq_scan.txt' using 1:6 w lines lc 1 lt 2 lw 2 title 'PLL_{/Symbol d}({/Symbol a}_L)=0.5', \
    'dqq_scan.txt' using 1:7 w lines lc 1 lt 2 lw 2 notitle, \
    'dqq_scan.txt' using 1:8 w lines lc 3 lt 2 title '{/Symbol \261 4%} lumi', \
    'dqq_scan.txt' using 1:9 w lines lc 3 lt 2 notitle, \
    'dqq_scan.txt' using 1:10 w lines lc 4 lt 2 title '{/Symbol \261 4% s}_{DY}', \
    'dqq_scan.txt' using 1:11 w lines lc 4 lt 2 notitle, \
    'dqq_scan.txt' using 1:12 w lines lc 5 lt 2 title '{/Symbol \261 4% s}_{ST}', \
    'dqq_scan.txt' using 1:13 w lines lc 5 lt 2 notitle, \
    'dqq_scan.txt' using 1:($2+6792851.8620) w lines ls 0 notitle
    #'dqq_scan.txt' using 1:($2<=-6792851.8620?$3:1/0):(0.2) w xerrorbars ps 1.5 pt 7 title 'best fit', \
    #'dqq_scan.txt' using 1:($2+6792851.87) w lines lc 1 title 'NLL'
    #'dqq_scan.txt' using 1:($3+$5) w lines lc 8 title 'PLL=2', \
    #'dqq_scan.txt' using 1:($3-$5) w lines lc 8 notitle, \
    #'dqq_scan.txt' using 1:($3+$4) w lines title 'alphaL+/-', \
    #'dqq_scan.txt' using 1:($3-$4) w lines lc 2 notitle, \
