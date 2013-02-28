set term postscript eps enhanced 22
set output 'dqq_scan.eps'

set grid
set xlabel '{/Symbol d}_{q@^{/=18-}q}'
set ylabel '{/Symbol a}_{L}'

set arrow from 0,5 to 0,-2 head filled lc 2 lw 2 lt 1

plot [] [-2:6]\
    1 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'dqq_scan.txt' using 1:3 w lines lc 1 lt 1 lw 3 title '{/Symbol a}_{L}', \
    'dqq_scan.txt' using 1:6 w lines lc 1 lt 2 lw 2 title '{/Symbol D} PLL=0.5', \
    'dqq_scan.txt' using 1:7 w lines lc 1 lt 2 lw 2 notitle, \
    'dqq_scan.txt' using 1:8 w lines lc 3 lt 2 title '{/Symbol \261 4%} lumi', \
    'dqq_scan.txt' using 1:9 w lines lc 3 lt 2 notitle, \
    'dqq_scan.txt' using 1:10 w lines lc 4 lt 2 title '{/Symbol \261 4% s}_{DY}', \
    'dqq_scan.txt' using 1:11 w lines lc 4 lt 2 notitle, \
    'dqq_scan.txt' using 1:12 w lines lc 5 lt 2 title '{/Symbol \261 4% s}_{ST}', \
    'dqq_scan.txt' using 1:13 w lines lc 5 lt 2 notitle, \
    'dqq_scan.txt' using 1:($2<=-6792851.8620?$3:1/0) ps 1.5 pt 7 title 'best fit', \
    'dqq_scan.txt' using 1:($2+6792851.8620) w lines ls 0 title '-LL_{min}'
    #'dqq_scan.txt' using 1:($2+6792851.87) w lines lc 1 title 'NLL'
    #'dqq_scan.txt' using 1:($3+$5) w lines lc 8 title 'PLL=2', \
    #'dqq_scan.txt' using 1:($3-$5) w lines lc 8 notitle, \
    #'dqq_scan.txt' using 1:($3+$4) w lines title 'alphaL+/-', \
    #'dqq_scan.txt' using 1:($3-$4) w lines lc 2 notitle, \
