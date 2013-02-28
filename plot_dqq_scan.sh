set grid
set xlabel 'd_qq'

plot [] [-2:6]\
    'dqq_scan.txt' using 1:3 w lines title 'alphaL', \
    'dqq_scan.txt' using 1:($3+$4) w lines title 'alphaL+/-', \
    'dqq_scan.txt' using 1:($3-$4) w lines lc 2 notitle, \
    'dqq_scan.txt' using 1:6 w lines lc 3 title 'PLL=0.5', \
    'dqq_scan.txt' using 1:7 w lines lc 3 notitle, \
    'dqq_scan.txt' using 1:($3+$5) w lines lc 8 title 'PLL=2', \
    'dqq_scan.txt' using 1:($3-$5) w lines lc 8 notitle, \
    'dqq_scan.txt' using 1:8 w lines lc 5 title 'lumi+/-', \
    'dqq_scan.txt' using 1:9 w lines lc 5 notitle, \
    'dqq_scan.txt' using 1:10 w lines lc 4 title 'xs_dy+', \
    'dqq_scan.txt' using 1:11 w lines lc 4 notitle, \
    'dqq_scan.txt' using 1:12 w lines lc 7 title 'xs_st+', \
    'dqq_scan.txt' using 1:13 w lines lc 7 notitle, \
    1 lc 2, \
    'dqq_scan.txt' using 1:($2+6792851.87) w lines lc 1 title 'NLL'
