set term postscript eps enhanced 22
set grid

set xlabel '{/Symbol d}_{q@^{/=18-}q}'
set ylabel 'Profile -logL({/Symbol d}_{q@^{/=18-}q})'
set arrow 4 from 0,5 to 0,-2 head filled lc 2 lw 2 lt 1

set output 'graphics/dqq_pll_fitTopPtOverSumPt.eps'
set title 'tt.pt/(t.pt+tbar.pt)'
pll_min=-7144679.2660
plot [] [-2:6]\
    1/0 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'data/dqq_scan_fitTopPtOverSumPt.txt' using 1:(-pll_min+$2) w lines ls 1 lw 2 lc 3 title 'PLL({/Symbol d}_{q@^{/=18-}q})', \
    0.5 lt 2 lc 0 notitle
    

set output 'graphics/dqq_pll_fitTopTanhAvgAbsSumRapidities.eps'
set title '(|y.t|+|y.tbar|)/2'
pll_min=-7062560.4882
plot [] [-2:6]\
    1/0 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'data/dqq_scan_fitTopTanhAvgAbsSumRapidities.txt' using 1:(-pll_min+$2) w lines ls 1 lw 2 lc 3 title 'PLL({/Symbol d}_{q@^{/=18-}q})', \
    0.5 lt 2 lc 0 notitle
    

set output 'graphics/dqq_pll_fitTopTanhRapiditySum.eps'
set title '|tt.y|'
pll_min=-7062113.4747
plot [] [-2:6]\
    1/0 lc 2 lw 3 lt 1 title 'POWHEG CT10', \
    'data/dqq_scan_fitTopTanhRapiditySum.txt' using 1:(-pll_min+$2) w lines ls 1 lw 2 lc 3 title 'PLL({/Symbol d}_{q@^{/=18-}q})', \
    0.5 lt 2 lc 0 notitle
    
