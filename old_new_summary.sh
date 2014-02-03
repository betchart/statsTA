set term postscript landscape enhanced color
set output 'graphics/old_new.ps'

set xlabel 'Reconstructed intrinsic A@_c^y   (%)'
set ylabel 'Reconstructed intrinsic A@_c^{/Symbol f}   (%)'

set grid

set size square
set key out left box

set title 'POWHEG Templates: q@^{\261}q --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttqq_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
  'ttqg_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
    'ttag_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle


set title 'POWHEG Templates: qg --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttqg_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
    'ttag_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
  'ttqq_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle

set title 'POWHEG Templates: @^{\261}qg --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttag_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
  'ttqg_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
  'ttqq_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle

ttqqelx = 1.27876312101
ttqqely = 0.340582398093
ttqqmux = 1.18611282378
ttqqmuy = 0.372053765215
        	 
ttqgelx = 0.860171410246
ttqgely = 0.702806808464
ttqgmux = 0.633168062179
ttqgmuy = 0.629287542245
        	 
        	 
ttagelx = -0.15455663305
ttagely = -0.502209020179
ttagmux = -0.155589022525
ttagmuy = -0.654336363133

set title 'Extended model asymmetries'
set object circle at ttqqelx,ttqqely size 0.02 fc rgb 'red'
set object circle at ttqqmux,ttqqmuy size 0.02 fc rgb 'red'
set object circle at ttqgelx,ttqgely size 0.02 fc rgb 'blue'
set object circle at ttqgmux,ttqgmuy size 0.02 fc rgb 'blue'
set object circle at ttagelx,ttagely size 0.02 fc rgb 'green'
set object circle at ttagmux,ttagmuy size 0.02 fc rgb 'green'

aL = 0.00114024
aT = 0.718295
set object circle at aL*ttqqelx,aL*ttqqely size 0.015 
set object circle at aL*ttqqmux,aL*ttqqmuy size 0.015 
set object circle at aT*ttqgelx,aT*ttqgely size 0.015 
set object circle at aT*ttqgmux,aT*ttqgmuy size 0.015 
set object circle at aT*ttagelx,aT*ttagely size 0.015 
set object circle at aT*ttagmux,aT*ttagmuy size 0.015 

plot [-2:2] [-2:2] \
    ttqqely/ttqqelx * x w lines lt 1 lc 1 title 'q@^{\261}q (e)', \
    ttqqmuy/ttqqmux * x w lines lt 2 lc 1 title 'q@^{\261}q ({/Symbol m})', \
    ttqgely/ttqgelx * x w lines lt 1 lc 3 title 'qg (e)', \
    ttqgmuy/ttqgmux * x w lines lt 2 lc 3 title 'qg ({/Symbol m})', \
    ttagely/ttagelx * x w lines lt 1 lc 2 title '@^{\261}qg (e)', \
    ttagmuy/ttagmux * x w lines lt 2 lc 2 title '.     @^{\261}qg ({/Symbol m})'

